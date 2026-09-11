#' Run DRUID
#'
#' Given a selection on the drug compendium, produces a data frame with the
#' results of DRUID.
#'
#' @param dge_matrix Nx2 matrix: column 1 log2 fold-change, column 2 p-value
#'   (e.g. from \code{limma} or \code{DESeq2}).
#' @param druid_direction \code{"pos"} mimics query phenotype, \code{"neg"}
#'   reverts it. Defaults to \code{"neg"}.
#' @param fold_thr Absolute log2FC threshold. Defaults to 0.
#' @param pvalue_thr P-value threshold. Defaults to 0.05.
#' @param entrez Entrez IDs for genes in \code{dge_matrix} (same order).
#' @param num_random Number of random gene sets for empirical null. Defaults to 1000.
#' @param selection Dataset id (1-5) or name (\code{cmap}, \code{lincs}, ...).
#'   Ignored if \code{tfidf_matrix} is supplied.
#' @param tfidf_matrix Optional drug x feature matrix (BYO corpus). When set,
#'   \code{tfidf_crossproduct} and \code{drugs} must also be supplied.
#' @param tfidf_crossproduct Row squared-norms of \code{tfidf_matrix}.
#' @param drugs Optional data frame with at least \code{name}, \code{concentration},
#'   \code{cell_line} columns (one row per drug profile).
#' @param data_source Label for results when using BYO matrices. Defaults to
#'   \code{"custom"}.
#' @param min_matches Minimum overlapping features to keep a drug profile. Defaults to 3.
#' @param tfidf_mode Corpus weighting: \code{"geom_tf"} (geometric mean of dual
#'   TF-IDF views; recommended), \code{"combined"} (Hadamard / historical),
#'   \code{"binary"}, or \code{"drug_tfidf"}.
#' @param n_cores Cores for batched null (fork). Defaults to 1 (reproducible).
#' @return A tibble sorted by decreasing DRUID score.
#' @importFrom dplyr arrange desc filter
#' @export
run_druid <- function(dge_matrix,
                      druid_direction = c("neg", "pos"),
                      fold_thr = 0,
                      pvalue_thr = 0.05,
                      entrez,
                      num_random = 1000,
                      selection = NULL,
                      tfidf_matrix = NULL,
                      tfidf_crossproduct = NULL,
                      drugs = NULL,
                      data_source = "custom",
                      min_matches = 3,
                      tfidf_mode = c("geom_tf", "combined", "binary", "drug_tfidf"),
                      n_cores = 1L) {

  if (missing(dge_matrix)) stop("Need differential expression data.")
  if (missing(entrez)) stop("Need EntrezIDs for genes in dge_matrix.")
  if (ncol(dge_matrix) != 2) stop("dge_matrix must be Nx2 (log2FC, p-value).")
  if (length(entrez) != nrow(dge_matrix)) {
    stop("entrez length must match nrow(dge_matrix).")
  }

  druid_direction <- match.arg(druid_direction)
  tfidf_mode <- match.arg(tfidf_mode)
  if (missing(min_matches)) min_matches <- 3
  if (missing(num_random)) num_random <- 1000
  if (missing(fold_thr)) fold_thr <- 0
  if (missing(pvalue_thr)) pvalue_thr <- 0.05

  if (!is.null(tfidf_matrix)) {
    if (is.null(tfidf_crossproduct)) stop("tfidf_crossproduct required with tfidf_matrix.")
    if (is.null(drugs)) stop("drugs data frame required with tfidf_matrix.")
    tfidf <- tfidf_matrix
    cpm <- as.numeric(tfidf_crossproduct)
    B <- binary_from_weighted(tfidf)
    source_name <- data_source
  } else {
    if (missing(selection) || is.null(selection)) {
      stop("Need dataset selection (e.g. 'cmap') or supply tfidf_matrix.")
    }
    corpus <- prepare_druid_corpus(selection = selection, tfidf_mode = tfidf_mode)
    tfidf <- corpus$tfidf
    cpm <- corpus$cpm
    drugs <- corpus$drugs
    B <- corpus$binary
    source_name <- corpus$name
    tfidf_mode <- corpus$tfidf_mode
  }
  gene_space <- colnames(tfidf)

  query_vector <- druid_geneset(
    dge_matrix = dge_matrix,
    desired_effect = druid_direction,
    fold_thr = fold_thr,
    pvalue_thr = pvalue_thr,
    entrez = entrez,
    gene_space = gene_space
  )

  if (sum(query_vector) != 0) {
    tt <- gene_space[which(query_vector != 0)]
    B_q <- B[, query_vector != 0, drop = FALSE]
    t2 <- as.integer(Matrix::rowSums(B_q != 0))

    # labels for query features (symbols if org.Hs.eg.db available)
    parts <- strsplit(tt, " ", fixed = TRUE)
    entrez_q <- vapply(parts, function(z) z[[1]], character(1))
    dir_q <- vapply(parts, function(z) if (length(z) >= 2) z[[2]] else NA_character_, character(1))
    if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      sym_q <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(entrez_q),
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      )
      feature_labels <- paste(unname(sym_q[entrez_q]), dir_q)
    } else {
      feature_labels <- tt
    }

    # matched gene strings via sparse summary (no per-row apply over full matrix)
    b1 <- character(nrow(B_q))
    sm <- Matrix::summary(B_q)
    if (nrow(sm) > 0) {
      labs <- feature_labels[sm$j]
      split_labs <- split(labs, sm$i)
      b1[as.integer(names(split_labs))] <- vapply(
        split_labs,
        function(z) paste(z, collapse = " | "),
        character(1)
      )
    }

    query_similarities <- cosine_similarity(
      tfidf_matrix = tfidf,
      tfidf_crossprod_mat = cpm,
      query_vector = query_vector
    )

    prandom <- random_probability(
      similarity_results = query_similarities,
      gs_size = sum(query_vector),
      num_sets = num_random,
      target_tfidf = tfidf,
      tfidf_crossprod_mat = cpm,
      n_cores = n_cores
    )
    prandom[which(t2 < min_matches)] <- 1

    dscore <- druid_score(
      similarity_results = query_similarities,
      random_probabilities = prandom,
      num_random = num_random
    )
    dscore[!is.finite(as.vector(dscore))] <- 1

    res <- tibble::tibble(
      data_source = source_name,
      tfidf_mode = tfidf_mode,
      drug_name = as.character(drugs$name),
      concentration = drugs$concentration,
      cell_line = as.character(drugs$cell_line),
      query_size = sum(query_vector),
      number_matches = t2,
      matched_genes = b1,
      cosine_similarity = as.vector(query_similarities),
      probability_random = prandom,
      druid_score = as.vector(dscore)
    )
    res <- dplyr::arrange(res, dplyr::desc(.data$druid_score))
    res <- dplyr::filter(res, .data$number_matches >= min_matches)
  } else {
    res <- tibble::tibble(
      data_source = source_name,
      tfidf_mode = tfidf_mode,
      drug_name = NA_character_,
      concentration = NA_real_,
      cell_line = NA_character_,
      query_size = 0,
      number_matches = 0,
      matched_genes = NA_character_,
      cosine_similarity = 0,
      probability_random = 1,
      druid_score = 1
    )
  }

  res
}
