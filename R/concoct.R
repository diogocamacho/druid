#' concoct: DRUID wrapper
#'
#' Start-to-finish wrapper to identify drug profiles that mimic or revert a
#' query phenotype. Computes cosine similarity to the drug corpus, an empirical
#' null, and DRUID scores.
#'
#' @param dge_matrix Nx2 matrix: column 1 log2 fold-change, column 2 p-value.
#' @param num_random Number of random sets for the empirical null. Defaults to 1000.
#' @param druid_direction \code{"pos"} or \code{"neg"}. Defaults to \code{"neg"}.
#' @param fold_thr Absolute log2FC gate. Defaults to 0.
#' @param pvalue_thr P-value gate. Defaults to 0.05.
#' @param entrez Entrez IDs aligned to \code{dge_matrix} rows.
#' @param dataset Dataset name or id (\code{cmap}, \code{lincs}, ..., \code{all}).
#'   If missing, falls back to interactive \code{\link{compendium_selection}}.
#' @param tfidf_mode Corpus weighting mode. Defaults to \code{"geom_tf"}
#'   (geometric mean of dual TF-IDF views). Also \code{"combined"},
#'   \code{"binary"}, \code{"drug_tfidf"}.
#' @param min_matches Minimum feature overlaps to report. Defaults to 3.
#' @param n_cores Cores for the empirical null (fork). Defaults to 1.
#' @return A tibble sorted by decreasing DRUID score.
#' @export
concoct <- function(dge_matrix,
                    num_random = 1000,
                    druid_direction = c("neg", "pos"),
                    fold_thr = 0,
                    pvalue_thr = 0.05,
                    entrez,
                    dataset = NULL,
                    tfidf_mode = c("geom_tf", "combined", "binary", "drug_tfidf"),
                    min_matches = 3,
                    n_cores = 1L) {

  message("Checks and balances...")
  if (missing(dge_matrix)) stop("Need differential expression data.")
  if (ncol(dge_matrix) != 2) stop("Differential expression data needs to be Nx2 matrix.")
  if (missing(entrez)) stop("Need EntrezIDs for genes in dge_matrix.")

  druid_direction <- match.arg(druid_direction)
  tfidf_mode <- match.arg(tfidf_mode)
  if (missing(num_random)) num_random <- 1000
  if (missing(fold_thr)) fold_thr <- 0
  if (missing(pvalue_thr)) pvalue_thr <- 0.05

  message("--- Concocting with DRUID ---")
  message(paste("tfidf_mode:", tfidf_mode))
  message("")

  if (is.null(dataset)) {
    message("No dataset specified; launching interactive selection.")
    druid_data <- compendium_selection()
    if (is.null(druid_data) || is.na(druid_data)) stop("No valid dataset selected.")
  } else {
    druid_data <- normalize_dataset_selection(dataset)$id
  }
  message("")

  run_one <- function(sel) {
    message(paste("Running DRUID on dataset selection", sel, "with tfidf_mode =", tfidf_mode))
    run_druid(
      dge_matrix = dge_matrix,
      druid_direction = druid_direction,
      fold_thr = fold_thr,
      pvalue_thr = pvalue_thr,
      entrez = entrez,
      num_random = num_random,
      selection = sel,
      min_matches = min_matches,
      tfidf_mode = tfidf_mode,
      n_cores = n_cores
    )
  }

  if (druid_data != 6) {
    res <- run_one(druid_data)
  } else {
    message(":: Running DRUID on all data sets ::")
    message("!!Warning: depending on processor speed, this could take a while.")
    message("")
    res <- vector(mode = "list", length = 5)
    for (i in seq_len(5)) {
      res[[i]] <- run_one(i)
      message("")
    }
    res <- dplyr::bind_rows(res)
  }

  message("DONE.")
  res
}
