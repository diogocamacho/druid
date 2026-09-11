#' DRUID geneset
#'
#' Build a binary query vector over the corpus gene-direction feature space.
#'
#' @param dge_matrix Nx2 matrix: column 1 log2 fold-change, column 2 p-value.
#' @param desired_effect \code{"pos"} or \code{"neg"} (reversal negates FC first).
#' @param fold_thr Absolute log2FC threshold. Defaults to 0.
#' @param pvalue_thr P-value threshold. Defaults to 0.05.
#' @param entrez Entrez IDs aligned to \code{dge_matrix} rows.
#' @param gene_space Column names of the corpus matrix (gene-direction tokens).
#' @return Integer vector of query weights (0/1) over \code{gene_space}.
#' @keywords internal
druid_geneset <- function(dge_matrix,
                          desired_effect = c("pos", "neg"),
                          fold_thr = 0,
                          pvalue_thr = 0.05,
                          entrez,
                          gene_space) {
  if (missing(dge_matrix) || ncol(dge_matrix) != 2) {
    stop("dge_matrix must be an Nx2 matrix (log2FC, p-value).")
  }
  if (missing(entrez)) stop("Need Entrez IDs for genes in dge_matrix.")
  if (missing(gene_space)) stop("Need gene_space (corpus column names).")
  if (length(entrez) != nrow(dge_matrix)) {
    stop("entrez length must match nrow(dge_matrix).")
  }

  desired_effect <- match.arg(desired_effect)
  if (missing(fold_thr)) fold_thr <- 0
  if (missing(pvalue_thr)) pvalue_thr <- 0.05

  query_vector <- integer(length(gene_space))

  if (desired_effect == "neg") {
    dge_matrix[, 1] <- -1 * dge_matrix[, 1]
  }

  up_genes <- which(dge_matrix[, 1] > 0)
  down_genes <- which(dge_matrix[, 1] < 0)
  gs_dir <- character(nrow(dge_matrix))
  gs_dir[up_genes] <- "up"
  gs_dir[down_genes] <- "down"
  gs_eff <- paste(entrez, gs_dir)

  b1 <- which(abs(dge_matrix[, 1]) > fold_thr)
  b2 <- which(dge_matrix[, 2] < pvalue_thr)

  if (length(b1) != 0 && length(b2) != 0) {
    x1 <- intersect(b1, b2)
    x2 <- gs_eff[x1]
    query_vector[which(gene_space %in% x2)] <- 1
  } else if (length(b1) != 0 && length(b2) == 0) {
    pvalue_thr <- 0.05
    b2 <- which(dge_matrix[, 2] < pvalue_thr)
    x1 <- intersect(b1, b2)
    if (length(x1) == 0) {
      message("No differentially expressed genes.")
    } else {
      message("q-value threshold yielded no statistically significant genes. Changing threshold to 0.05.")
      x2 <- gs_eff[x1]
      query_vector[which(gene_space %in% x2)] <- 1
    }
  } else if (length(b1) == 0 && length(b2) != 0) {
    fold_thr <- 0
    b1 <- which(abs(dge_matrix[, 1]) > fold_thr)
    x1 <- intersect(b1, b2)
    if (length(x1) == 0) {
      message("No differentially expressed genes.")
    } else {
      message("Fold-change threshold yielded no statistically significant genes. Changing threshold to 0.")
      x2 <- gs_eff[x1]
      query_vector[which(gene_space %in% x2)] <- 1
    }
  } else {
    message("No differentially expressed genes.")
  }

  if (sum(query_vector) == 0) {
    message("No genes matched in drug profiles.")
  }

  query_vector
}
