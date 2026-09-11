#' Recover binary drug x gene-direction matrix from a weighted sparse matrix
#'
#' Harmonizome-style corpora are ternary (-1/0/+1) encoded as one-hot
#' gene-direction columns. Nonzero pattern of a TF-IDF matrix recovers that
#' binary support.
#'
#' @param weighted_matrix Sparse matrix (e.g. existing ctfidf)
#' @return Binary dgCMatrix with the same dimnames
#' @keywords internal
binary_from_weighted <- function(weighted_matrix) {
  B <- weighted_matrix
  B@x <- rep(1, length(B@x))
  B
}

#' Drug-as-document TF-IDF
#'
#' For binary input: tfidf[d,g] = (1/s_d) * log(D / df_g)
#'
#' @param data_matrix Binary (or count) sparse drug x feature matrix
#' @return Sparse TF-IDF matrix
#' @keywords internal
drug_tfidf <- function(data_matrix) {
  D <- nrow(data_matrix)
  s <- Matrix::rowSums(data_matrix)
  df <- Matrix::colSums(data_matrix)
  if (any(s == 0)) stop("Empty drug signatures present (rowSums == 0).")
  df_safe <- as.numeric(df)
  df_safe[df_safe == 0] <- NA_real_
  idf <- log(D / df_safe)
  idf[is.na(idf)] <- 0
  out <- Matrix::Diagonal(x = 1 / as.numeric(s)) %*% data_matrix %*% Matrix::Diagonal(x = idf)
  dimnames(out) <- dimnames(data_matrix)
  out
}

#' Gene-as-document TF-IDF, aligned to drug x gene layout
#'
#' For binary input: tfidf[d,g] = (1/df_g) * log(T / s_d)
#'
#' @param data_matrix Binary sparse drug x feature matrix
#' @return Sparse TF-IDF matrix (same dims as data_matrix)
#' @keywords internal
gene_tfidf <- function(data_matrix) {
  Tm <- ncol(data_matrix)
  s <- as.numeric(Matrix::rowSums(data_matrix))
  df <- as.numeric(Matrix::colSums(data_matrix))
  if (any(s == 0)) stop("Empty drug signatures present (rowSums == 0).")
  df_safe <- df
  df_safe[df_safe == 0] <- NA_real_
  w_row <- log(Tm / s)
  w_col <- 1 / df_safe
  w_col[is.na(w_col)] <- 0
  out <- Matrix::Diagonal(x = w_row) %*% data_matrix %*% Matrix::Diagonal(x = w_col)
  dimnames(out) <- dimnames(data_matrix)
  out
}

#' Combined (bidirectional) TF-IDF via Hadamard product of dual views
#'
#' For binary input:
#' ctfidf[d,g] = (1/(s_d * df_g)) * log(D/df_g) * log(T/s_d)
#'
#' This is element-wise, not a matrix product of tfidf(dtm) and t(tfidf(dtm)).
#'
#' @param data_matrix Binary (or count) sparse drug x feature matrix
#' @return Sparse combined TF-IDF matrix
#' @keywords internal
combined_tfidf <- function(data_matrix) {
  drug_tfidf(data_matrix) * gene_tfidf(data_matrix)
}

#' Geometric mean of the two TF-IDF views
#'
#' S = sqrt(W_D * W_G); softer than Hadamard product, still multiplicative.
#'
#' @param data_matrix Binary sparse drug x feature matrix
#' @return Sparse matrix
#' @keywords internal
geom_mean_tfidf <- function(data_matrix) {
  Wd <- drug_tfidf(data_matrix)
  Wg <- gene_tfidf(data_matrix)
  # element-wise sqrt of product on shared support
  out <- Wd
  out@x <- sqrt(Wd@x * Wg@x)
  out
}

#' Arithmetic mean of the two raw TF-IDF views
#'
#' S = (W_D + W_G) / 2
#'
#' @param data_matrix Binary sparse drug x feature matrix
#' @return Sparse matrix
#' @keywords internal
arith_mean_tfidf <- function(data_matrix) {
  0.5 * (drug_tfidf(data_matrix) + gene_tfidf(data_matrix))
}

#' Row-contextual z-score among observed (nonzero) weights
#'
#' For each drug row, mu/sd are computed over that row's nonzero TF-IDF values.
#' Structural zeros stay zero (no association). This avoids densifying while
#' still ranking how strong an edge is relative to other edges of the same drug.
#'
#' @param W Sparse numeric matrix
#' @return Sparse matrix of row-contextual z-scores on original support
#' @keywords internal
row_zscore_sparse <- function(W) {
  sm <- Matrix::summary(W)
  # per-row mean and sd of nonzero entries
  rs <- as.numeric(Matrix::rowSums(W))
  rss <- as.numeric(Matrix::rowSums(W * W))
  nnz <- as.numeric(Matrix::rowSums(W != 0))
  mu <- ifelse(nnz > 0, rs / nnz, 0)
  var <- ifelse(nnz > 1, pmax(rss / nnz - mu * mu, 0), 0)
  # sample-ish: use population over nnz; if nnz==1, sigma=1 so z=0
  sigma <- sqrt(var)
  sigma[nnz <= 1 | sigma == 0] <- 1
  z <- (sm$x - mu[sm$i]) / sigma[sm$i]
  z[nnz[sm$i] <= 1] <- 0
  Matrix::sparseMatrix(
    i = sm$i, j = sm$j, x = z,
    dims = dim(W), dimnames = dimnames(W),
    index1 = TRUE
  )
}

#' Column-contextual z-score among observed (nonzero) weights
#'
#' For each gene-direction column, mu/sd over nonzero drug weights.
#'
#' @param W Sparse numeric matrix
#' @return Sparse matrix of column-contextual z-scores on original support
#' @keywords internal
col_zscore_sparse <- function(W) {
  sm <- Matrix::summary(W)
  cs <- as.numeric(Matrix::colSums(W))
  css <- as.numeric(Matrix::colSums(W * W))
  nnz <- as.numeric(Matrix::colSums(W != 0))
  mu <- ifelse(nnz > 0, cs / nnz, 0)
  var <- ifelse(nnz > 1, pmax(css / nnz - mu * mu, 0), 0)
  sigma <- sqrt(var)
  sigma[nnz <= 1 | sigma == 0] <- 1
  z <- (sm$x - mu[sm$j]) / sigma[sm$j]
  z[nnz[sm$j] <= 1] <- 0
  Matrix::sparseMatrix(
    i = sm$i, j = sm$j, x = z,
    dims = dim(W), dimnames = dimnames(W),
    index1 = TRUE
  )
}

#' Global z-score over all nonzero TF-IDF weights in the matrix
#'
#' @param W Sparse numeric matrix
#' @return Sparse matrix with globally z-scored nonzero values
#' @keywords internal
global_zscore_sparse <- function(W) {
  x <- W@x
  mu <- mean(x)
  sigma <- stats::sd(x)
  if (is.na(sigma) || sigma == 0) sigma <- 1
  out <- W
  out@x <- (x - mu) / sigma
  out
}

#' Stouffer fuse of dual TF-IDF views with row/column contextual z-scores
#'
#' z_D = row-zscore(drug TF-IDF), z_G = col-zscore(gene TF-IDF),
#' S = (z_D + z_G) / sqrt(2). NLP association strengths; Stouffer combine.
#'
#' @param data_matrix Binary sparse drug x feature matrix
#' @return Sparse Stouffer-combined score matrix
#' @keywords internal
stouffer_tfidf <- function(data_matrix) {
  z_d <- row_zscore_sparse(drug_tfidf(data_matrix))
  z_g <- col_zscore_sparse(gene_tfidf(data_matrix))
  (z_d + z_g) / sqrt(2)
}

#' Stouffer fuse after global z-score of each TF-IDF view
#'
#' Softer / scale-matched variant: one global mu,sigma per view, then Stouffer.
#'
#' @param data_matrix Binary sparse drug x feature matrix
#' @return Sparse matrix
#' @keywords internal
stouffer_global_tfidf <- function(data_matrix) {
  z_d <- global_zscore_sparse(drug_tfidf(data_matrix))
  z_g <- global_zscore_sparse(gene_tfidf(data_matrix))
  (z_d + z_g) / sqrt(2)
}

#' Build a named list of corpus weightings for ablation
#'
#' @param data_matrix Binary sparse matrix (or weighted matrix whose support is used)
#' @param from_weighted If TRUE, treat data_matrix as weighted and binarize first
#' @return Named list of list(tfidf=, cpm=) per mode
#' @keywords internal
build_corpus_modes <- function(data_matrix, from_weighted = FALSE) {
  B <- if (from_weighted) binary_from_weighted(data_matrix) else data_matrix
  modes <- list(
    binary = B,
    drug_tfidf = drug_tfidf(B),
    gene_tfidf = gene_tfidf(B),
    combined = combined_tfidf(B),
    geom_mean = geom_mean_tfidf(B),
    arith_mean = arith_mean_tfidf(B),
    stouffer_tfidf = stouffer_tfidf(B),
    stouffer_global = stouffer_global_tfidf(B)
  )
  lapply(modes, function(W) {
    list(
      tfidf = W,
      cpm = as.numeric(Matrix::rowSums(W * W))
    )
  })
}
