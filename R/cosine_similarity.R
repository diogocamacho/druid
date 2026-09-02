#' Vector cosine similarity
#'
#' Cosine similarity between a query vector and each *row* of a TF-IDF matrix.
#'
#' @param query_vector Numeric vector (or 1 x G matrix) of length ncol(tfidf_matrix)
#' @param tfidf_matrix Sparse (or dense) drug x feature matrix
#' @param tfidf_crossprod_mat Row squared-norms of \code{tfidf_matrix}
#' @return Numeric vector of cosine similarities (one per row)
cosine_similarity <- function(query_vector, tfidf_matrix, tfidf_crossprod_mat) {
  if (is.matrix(query_vector)) {
    x2 <- as.numeric(tcrossprod(query_vector))
    y1 <- as.numeric(tfidf_matrix %*% t(query_vector))
  } else {
    query_vector <- as.numeric(query_vector)
    x2 <- as.numeric(crossprod(query_vector))
    y1 <- as.numeric(tfidf_matrix %*% query_vector)
  }

  x3 <- sqrt(as.numeric(tfidf_crossprod_mat) * x2)
  cs <- y1 / x3
  cs[!is.finite(cs)] <- 0
  cs
}
