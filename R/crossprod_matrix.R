#' Cross-product matrix for tf-idf
#'
#' Returns the cross-product for an input tf-idf matrix
#'
#' @param tfidf_matrix Sparse tf-idf matrix
#' @return cross-product of tf-idf matrix
#' @export
crossprod_matrix <- function(tfidf_matrix) {
  as.numeric(Matrix::rowSums(tfidf_matrix * tfidf_matrix))
}