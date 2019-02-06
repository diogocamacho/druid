#' Cross-product matrix for tf-idf
#'
#' Returns the cross-product for an input tf-idf matrix
#'
#' @param tfidf_matrix Sparse tf-idf matrix
#' @return cross-product of tf-idf matrix
crossprod_matrix <- function(tfidf_matrix) {
  
  cpm <- apply(tfidf_matrix, 1, crossprod) # <-- crossproduct matrix for tfidf
  return(cpm)
  
}