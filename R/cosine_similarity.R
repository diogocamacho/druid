#' Vector cosine similarity
#'
#' Returns the cosine similarity between 2 vectors 
#'
#' @param tfidf_matrix Sparse tf-idf matrix (computed from \code{\link{ctfidf}})
#' @param query_vector A vector to compare to the tf-idf matrix. Length of vector is equal to number of columns of tf-idf matrix
#' @param tfidf_crossprod_mat The cross-product matrix for the tf-idf matrix (computed with \code{\link{crossprod_matrix}})
#' @return Cosine similarity between query vector and each *row* in tf-idf matrix
cosine_similarity <- function(query_vector, tfidf_matrix, tfidf_crossprod_mat) {
  
  if(class(query_vector) == "matrix") {
    x2 <- tcrossprod(query_vector) # <-- query_vector is a 1xG vector  
    y1 <- tfidf_matrix %*% t(query_vector)
  } else if(class(query_vector) == "numeric") {
    x2 <- crossprod(query_vector) # <-- query_vector is a 1xG vector
    y1 <- tfidf_matrix %*% query_vector
  }
  
  x3 <- sqrt(tfidf_crossprod_mat * as.vector(x2))
  cs <- y1 / x3 # <-- cosine similarity
  
  return(cs)
  
}



