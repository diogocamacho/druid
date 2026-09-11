#' Combined tf-idf matrix
#'
#' Returns a combined (bidirectional) tf-idf sparse matrix for a condition-gene
#' matrix. For each nonzero (drug, gene-direction) cell this is the Hadamard
#' product of drug-as-document TF-IDF and gene-as-document TF-IDF:
#'
#' \deqn{ctfidf_{d,g} = tfidf^{(D)}_{d,g} \times tfidf^{(G)}_{g,d}}
#'
#' For binary one-hot input that equals the closed form
#' (1/(s_d df_g)) * log(D/df_g) * log(T/s_d).
#' This is not the matrix product tfidf(dtm) %*% t(tfidf(dtm)).
#'
#' @param data_matrix Sparse document-term matrix (condition by gene matrix), named
#' @return A sparse matrix with the computed combined tf-idf
#' @export
ctfidf <- function(data_matrix) {
  combined_tfidf(data_matrix)
}
