#' Combined tf-idf matrix
#'
#' Returns a tf-idf sparse matrix given a condition-gene matrix.
#' This matrix is a DxT matrix, where D is the set of conditions (disease indications) T are the genes affecting that condition.
#' To correct effect/prevalence of genes in a given condition or to correct conditions affected by many genes, we compute the tf-idf matrix twice, to get the combined tf-idf: 
#' 
#' \eqn{ctfidf = tfidf(dtm) \times t(tfidf(dtm))}
#' 
#' @param data_matrix Sparse document-term matrix (condition by gene matrix), named (rows and columns)
#' @return A sparse matrix with the computed combined tf-idf
ctfidf <- function(data_matrix)
{
  # tidy data ----
  dtm1 <- tidytext::tidy(data_matrix)
  dtm2 <- tidytext::tidy(t(data_matrix))
  
  # count words in both dtms ----
  words1 <- dtm1 %>% 
    dplyr::count(row, column) %>% 
    dplyr::ungroup()

  words2 <- dtm2 %>% 
    dplyr::count(row, column) %>% 
    dplyr::ungroup()
  
  # count all terms ----
  terms1 <- words1 %>% 
    dplyr::group_by(row) %>% 
    dplyr::summarize(total = sum(n))

  terms2 <- words2 %>% 
    dplyr::group_by(row) %>% 
    dplyr::summarize(total = sum(n))
  
  # calculate tf-idf ----
  tfidf1 <- words1 %>% 
    tidytext::bind_tf_idf(column, row, n)

  tfidf2 <- words2 %>% 
    tidytext::bind_tf_idf(column, row, n)
  
  # combine into 1 tfidf matrix ----
  # output matrix is tfidf1 * t(tfidf2)
  x1 <- tfidf1 %>% dplyr::arrange(., row)
  x2 <- tfidf2 %>% dplyr::arrange(., column)
  x3 <- x1$tf_idf * x2$tf_idf
  
  ctfidf <- matrix(0, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  for (i in seq(1, nrow(data_matrix))) {
    b1 <- which(x1$row == rownames(data_matrix)[i])
    b2 <- match(x1$column[b1], colnames(data_matrix))
    ctfidf[i, b2] <- x3[b1]
  }
  rownames(ctfidf) <- rownames(data_matrix)
  colnames(ctfidf) <- colnames(data_matrix)
  
  # make sparse matrix ----
  ctfidf <- Matrix::Matrix(ctfidf, sparse = TRUE)

  return(ctfidf)

}