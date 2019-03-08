#' DRUID Score
#'
#' This function computes the DRUID score for a given query vector. The DRUID score is dependent on the cosine similarity value as well as the random probability for any given drug against the query vector.
#'
#' @param similarity_results Vector of cosine similarities
#' @param gs_size Size of the query set
#' @param num_sets Number of random sets to generate
#' @param target_tfidf Drug tf-idf matrix.
#' @param tfidf_crossprod_mat Cross-product matrix for drug tf-idf.
#' @return A vector of random probabilities for each drug given the gene set size.
druid_score <- function(similarity_results, random_probabilities, num_random) {
  score <- 1 + (similarity_results / max(similarity_results))  + (-log10(random_probabilities))
  score[random_probabilities == 0] <- 1 + (similarity_results[random_probabilities == 0] / max(similarity_results)) - log10(1 / num_random)
  return(score)
}