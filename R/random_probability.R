#' Random probability
#'
#' This function computes the cosine similarity between the set of randomized gene sets and the pathway tf-idf matrix. It calls the cosine_similarity function and returns the probability of a given cosine similarity being random.
#'
#' @param similarity_results Vector of cosine similarities
#' @param gs_size Size of the query set
#' @param num_sets Number of random sets to generate
#' @param target_tfidf Pathway tf-idf matrix.
#' @param tfidf_crossprod_mat Cross-product matrix for pathway tf-idf.
#' @return A vector of random probabilities for each pathway given the gene set size.
random_probability <- function(similarity_results, gs_size, num_sets, target_tfidf, tfidf_crossprod_mat) {
  usize <- ncol(target_tfidf)
  rnd_sim <- vector(mode = "list", length = num_sets)
  for (i in seq(1, num_sets)) {
    qv <- integer(length = usize)
    wm <- integer(length = nrow(target_tfidf))
    qv[sample(usize, size = gs_size, replace = FALSE)] <- 1
    a1 <- crossprod(qv) # <-- query_vector is a 1xG vector
    a2 <- sqrt(tfidf_crossprod_mat * as.vector(a1))
    a3 <- target_tfidf %*% qv
    cs <- a3 / a2 # <-- cosine similarity
    wm[which(as.vector(cs > similarity_results))] <- 1  # <--- which ones are *better* than the real cosine similarity
    rnd_sim[[i]] <- wm
  }
  rnd_sim <- do.call(rbind, rnd_sim)
  rnd_sim <- colSums(rnd_sim) / num_sets # <--- probability that it IS random
  return(rnd_sim)
}
