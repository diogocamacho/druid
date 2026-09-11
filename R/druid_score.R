#' DRUID Score
#'
#' Composite score from normalized cosine similarity and empirical null
#' probability.
#'
#' @param similarity_results Vector of cosine similarities (one per drug).
#' @param random_probabilities Empirical null probabilities from
#'   \code{\link{random_probability}}.
#' @param num_random Number of random draws used for the null (for p=0 floor).
#' @return Numeric vector of DRUID scores.
#' @keywords internal
druid_score <- function(similarity_results, random_probabilities, num_random) {
  similarity_results <- as.numeric(similarity_results)
  random_probabilities <- as.numeric(random_probabilities)
  n <- length(similarity_results)
  if (length(random_probabilities) != n) {
    stop("random_probabilities must have the same length as similarity_results.")
  }

  max_sim <- max(similarity_results, na.rm = TRUE)
  if (!is.finite(max_sim) || max_sim <= 0) {
    return(rep(1, n))
  }

  norm_sim <- similarity_results / max_sim
  score <- 1 + norm_sim + (-log10(random_probabilities))

  p_floor <- 1 / num_random
  zero_p <- random_probabilities <= 0 | !is.finite(random_probabilities)
  score[zero_p] <- 1 + norm_sim[zero_p] - log10(p_floor)

  score[!is.finite(score)] <- 1
  score
}
