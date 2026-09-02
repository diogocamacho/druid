#' Cosine similarities for one query against a corpus mode
#'
#' @param query_vector Numeric vector length ncol(tfidf)
#' @param tfidf Sparse corpus matrix
#' @param cpm Row squared-norms of tfidf
#' @return Numeric vector of cosine similarities (one per drug row)
ablation_cosine <- function(query_vector, tfidf, cpm) {
  x2 <- as.numeric(crossprod(query_vector))
  if (x2 == 0) return(rep(0, nrow(tfidf)))
  y1 <- as.numeric(tfidf %*% query_vector)
  denom <- sqrt(cpm * x2)
  cs <- y1 / denom
  cs[!is.finite(cs)] <- 0
  cs
}

#' Leave-one-out recovery metrics for one query profile
#'
#' @param query_idx Row index used to build the query
#' @param drug_names Character vector of drug names (length = nrow corpus)
#' @param cosine Numeric cosine vector against full corpus
#' @param ks Integer vector of top-k cutoffs
#' @return Named list of metrics
ablation_recovery_metrics <- function(query_idx, drug_names, cosine, ks = c(1, 5, 25)) {
  self_name <- drug_names[query_idx]
  keep <- seq_along(cosine) != query_idx
  cos_o <- cosine[keep]
  names_o <- drug_names[keep]
  ord <- order(cos_o, decreasing = TRUE)
  ranked_names <- names_o[ord]
  same <- ranked_names == self_name
  if (!any(same)) {
    best_rank <- NA_integer_
  } else {
    best_rank <- which(same)[1]
  }
  out <- list(
    query_idx = query_idx,
    drug_name = self_name,
    best_same_name_rank = best_rank,
    n_same_name_others = sum(drug_names == self_name) - 1L
  )
  for (k in ks) {
    out[[paste0("hit_top", k)]] <- as.integer(!is.na(best_rank) && best_rank <= k)
  }
  out
}
