#' Random probability (empirical null for cosine similarity)
#'
#' For each drug profile, estimates the fraction of random binary query vectors
#' (same size as the real query) whose cosine similarity exceeds the observed
#' similarity.
#'
#' Implementation notes:
#' - Random queries are evaluated in batches via one sparse matrix multiply
#'   per batch (avoids per-draw R loop + \code{rbind} of indicator matrices).
#' - Optional multicore batching via \code{parallel::mclapply} when
#'   \code{n_cores > 1} (fork backend; ignored on Windows).
#'
#' @param similarity_results Numeric vector of observed cosine similarities
#'   (one per drug / row of \code{target_tfidf}).
#' @param gs_size Size of the query set (number of active features).
#' @param num_sets Number of random sets to generate.
#' @param target_tfidf Drug x feature weight matrix.
#' @param tfidf_crossprod_mat Row squared-norms of \code{target_tfidf}.
#' @param batch_size Random queries per matrix multiply. Defaults to 1000.
#' @param n_cores Number of fork workers for batch parallelism. Defaults to 1.
#' @return Numeric vector of empirical probabilities (one per drug).
random_probability <- function(similarity_results,
                               gs_size,
                               num_sets,
                               target_tfidf,
                               tfidf_crossprod_mat,
                               batch_size = 1000L,
                               n_cores = 1L) {
  similarity_results <- as.numeric(similarity_results)
  ndrug <- nrow(target_tfidf)
  usize <- ncol(target_tfidf)
  gs_size <- as.integer(gs_size)
  num_sets <- as.integer(num_sets)
  batch_size <- as.integer(batch_size)
  n_cores <- as.integer(n_cores)

  if (length(similarity_results) != ndrug) {
    stop("similarity_results length must equal nrow(target_tfidf).")
  }
  if (num_sets <= 0) stop("num_sets must be positive.")
  if (gs_size <= 0) return(rep(1, ndrug))
  if (gs_size > usize) stop("gs_size exceeds number of features in target_tfidf.")

  denom <- sqrt(as.numeric(tfidf_crossprod_mat) * gs_size)
  denom[!is.finite(denom) | denom <= 0] <- Inf

  # batch sizes that sum to num_sets
  n_batches <- as.integer(ceiling(num_sets / batch_size))
  batch_sizes <- rep.int(batch_size, n_batches)
  batch_sizes[n_batches] <- num_sets - batch_size * (n_batches - 1L)

  eval_batch <- function(b, seed_offset = 0L) {
    # optional seed offset for reproducible parallel streams
    if (seed_offset > 0L) set.seed(seed_offset)
    rows <- vapply(
      seq_len(b),
      function(i) sample.int(usize, gs_size, replace = FALSE),
      integer(gs_size)
    )
    Q <- Matrix::sparseMatrix(
      i = as.integer(rows),
      j = rep.int(seq_len(b), times = gs_size),
      x = 1,
      dims = c(usize, b)
    )
    scores <- as.matrix(target_tfidf %*% Q)
    scores <- scores / denom
    rowSums(scores > similarity_results, na.rm = TRUE)
  }

  use_parallel <- n_cores > 1L &&
    n_batches > 1L &&
    .Platform$OS.type != "windows" &&
    requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    # independent RNG streams per batch
    base_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      sum(as.integer(.Random.seed[seq_len(min(3L, length(.Random.seed)))]) )
    } else {
      sample.int(.Machine$integer.max, 1L)
    }
    beats_list <- parallel::mclapply(
      seq_len(n_batches),
      function(k) eval_batch(batch_sizes[[k]], seed_offset = as.integer(base_seed + k * 10007L)),
      mc.cores = min(n_cores, n_batches)
    )
    beats <- Reduce(`+`, beats_list)
  } else {
    beats <- numeric(ndrug)
    for (k in seq_len(n_batches)) {
      beats <- beats + eval_batch(batch_sizes[[k]])
    }
  }

  as.numeric(beats) / num_sets
}
