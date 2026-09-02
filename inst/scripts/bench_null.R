#!/usr/bin/env Rscript
# Benchmark: legacy loop null vs batched sparse null; full concoct timing.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
})

root <- if (file.exists("R/random_probability.R")) getwd() else normalizePath(".")
setwd(root)

for (f in c(
  "corpus_weighting.R", "prepare_corpus.R", "crossprod_matrix.R",
  "cosine_similarity.R", "druid_geneset.R", "random_probability.R",
  "druid_score.R", "run_druid.R", "concoct.R", "compendium_selection.R"
)) {
  source(file.path("R", f))
}

load("data/cmap_druid.RData")

# ---- legacy null (for accuracy + timing comparison) ----
random_probability_legacy <- function(similarity_results, gs_size, num_sets, target_tfidf, tfidf_crossprod_mat) {
  usize <- ncol(target_tfidf)
  rnd_sim <- vector(mode = "list", length = num_sets)
  for (i in seq_len(num_sets)) {
    qv <- integer(usize)
    wm <- integer(nrow(target_tfidf))
    qv[sample(usize, size = gs_size, replace = FALSE)] <- 1L
    a1 <- crossprod(qv)
    a2 <- sqrt(tfidf_crossprod_mat * as.vector(a1))
    a3 <- target_tfidf %*% qv
    cs <- a3 / a2
    wm[which(as.vector(cs > similarity_results))] <- 1L
    rnd_sim[[i]] <- wm
  }
  colSums(do.call(rbind, rnd_sim)) / num_sets
}

corpus <- prepare_druid_corpus("cmap", "geom_tf")
W <- corpus$tfidf
cpm <- corpus$cpm
B <- corpus$binary

# mid-size query from a real drug profile
feat_idx <- which(as.numeric(B[100, ]) != 0)
gs_size <- min(200L, length(feat_idx))
qv <- numeric(ncol(B))
qv[feat_idx[seq_len(gs_size)]] <- 1
cs <- cosine_similarity(qv, W, cpm)

num_sets <- 1000L
set.seed(1)
t_new <- system.time({
  p_new <- random_probability(cs, gs_size, num_sets, W, cpm, batch_size = 500L)
})[["elapsed"]]

set.seed(1)
t_old <- system.time({
  p_old <- random_probability_legacy(cs, gs_size, num_sets, W, cpm)
})[["elapsed"]]

# Same seed => identical samples if we reshuffle the same way... actually
# legacy samples inside loop and new samples in batches — different stream use.
# Check statistical agreement: mean abs diff should be small (Monte Carlo noise)
# For identity check, run new twice with same seed:
set.seed(42)
p_a <- random_probability(cs, gs_size, num_sets, W, cpm)
set.seed(42)
p_b <- random_probability(cs, gs_size, num_sets, W, cpm)

cat(sprintf("Null gs_size=%d num_sets=%d drugs=%d features=%d\n", gs_size, num_sets, nrow(W), ncol(W)))
cat(sprintf("legacy elapsed: %.2fs\n", t_old))
cat(sprintf("batched elapsed: %.2fs\n", t_new))
cat(sprintf("speedup: %.1fx\n", t_old / max(t_new, 1e-6)))
cat(sprintf("deterministic replay max|diff|: %.3g\n", max(abs(p_a - p_b))))
cat(sprintf("legacy vs batched mean|diff|: %.4f  max|diff|: %.4f\n", mean(abs(p_old - p_new)), max(abs(p_old - p_new))))

# Full concoct (includes corpus rebuild + null)
set.seed(1)
entrez <- sub(" (up|down)$", "", colnames(B)[feat_idx[seq_len(gs_size)]])
dir_up <- grepl(" up$", colnames(B)[feat_idx[seq_len(gs_size)]])
dge <- cbind(ifelse(dir_up, 1, -1), 0.001)
t_full <- system.time({
  res <- concoct(
    dge_matrix = dge, num_random = num_sets, druid_direction = "pos",
    fold_thr = 0, pvalue_thr = 0.05, entrez = entrez,
    dataset = "cmap", tfidf_mode = "geom_tf", min_matches = 3
  )
})[["elapsed"]]
cat(sprintf("full concoct (geom_tf, num_random=%d): %.2fs  nrow=%d\n", num_sets, t_full, nrow(res)))
