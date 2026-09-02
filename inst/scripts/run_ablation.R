#!/usr/bin/env Rscript
# Primary corpus ablation on Harmonizome CMap (bundled cmap_druid).
# Leave-one-out same-drug-name recovery by cosine rank.
# Query TF is binary (matched to ternary Harmonizome signatures).

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
n_queries <- if (length(args) >= 1) as.integer(args[[1]]) else 400L
seed <- if (length(args) >= 2) as.integer(args[[2]]) else 42L
out_dir <- if (length(args) >= 3) args[[3]] else "ablation_results"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- normalizePath(sub("^--file=", "", file_arg))
  root <- normalizePath(file.path(dirname(script_path), "../.."))
} else {
  root <- getwd()
}

source(file.path(root, "R/corpus_weighting.R"))
source(file.path(root, "R/ablation_score.R"))

message("Loading cmap_druid...")
load(file.path(root, "data/cmap_druid.RData"))

tfidf_legacy <- cmap_druid$tfidf
drugs <- cmap_druid$drugs
drug_names <- as.character(drugs$name)

message("Building corpus modes from binary support of legacy ctfidf...")
t0 <- proc.time()[["elapsed"]]
modes <- build_corpus_modes(tfidf_legacy, from_weighted = TRUE)
message(sprintf("Corpus modes built in %.1fs", proc.time()[["elapsed"]] - t0))

# Sanity: legacy combined vs recomputed combined (allow small numeric drift)
B <- modes$binary$tfidf
recombined <- modes$combined$tfidf
# Compare on shared nonzero pattern correlation of values
common <- (B != 0)
legacy_vals <- tfidf_legacy[common]
recomb_vals <- recombined[common]
corr <- suppressWarnings(cor(as.numeric(legacy_vals), as.numeric(recomb_vals)))
message(sprintf("Correlation legacy ctfidf vs recomputed combined (nnz): %.4f", corr))

set.seed(seed)
# Prefer drugs that have at least one other profile with the same name
name_counts <- table(drug_names)
eligible <- which(as.integer(name_counts[drug_names]) >= 2L)
if (length(eligible) < n_queries) {
  warning("Fewer eligible queries than requested; using all eligible.")
  query_idx <- eligible
} else {
  query_idx <- sort(sample(eligible, n_queries))
}
message(sprintf("Running %d leave-one-out queries (seed=%d)", length(query_idx), seed))

ks <- c(1L, 5L, 25L)
rows <- list()
row_i <- 1L

for (mode_name in names(modes)) {
  message("Mode: ", mode_name)
  W <- modes[[mode_name]]$tfidf
  cpm <- modes[[mode_name]]$cpm
  t_mode <- proc.time()[["elapsed"]]

  for (qi in query_idx) {
    q <- as.numeric(B[qi, ])  # binary query from ternary support
    cs <- ablation_cosine(q, W, cpm)
    m <- ablation_recovery_metrics(qi, drug_names, cs, ks = ks)
    rows[[row_i]] <- data.frame(
      mode = mode_name,
      query_idx = m$query_idx,
      drug_name = m$drug_name,
      n_same_name_others = m$n_same_name_others,
      best_same_name_rank = m$best_same_name_rank,
      hit_top1 = m$hit_top1,
      hit_top5 = m$hit_top5,
      hit_top25 = m$hit_top25,
      self_cosine = cs[qi],
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
  message(sprintf("  done in %.1fs", proc.time()[["elapsed"]] - t_mode))
}

results <- do.call(rbind, rows)
results_path <- file.path(out_dir, "ablation_per_query.tsv")
write.table(results, results_path, sep = "\t", quote = FALSE, row.names = FALSE)

summarize_mode <- function(df) {
  ranks <- df$best_same_name_rank
  data.frame(
    mode = df$mode[1],
    n_queries = nrow(df),
    median_best_rank = median(ranks, na.rm = TRUE),
    mean_best_rank = mean(ranks, na.rm = TRUE),
    pct_hit_top1 = 100 * mean(df$hit_top1),
    pct_hit_top5 = 100 * mean(df$hit_top5),
    pct_hit_top25 = 100 * mean(df$hit_top25),
    pct_any_same_name = 100 * mean(!is.na(ranks)),
    median_self_cosine = median(df$self_cosine),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(split(results, results$mode), summarize_mode))
mode_order <- c(
  "binary", "drug_tfidf", "gene_tfidf", "combined",
  "geom_mean", "arith_mean", "stouffer_tfidf", "stouffer_global"
)
summary_df <- summary_df[match(mode_order, summary_df$mode), ]
summary_df <- summary_df[!is.na(summary_df$mode), ]
summary_path <- file.path(out_dir, "ablation_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Decision: maximize top-25 then top-5 then top-1; tie-break lower median rank
score <- summary_df$pct_hit_top25 * 1e6 + summary_df$pct_hit_top5 * 1e3 +
  summary_df$pct_hit_top1 - summary_df$median_best_rank
winner <- summary_df$mode[which.max(score)]

report <- c(
  "# DRUID corpus ablation report",
  "",
  sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Corpus: bundled Harmonizome CMap (`cmap_druid`), %d profiles x %d features", nrow(B), ncol(B)),
  sprintf("Queries: %d leave-one-out profiles with >=1 same-name partner (seed=%d)", length(query_idx), seed),
  "Query encoding: binary gene-direction (matched to ternary Harmonizome signatures)",
  "Ranking metric: cosine similarity; self profile excluded from rank",
  sprintf("Legacy vs recomputed Hadamard combined TF-IDF correlation (nnz): %.4f", corr),
  "",
  "## Modes",
  "",
  "- `binary`: one-hot support",
  "- `drug_tfidf` / `gene_tfidf`: single-view TF-IDF",
  "- `combined`: Hadamard product W_D * W_G (historical)",
  "- `geom_mean`: sqrt(W_D * W_G)",
  "- `arith_mean`: (W_D + W_G) / 2",
  "- `stouffer_tfidf`: row-z(W_D) + col-z(W_G), then /sqrt(2) (nonzero backgrounds)",
  "- `stouffer_global`: global-z(W_D) + global-z(W_G), then /sqrt(2)",
  "",
  "## Summary",
  "",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("**Recommended default `tfidf_mode`: `%s`**", winner),
  "",
  "## Notes",
  "",
  "- Asymmetric continuous query TF (|log2FC|) was not run: corpus is Harmonizome -1/0/+1.",
  "- Monte Carlo DRUID null was not used for ranking; cosine recovery is the primary ablation signal.",
  "- Stouffer z-scores use nonzero-only backgrounds (sparse-safe); structural zeros stay 0.",
  sprintf("- Per-query results: `%s`", results_path),
  sprintf("- Summary table: `%s`", summary_path),
  ""
)

report_path <- file.path(out_dir, "ablation_report.md")
writeLines(report, report_path)

message("Winner: ", winner)
message("Wrote ", summary_path)
message("Wrote ", report_path)
print(summary_df, row.names = FALSE)
