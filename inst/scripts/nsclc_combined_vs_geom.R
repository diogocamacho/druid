suppressPackageStartupMessages({
  library(Matrix)
  library(limma)
})
source("R/corpus_weighting.R")
source("R/ablation_score.R")
source("R/druid_geneset.R")

load("/Users/dcamacho/Downloads/camacho_druid/data/lung_cancer_gse19804.RData")
tumor_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[grep("tumor", samples_lung$clinical_diagnosis)])
healthy_lung <- which(colnames(expression_lung) %in% samples_lung$cel_file[grep("normal", samples_lung$clinical_diagnosis)])
design <- model.matrix(~ 0 + factor(c(rep("case", length(tumor_lung)), rep("ctr", length(healthy_lung)))))
colnames(design) <- c("case", "ctr")
expr <- cbind(expression_lung[, tumor_lung], expression_lung[, healthy_lung])
fit2 <- eBayes(contrasts.fit(lmFit(expr, design), makeContrasts(case - ctr, levels = design)))
lung_res <- topTable(fit2, number = Inf, sort.by = "none")
cat(sprintf("DE |logFC|>1 & adj.P<0.01: %d\n", sum(abs(lung_res$logFC) > 1 & lung_res$adj.P.Val < 0.01)))

load("data/cmap_druid.RData")
modes <- build_corpus_modes(cmap_druid$tfidf, from_weighted = TRUE)
B <- modes$binary$tfidf
drugs <- cmap_druid$drugs
qv <- druid_geneset(cbind(lung_res$logFC, lung_res$adj.P.Val), "neg", 1, 0.01,
                    as.character(genes_lung$ENTREZID), colnames(B))
cat(sprintf("Query features: %d\n", sum(qv)))
nmatch <- as.integer(Matrix::rowSums(B[, qv != 0, drop = FALSE] != 0))

score_mode <- function(mode_name) {
  cs <- as.numeric(ablation_cosine(qv, modes[[mode_name]]$tfidf, modes[[mode_name]]$cpm))
  data.frame(mode = mode_name, drug_name = as.character(drugs$name),
             concentration = drugs$concentration, cell_line = as.character(drugs$cell_line),
             number_matches = nmatch, cosine = cs, stringsAsFactors = FALSE)
}
res_c <- score_mode("combined")
res_g <- score_mode("geom_mean")
res_c <- res_c[res_c$number_matches >= 3, ]
res_g <- res_g[res_g$number_matches >= 3, ]

rankify <- function(df) {
  df <- df[order(df$cosine, decreasing = TRUE), ]
  df$rank <- seq_len(nrow(df))
  df
}
res_c <- rankify(res_c)
res_g <- rankify(res_g)
u_c <- res_c[!duplicated(res_c$drug_name), ]
u_g <- res_g[!duplicated(res_g$drug_name), ]

top_c_prof <- head(res_c, 25)
top_g_prof <- head(res_g, 25)
top_c_drug <- head(u_c, 25)
top_g_drug <- head(u_g, 25)

rg <- setNames(u_g$rank, u_g$drug_name)
rc <- setNames(u_c$rank, u_c$drug_name)
cg <- setNames(u_g$cosine, u_g$drug_name)
cc <- setNames(u_c$cosine, u_c$drug_name)

tab1 <- data.frame(
  drug_name = top_c_drug$drug_name,
  rank_combined = top_c_drug$rank,
  rank_geom = as.integer(rg[top_c_drug$drug_name]),
  delta = as.integer(rg[top_c_drug$drug_name]) - top_c_drug$rank,
  cos_combined = round(top_c_drug$cosine, 5),
  cos_geom = round(as.numeric(cg[top_c_drug$drug_name]), 5),
  cell = top_c_drug$cell_line,
  stringsAsFactors = FALSE
)
tab2 <- data.frame(
  drug_name = top_g_drug$drug_name,
  rank_geom = top_g_drug$rank,
  rank_combined = as.integer(rc[top_g_drug$drug_name]),
  delta = top_g_drug$rank - as.integer(rc[top_g_drug$drug_name]),
  cos_geom = round(top_g_drug$cosine, 5),
  cos_combined = round(as.numeric(cc[top_g_drug$drug_name]), 5),
  cell = top_g_drug$cell_line,
  stringsAsFactors = FALSE
)

overlap <- length(intersect(top_c_drug$drug_name, top_g_drug$drug_name))
only_c <- setdiff(top_c_drug$drug_name, top_g_drug$drug_name)
only_g <- setdiff(top_g_drug$drug_name, top_c_drug$drug_name)

out <- "ablation_results/nsclc_gse19804"
dir.create(out, FALSE, TRUE)
write.table(top_c_prof, file.path(out, "top25_profiles_combined_cosine.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_g_prof, file.path(out, "top25_profiles_geom_mean_cosine.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab1, file.path(out, "top25_unique_from_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab2, file.path(out, "top25_unique_from_geom.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== GSE19804 NSCLC: tumor vs normal, neg, |FC|>1, adjP<0.01; CMap; cosine rank ===\n")
cat(sprintf("Unique-drug top25 overlap: %d/25\n", overlap))
cat("Only in combined top25:", paste(only_c, collapse = ", "), "\n")
cat("Only in geom_mean top25:", paste(only_g, collapse = ", "), "\n")
cat(sprintf("Spearman unique-drug ranks: %.3f\n", cor(u_c$rank, as.integer(rg[u_c$drug_name]), method = "spearman")))
cat(sprintf("Median |delta| for combined top25 unique: %.1f\n", median(abs(tab1$delta))))

cat("\n--- Top 25 UNIQUE drugs by combined (geom rank + delta) ---\n")
print(tab1, row.names = FALSE)
cat("\n--- Top 25 UNIQUE drugs by geom_mean (combined rank + delta) ---\n")
print(tab2, row.names = FALSE)
cat("\n--- Top 25 PROFILES combined ---\n")
print(top_c_prof[, c("rank", "drug_name", "cell_line", "cosine", "number_matches")], row.names = FALSE)
cat("\n--- Top 25 PROFILES geom_mean ---\n")
print(top_g_prof[, c("rank", "drug_name", "cell_line", "cosine", "number_matches")], row.names = FALSE)
