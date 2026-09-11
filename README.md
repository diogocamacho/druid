# DRUID: DRUg Indication Discoverer

DRUID ranks drug transcriptional profiles that **mimic** or **revert** a query gene expression signature. It uses a Harmonizome-style drug corpus (gene-direction tokens), dual TF-IDF weighting (default: geometric mean of drug- and gene-as-document views), cosine similarity, and an empirical null model.

## Installation

```r
# install.packages("devtools")
devtools::install_github("diogocamacho/druid")
```

Optional: install [`cauldron`](https://github.com/diogocamacho/cauldron) for LINCS, CTD, and other compendia beyond bundled CMAP.

## Quick start

```r
library(DRUID)

# dge_matrix: column 1 = log2 fold-change, column 2 = p-value
# entrez: Entrez IDs in the same order as dge_matrix rows

res <- concoct(
  dge_matrix = dge_matrix,
  entrez = entrez_ids,
  dataset = "cmap",           # headless; omit for interactive picker
  tfidf_mode = "geom_tf",     # recommended (default)
  druid_direction = "neg",    # revert disease signature
  fold_thr = 1,
  pvalue_thr = 0.01,
  num_random = 1000,
  n_cores = 1                 # increase for faster null (fork, non-Windows)
)
```

Results are a tibble sorted by `druid_score`, with drug metadata, overlap counts, cosine similarity, and empirical p-values.

### Corpus weighting modes (`tfidf_mode`)

| Mode | Description |
|------|-------------|
| `geom_tf` | **Default.** Geometric mean of drug-view and gene-view TF-IDF |
| `combined` | Hadamard product (historical DRUID) |
| `binary` | One-hot support only |
| `drug_tfidf` | Drug-as-document TF-IDF only |

### Datasets (`dataset`)

| Value | Source |
|-------|--------|
| `cmap` | Bundled Harmonizome CMAP (no `cauldron` required) |
| `lincs`, `small_molecules`, `natural_products`, `ctd` | Requires `cauldron` |
| `all` | Run all five compendia and bind results |

## Bring your own corpus

Build a binary drug × gene-direction matrix, then compute weights:

```r
library(Matrix)
# data_matrix: sparse binary, rows = drug profiles, cols = "ENTREZ up/down"

W <- ctfidf(data_matrix)              # combined (Hadamard)
cpm <- crossprod_matrix(W)

res <- run_druid(
  dge_matrix = dge_matrix,
  entrez = entrez_ids,
  tfidf_matrix = W,
  tfidf_crossproduct = cpm,
  drugs = drug_metadata_df,
  data_source = "my_corpus"
)
```

Or use `prepare_druid_corpus()` internally via `selection` + `tfidf_mode` on bundled/cauldron data.

## NSCLC example (GSE19804)

See `inst/scripts/nsclc_combined_vs_geom.R` and `ablation_results/nsclc_gse19804/` for a manuscript-style comparison of `combined` vs `geom_tf` on lung tumor vs normal.

Typical manuscript settings:

```r
res <- concoct(
  dge_matrix = cbind(lung_res$logFC, lung_res$adj.P.Val),
  entrez = genes_lung$ENTREZID,
  dataset = "cmap",
  tfidf_mode = "geom_tf",
  druid_direction = "neg",
  fold_thr = 1,
  pvalue_thr = 0.01,
  num_random = 10000
)
```

## Data notes

Harmonizome CMap/LINCS signatures are **ternary** (−1 / 0 / +1) encoded as discrete gene-direction features. Query signatures should use **binary direction** after FC/p-value gating — not continuous fold-change weights. See `inst/extdata/DATA_SOURCES.md`.

## Ablation results

Corpus mode comparison (6052 leave-one-out CMAP queries): summaries in `inst/extdata/ablation_summary.tsv` and `ablation_results/ablation_report.md`.

## Citation

Camacho et al., DRUID — Drug Indication Discoverer.
