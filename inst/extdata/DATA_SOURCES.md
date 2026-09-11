# DRUID data sources and versioning

## Bundled CMAP (`cmap_druid`)

- **Object:** `DRUID::cmap_druid` (LazyData)
- **File:** `data/cmap_druid.RData`
- **Origin:** Connectivity Map 2.0 via [Harmonizome](https://maayanlab.cloud/Harmonizome/)
- **Shape:** 6,100 drug profiles × 19,672 gene-direction features (`"ENTREZ up"` / `"ENTREZ down"`)
- **Stored weights:** Historical Hadamard `combined` TF-IDF (recomputable from binary support)

At runtime, DRUID **rebuilds** corpus weights from the binary support of the stored matrix using `tfidf_mode` (default `geom_tf`). The bundled Hadamard weights are not used directly unless `tfidf_mode = "combined"` after recomputation.

## Optional compendia (`cauldron`)

Install `cauldron` for additional datasets:

| ID | Name | Notes |
|----|------|-------|
| 1 | CMAP | Same family as bundled data |
| 2 | LINCS | L1000-style Harmonizome signatures |
| 3 | Small molecules | |
| 4 | Natural products | |
| 5 | CTD | Comparative Toxicogenomics Database |

```r
# Remotes::install_github("diogocamacho/cauldron")
res <- concoct(..., dataset = "lincs")
```

## Representation constraints

- Drug side: **discrete** gene-direction tokens (Harmonizome −1/0/+1 → one-hot)
- Query side: **binary** direction after FC/p thresholds (matched representation)
- Continuous query TF (|log2FC| weights) is asymmetric with this corpus and not recommended as default

## Versioning policy

1. **Code** (this package): version in `DESCRIPTION`; corpus math in `R/corpus_weighting.R`
2. **Bundled data:** bump minor version when `cmap_druid.RData` changes; document in `NEWS.md`
3. **External (`cauldron`):** pin via `Remotes:` SHA or release tag for reproducibility

## Rebuilding corpora

From a binary drug × feature matrix `B`:

```r
source(system.file("..", "R/corpus_weighting.R", package = "DRUID")) # dev only
# or use exported helpers after install via internal API

W_geom <- geom_mean_tfidf(B)   # recommended
W_comb <- combined_tfidf(B)    # historical
cpm    <- crossprod_matrix(W_geom)
```

Ablation scripts: `inst/scripts/run_ablation.R`, `inst/scripts/nsclc_combined_vs_geom.R`
