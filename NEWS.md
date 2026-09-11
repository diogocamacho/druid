# DRUID 1.2.0

* `concoct()` / `run_druid()` take `tfidf_mode`; default `geom_tf` (geometric mean of dual TF-IDF views); also `combined`, `binary`, `drug_tfidf`
* vectorized `random_probability()` via batched sparse multiplies; optional `n_cores` fork parallelism
* headless `dataset=` on `concoct()`; CMAP falls back to bundled `cmap_druid` when `cauldron` is absent
* correctness guards in `druid_score`, `druid_geneset`, and input validation throughout
* package contract: DESCRIPTION Imports, explicit NAMESPACE exports, `.Rbuildignore`
* README + getting-started vignette; testthat + GitHub Actions R-CMD-check
* data versioning doc: `inst/extdata/DATA_SOURCES.md`
* fixed `ctfidf` docs (Hadamard dual TF-IDF, not matrix product)

# DRUID 1.1.1 (2020-09-14)

* fixed a redundancy on reporting when a signature yielded no significant genes

# DRUID 1.1.0 (2020-09-08)

* filtering out drug matches with less than `min_matches`
* added `AnnotationDbi` and `org.Hs.eg.db` dependencies
* added gene symbols to matched genes in query–drug overlap

# DRUID 1.0.0 (2020-09-08)

* re-versioning for easier following of updates
* updated DESCRIPTION for Bioconductor `biocViews`
* fixed bug where signatures with 0 matches could get good DRUID scores
