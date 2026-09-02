# NEWS

Version 1.1.1 (09/14/2020)
+ fixed a redundancy on reporting when a signature yielded no significant genes.

Version 1.1.0 (09/08/2020)
+ filtering out drug matches with less than `min_matches`
+ added `AnnotationDbi` dependency
+ added `org.Hs.eg.db` dependency
+ added gene symbols to genes that are matched from the query to a given drug signature

Version 1.0.0 (09/08/2020)
+ re-versioning number for easier following of updates;
+ updated DESCRIPTION file to include the `biocViews:` call-out for Bioconductor installs;
+ fixed a bug on score function where signatures with 0 matches where getting good DRUID scores. 

Version 1.2.0 (dev)
+ vectorized `random_probability()` via batched sparse multiplies; optional `n_cores` fork parallelism
+ typically ~2–4x vs legacy per-draw loop+rbind; ~2x more with `n_cores=4` on CMap-sized corpora
+ faster matched-gene strings in `run_druid()` (sparse summary, no full-matrix apply)
+ `concoct()` / `run_druid()` take `tfidf_mode`; default `geom_tf` (geometric mean of dual TF-IDF views)
+ also `combined` (historical Hadamard), `binary`, `drug_tfidf`
+ headless `dataset=` on `concoct()` (cmap/lincs/.../all); interactive picker only if omitted
+ CMAP falls back to bundled `cmap_druid` when `cauldron` is not installed
+ sparse corpus builders in `R/corpus_weighting.R`; ablation + NSCLC scripts under `inst/scripts/`
+ fixed ctfidf docs/math (Hadamard dual TF-IDF, not matrix product)
+ `compendium_selection()` uses `prompt=` and `stop()` on invalid input
