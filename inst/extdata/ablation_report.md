# DRUID corpus ablation report (Stouffer + variants)

Date: 2026-09-02

## Setup

- **Corpus:** Harmonizome CMap (`cmap_druid`), 6100 × 19672
- **Queries:** 6052 leave-one-out profiles with ≥1 same-name partner
- **Query:** binary gene-direction (matched to ternary −1/0/+1)
- **Rank:** cosine; self excluded
- Legacy Hadamard `combined` vs recomputed: correlation 1.0000 on nnz

## Modes

| mode | definition |
|------|------------|
| `binary` | one-hot support |
| `drug_tfidf` | \((1/s_d)\log(D/\mathrm{df}_g)\) |
| `gene_tfidf` | \((1/\mathrm{df}_g)\log(T/s_d)\) |
| `combined` | Hadamard \(W_D \times W_G\) (historical) |
| `geom_mean` | \(\sqrt{W_D W_G}\) |
| `arith_mean` | \((W_D + W_G)/2\) |
| `stouffer_tfidf` | \(\big(\mathrm{row\text{-}z}(W_D)+\mathrm{col\text{-}z}(W_G)\big)/\sqrt{2}\) (nonzero backgrounds) |
| `stouffer_global` | \(\big(\mathrm{global\text{-}z}(W_D)+\mathrm{global\text{-}z}(W_G)\big)/\sqrt{2}\) |

## Summary (n = 6052)

| mode | median rank | mean rank | top1 % | top5 % | top25 % | median self-cos |
|------|------------:|----------:|-------:|-------:|--------:|----------------:|
| binary | 570 | 919 | 7.27 | 12.03 | **18.46** | 1.000 |
| drug_tfidf | 570 | 923 | 7.47 | 12.11 | 18.36 | 0.976 |
| geom_mean | **567** | 922 | 6.92 | 11.67 | 17.94 | 0.806 |
| arith_mean | 603 | 943 | 6.53 | 10.72 | 16.74 | 0.602 |
| gene_tfidf | 609 | 953 | 6.41 | 10.59 | 16.24 | 0.509 |
| combined | 622 | 967 | 6.13 | 9.90 | 15.52 | 0.399 |
| stouffer_tfidf | 732 | 1223 | 5.83 | 9.70 | 13.93 | **−0.046** |
| stouffer_global | 1220 | 1534 | 2.07 | 4.15 | 7.37 | **−0.098** |

**Winner by decision rule (top25 → top5 → top1 → median rank): `binary`.**

`geom_mean` has the best median rank but loses slightly on top-25 vs binary/`drug_tfidf`.

## Interpretation

1. **Hadamard `combined` still underperforms** single-view `drug_tfidf` and binary.
2. **`geom_mean` is the best dual-view NLP fuse** here — softer than product, competitive with binary.
3. **`arith_mean` sits between** geom_mean and Hadamard; better than product, worse than geom_mean/binary.
4. **`stouffer_tfidf` / `stouffer_global` hurt recovery.** Z-scoring centers edge weights so a flat binary query is no longer aligned with the drug row (median self-cosine ≈ 0 or negative). Stouffer fuse is fine for *edge ranking*; it is a poor corpus transform for *cosine against an unstandardized binary query*.

## Implication for the CLR/Stouffer premise

The premise (combine both bipartite contexts without privileging one axis) is still sound. On this retrieval setup it is better realized by **`geom_mean` or `arith_mean` of raw TF-IDF views** than by z-score + Stouffer then cosine.

A fairer test of Stouffer for connectivity would need a query transform matched to the corpus (e.g. z-score the query the same way, or use Stouffer only to filter/rank edges before building a binary/weighted query) — not cosine of binary \(q\) against Stouffer-z rows.

## Artifacts

- `ablation_results/ablation_summary.tsv`
- `ablation_results/ablation_per_query.tsv`
- `inst/scripts/run_ablation.R`
- `R/corpus_weighting.R`
