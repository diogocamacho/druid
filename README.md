# DRUID: DRUg Indication Discoverer

DRUID, or DRUg Indication Discoverer, is an algorithm that identifies drug profiles that revert or mimic a condition of interest.  For example, DRUID can be used to repurpose compounds for novel indications given a gene expression profile, it can be used to prioritize compounds for a compendium of disease states, and it can be incorporated into computational chemistry pipelines for the identification and characterization of drug properties for novel design. 

## Installation

Install from GitHub using `devtools` as:

```
devtools::install_github("diogocamacho/druid")
```

## Running
The easiest way to run DRUID is to use its wrapper `concoct` as:

```
library(DRUID)
res <- concoct(dge_matrix, tfidf_matrix, crossproduct_vector, number_random, effect_direction, fold_thr, pvalue_thr, entrez_ids)
```

where `dge_matrix` is a 2-column matrix for the query gene expression signature (column 1: fold-changes; column 2: p-values); `tfidf_matrix` is the calculated corrected TF-IDF matrix (see `ctfidf` function); `crossproduct_vector` is the crossproduct of the TF-IDF matrix (see `crossprod_matrix` function); `number_random` is the number of random simulations to be run to assess significance of scores (defaults to 1,000 - see `random_probability` function); `effect_direction` is the desired effect of the drug ("neg" for a reversal of the query phenotype, "pos" for a mimicking of the phenotype. Defaults to "neg"); `fold_thr` is a threshold for the fold changes (defaults to log2 = 0); `pvalue_thr` is the threshold for expression change significance (defaults to 0.05); and `entrez_ids` are the EntrezIDs for the genes in the query signature.

DRUID comes pre-packaged with a TF-IDF matrix and cross-product vector that were derived from the Connectivity Map data (as available in the [Harmonizome](http://amp.pharm.mssm.edu/Harmonizome/)). 

## Example
As an example, I will use the CMAP TF-IDF to generate a query vector and run DRUID on it.

```
gset <- unique(gsub(" down", "", gsub(" up", "", sample(colnames(DRUID::cmap_druid$tfidf), 100))))

query_matrix <- matrix(1, ncol = 2, nrow = length(gset))
query_matrix[, 2] <- 0
query_matrix[sample(x = seq(1, length(gset)), size = 0.25 * length(gset)), 1] <- -1
```

With the generated query matrix, we can now run DRUID on it using the `concoct` wrapper:

```
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
```

The output of DRUID is a `tibble` data frame with all the scores for all the drugs.  Specifically, the columns in this data frame are:

  * number_matches: number of genes in query signature found in drug signature;
  * cosine_similarity: similarity of query signature to drug signature;
  * probability_random: probability of cosine similarity being better than random;
  * druid_score: calculated DRUID score, taking into account the cosine similarity and the random probability.
  
We can expand this data frame using `magrittr` and `tibble` together with the information on the drugs as:

```
example_druid <- example_druid %>% 
  tibble::add_column(., drug_name = DRUID::cmap_druid$drugs$name, .before = 1) %>%
  tibble::add_column(., concentration = DRUID::cmap_druid$drugs$concentration, .before = 2) %>%
  tibble::add_column(., cell_line = DRUID::cmap_druid$drugs$cell_line, .before = 3)
```

We can now use `ggplot2` to visualize the results:

```
example_druid %>% dplyr::filter(., cosine_similarity == 0) %>% ggplot() + geom_point(aes(x = drug_name, y = druid_score, color = cell_line), alpha = 0.5) + facet_grid(. ~ cell_line, scales = "free") + theme_bw() + theme(axis.text.x = element_blank())
```

## Using other sources
DRUID comes pre-packaged with the TF-IDF for the Connectivity Map data, but it's simple to generate a TF-IDF to meet your needs.  For that, the `ctfidf` function will be used, which uses the `tidytext` package by Julia Silge. Inputs to this function are a data matrix where the columns are the words (eg, genes and their direction) and the rows are the documents (eg, drugs).  

A drug profile will generate a different response on the transcriptome, with genes being differentially expressed.  As such, we can generate a one-hot encoded vector in which all possibilities of change are represented for each gene in a given condition.  These are the vectors that will be present in the data matrix that serves as input to `ctfidf`, where each "word" represents a gene (as an Entrez ID) and the direction of change that the "document" (drug) caused.

(NOTE: a one-hot encoding function was not included in DRUID.  Will be included in future releases.)

With a one-hot encoded matrix, we can generate a TF-IDF matrix as:

```
ex_tfidf <- ctfidf(data_matrix)
```

and the corresponding cross-product vector as:

```
ex_cp <- crossprod_matrix(ex_tfidf)
```
