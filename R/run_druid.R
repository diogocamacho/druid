#' Run DRUID
#' 
#' Given a selection on the drug compendium, produces a data frame with the results of DRUID. Makes the `concoct` function somewhat redundant (to be revised later.)
#'
#' @param dge_matrix This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.
#' @param tfidf_matrix tf-idf matrix drug-gene matrix. Column names are Entrez IDs. Computed with \code{\link{ctfidf}}
#' @param num_random Number of random sets to be generated to calculate significance of enrichment.  Defaults to 1,000.
#' @param druid_direction Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".
#' @param fold_thr Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)
#' @param pvalue_thr Threshold for the p-value of the fold change to be considered. Defaults to 0.05.
#' @param entrez EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.
#' @param selection EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.
#' @param min_matches Minimal number of matches needed to compute a DRUID score. Defaults to 3.
#' @return A data frame that is sorted on the DRUID score.
run_druid <- function(dge_matrix, druid_direction, fold_thr, pvalue_thr, entrez, num_random, selection, min_matches) {
  
  if(missing(min_matches)) min_matches <- 3
  
  # generate query vector ----
  # message("Generating query vector...")
  query_vector <- druid_geneset(dge_matrix = dge_matrix, 
                                desired_effect = druid_direction, 
                                fold_thr = fold_thr, 
                                pvalue_thr = pvalue_thr, 
                                entrez = entrez, 
                                gene_space = colnames(cauldron::druid_potion[[selection]]$tfidf))
  
  if (sum(query_vector) != 0) {
    # count number of matches on query vector to drug profile
    tt <- colnames(cauldron::druid_potion[[selection]]$tfidf)[which(query_vector != 0)]
    t2 <- apply(cauldron::druid_potion[[selection]]$tfidf, 1, function(y) length(intersect(tt, names(which(y != 0)))))
    
    # get symbol mappings
    a1 <- colnames(cauldron::druid_potion[[selection]]$tfidf)
    a2 <- sapply(sapply(a1, strsplit, " "), "[", 2)
    a3 <- sapply(sapply(a1, strsplit, " "), "[", 1)
    a4 <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = a3, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
    a5 <- tibble::tibble(entrez_direction = a1, symbol_direction = paste(a4, a2))
    b1 <- apply(cauldron::druid_potion[[selection]]$tfidf, 1, 
                function(y) { 
                  z1 <- which(a5$entrez_direction %in% intersect(tt, names(which(y != 0))))
                  z2 <- paste(a5$symbol_direction[z1], collapse = " | ") 
                  return(z2)})
    
    
    # compute cosine similarities against query vector ----
    # message("Computing cosine similarity on query vector...")
    query_similarities <- cosine_similarity(tfidf_matrix = cauldron::druid_potion[[selection]]$tfidf,
                                            tfidf_crossprod_mat = cauldron::druid_potion[[selection]]$cpm,
                                            query_vector = query_vector)
    
    # random probabilities ----
    # message("Computing cosine similarity on random vectors...")
    prandom <- random_probability(similarity_results = query_similarities,
                                  gs_size = sum(query_vector),
                                  num_sets = num_random,
                                  target_tfidf = cauldron::druid_potion[[selection]]$tfidf,
                                  tfidf_crossprod_mat = cauldron::druid_potion[[selection]]$cpm)
    
    prandom[which(t2 < min_matches)] <- 1
    
    # DRUID Score ----
    # message("Computing cosine similarity on random vectors...")
    dscore <- druid_score(similarity_results = query_similarities, 
                          random_probabilities = prandom, 
                          num_random = num_random)
    
    # results ----
    # message("Building results data frame...")
    res <- tibble(cosine_similarity = as.vector(query_similarities),
                  probability_random = prandom,
                  druid_score = as.vector(dscore))
    
    res <- res %>% 
      tibble::add_column(., number_matches = t2, .before = 1) %>%
      tibble::add_column(., matched_genes = b1, .before = 2)
    
    res <- res %>% 
      tibble::add_column(., drug_name = as.character(cauldron::druid_potion[[selection]]$drugs$name), .before = 1) %>%
      tibble::add_column(., concentration = cauldron::druid_potion[[selection]]$drugs$concentration, .before = 2) %>%
      tibble::add_column(., cell_line = as.character(cauldron::druid_potion[[selection]]$drugs$cell_line), .before = 3) %>%
      tibble::add_column(., data_source = names(cauldron::druid_potion)[selection], .before = 1) %>%
      dplyr::arrange(., desc(druid_score)) %>%
      dplyr::filter(., number_matches >= min_matches)
                
  } else {
    
    res <- tibble::tibble(cosine_similarity = 0,
                  probability_random = 1,
                  druid_score = 1)
    
    res <- res %>% 
      tibble::add_column(., number_matches = 0, .before = 1) %>%
      tibble::add_column(., matched_genes = NA, .before = 2)
    
    res <- res %>% 
      tibble::add_column(., drug_name = NA, .before = 1) %>%
      tibble::add_column(., concentration = NA, .before = 2) %>%
      tibble::add_column(., cell_line = NA, .before = 3) %>%
      tibble::add_column(., data_source = names(cauldron::druid_potion)[selection], .before = 1) %>%
      dplyr::arrange(., desc(druid_score)) %>%
      dplyr::filter(., number_matches >= min_matches)
  }
  
  return(res)
}
