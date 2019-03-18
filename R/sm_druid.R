sm_druid <- function(dge_matrix, druid_direction, fold_thr, pvalue_thr, entrez, num_random) {
  # generate query vector ----
  message("Generating query vector...")
  query_vector <- druid_geneset(dge_matrix = dge_matrix, 
                                desired_effect = druid_direction, 
                                fold_thr = fold_thr, 
                                pvalue_thr = pvalue_thr, 
                                entrez = entrez, 
                                gene_space = colnames(cauldron::druid_potion$small_molecules$tfidf))
  
  # compute cosine similarities against query vector ----
  message("Computing cosine similarity on query vector...")
  query_similarities <- cosine_similarity(tfidf_matrix = cauldron::druid_potion$small_molecules$tfidf,
                                          tfidf_crossprod_mat = cauldron::druid_potion$small_molecules$cpm,
                                          query_vector = query_vector)
  
  # random probabilities ----
  message("Computing cosine similarity on random vectors...")
  prandom <- random_probability(similarity_results = query_similarities,
                                gs_size = sum(query_vector),
                                num_sets = num_random,
                                target_tfidf = cauldron::druid_potion$small_molecules$tfidf,
                                tfidf_crossprod_mat = cauldron::druid_potion$small_molecules$cpm)
  
  # DRUID Score ----
  message("Computing cosine similarity on random vectors...")
  dscore <- druid_score(similarity_results = query_similarities, 
                        random_probabilities = prandom, 
                        num_random = num_random)
  
  # results ----
  message("Building results data frame...")
  res <- tibble(cosine_similarity = as.vector(query_similarities),
                probability_random = prandom,
                druid_score = as.vector(dscore))
  
  # count number of matches to query vector
  tt <- colnames(cauldron::druid_potion$small_molecules$tfidf)[which(query_vector != 0)]
  t2 <- apply(cauldron::druid_potion$small_molecules$tfidf, 1, function(y) length(intersect(tt, names(which(y != 0)))))
  
  res <- res %>% 
    tibble::add_column(., number_matches = t2, .before = 1)
  
  res <- res %>% 
    tibble::add_column(., drug_name = cauldron::druid_potion$small_molecules$drugs, .before = 1) %>%
    tibble::add_column(., concentration = NA, .before = 2) %>%
    tibble::add_column(., cell_line = NA, .before = 3) %>%
    tibble::add_column(., data_source = "small_molecules", .before = 1) %>%
    dplyr::arrange(., desc(druid_score))
  
  return(res)
}