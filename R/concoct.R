#' concoct: DRUID wrapper
#'
#' This function is a wrapper for a start to finish to identify drugs/frug profiles that mimic or revert a given profile of interest.
#' For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of pathways (for example, KEGG pathways.)
#' Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given enrichment being random.
#' The output of the wrapper is a tidy data frame with all the enrichment scores that are better than 0.
#'
#' @param dge_matrix This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.
#' @param tfidf_matrix tf-idf matrix drug-gene matrix. Column names are Entrez IDs. Computed with \code{\link{ctfidf}}
#' @param tfidf_crossproduct Cross-product of the tf-idf matrix. Computed with \code{\link{crossprod_matrix}}.
#' @param num_random Number of random sets to be generated to calculate significance of enrichment.  Defaults to 1,000.
#' @param druid_direction Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".
#' @param fold_thr Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)
#' @param pvalue_thr Threshold for the p-value of the fold change to be considered. Defaults to 0.05.
#' @param entrez EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.
#' @return A data frame.
concoct <- function(dge_matrix, tfidf_matrix, tfidf_crossproduct, num_random, druid_direction, fold_thr, pvalue_thr, entrez)
{
  # message("\014")
  # message("+----------------------------------------------+")
  # message("| DRUID: Drug Indication Discoverer            |")
  # message("|                                              |")
  # message("| Diogo M. Camacho, Ph.D.                      |")
  # message("| Wyss Institute @ Harvard University          |")
  # message("+----------------------------------------------+")
  # message("")
  message("--- Concocting with DRUID ---")
  
  # run checks ----
  message("Checks and balances...")
  if(missing(dge_matrix)) stop("Need differential expression data.")
  if(ncol(dge_matrix) != 2) stop("Differential expression data needs to be Nx2 matrix.")
  if(missing(tfidf_matrix)) stop("Need TF-IDF matrix.")
  if(missing(tfidf_crossproduct)) stop("Need cross-product vector.")
  if(missing(druid_direction)) druid_direction <- "neg"
  if(missing(num_random)) num_random <- 1000
  if(missing(fold_thr)) fold_thr <- 0
  if(missing(pvalue_thr)) pvalue_thr <- 0.05
  
  # generate query vector ----
  message("Generating query vector...")
  query_vector <- druid_geneset(dge_matrix = dge_matrix, 
                                desired_effect = druid_direction, 
                                fold_thr = fold_thr, 
                                pvalue_thr = pvalue_thr, 
                                entrez = entrez, 
                                gene_space = colnames(tfidf_matrix))

  # compute cosine similarities against query vector ----
  message("Computing cosine similarity on query vector...")
  query_similarities <- cosine_similarity(tfidf_matrix = tfidf_matrix,
                                          query_vector = query_vector,
                                          tfidf_crossprod_mat = tfidf_crossproduct)
  
  # random probabilities ----
  message("Computing cosine similarity on random vectors...")
  prandom <- random_probability(similarity_results = query_similarities,
                                gs_size = sum(query_vector),
                                num_sets = num_random,
                                target_tfidf = tfidf_matrix,
                                tfidf_crossprod_mat = tfidf_crossproduct)
  
  # DRUID Score ----
  message("Computing cosine similarity on random vectors...")
  dscore <- druid_score(similarity_results = query_similarities, random_probabilities = prandom, num_random = num_random)
  
  # results ----
  message("Building results data frame...")
  res <- data_frame(cosine_similarity = as.vector(query_similarities),
                    probability_random = prandom,
                    druid_score = as.vector(dscore))
  
  # count number of matches to query vector
  tt <- colnames(tfidf_matrix)[which(query_vector != 0)]
  t2 <- apply(tfidf_matrix, 1, function(y) length(intersect(tt, names(which(y != 0)))))
  
  res <- res %>% 
    tibble::add_column(., number_matches = t2, .before = 1)
  
  message("DONE.")
  return(res)
  
}