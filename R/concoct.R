#' concoct: DRUID wrapper
#'
#' This function is a wrapper for a start to finish to identify drugs/drug profiles that mimic or revert a given profile of interest.
#' For a given gene set of interest, the algorithm will compute the cosine similarity between this vector and the gene sets in a given collection of drug profiles.
#' Additionally, the algorithm will generate a set of random gene sets that will be used to determine the probability of a given profile similarity being random.
#' The output of the wrapper is a data frame with all the (unranked) DRUID scores.
#'
#' @param dge_matrix This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.
#' @param tfidf_matrix tf-idf matrix drug-gene matrix. Column names are Entrez IDs. Computed with \code{\link{ctfidf}}
#' @param num_random Number of random sets to be generated to calculate significance of enrichment.  Defaults to 1,000.
#' @param druid_direction Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".
#' @param fold_thr Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)
#' @param pvalue_thr Threshold for the p-value of the fold change to be considered. Defaults to 0.05.
#' @param entrez EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.
#' @return A data frame that is sorted on the DRUID score.
concoct <- function(dge_matrix, 
                    # tfidf_matrix, 
                    num_random, 
                    druid_direction, 
                    fold_thr, 
                    pvalue_thr, 
                    entrez)
{
  # run checks ----
  message("Checks and balances...")
  if(missing(dge_matrix)) stop("Need differential expression data.")
  if(ncol(dge_matrix) != 2) stop("Differential expression data needs to be Nx2 matrix.")
  # if(missing(tfidf_matrix)) stop("Need TF-IDF matrix.")
  # if(class(tfidf_matrix) != "dgCMatrix" | class(tfidf_matrix) != "matrix") stop("TF-IDF matrix is of wrong type. Please revise.")
  # if(sum(tfidf_matrix) == 0) stop("Stopping: TF-IDF matrix is zero. Please revise.")
  # if(missing(tfidf_crossproduct)) stop("Need cross-product vector.")
  if(missing(druid_direction)) druid_direction <- "neg"
  if(missing(num_random)) num_random <- 1000
  if(missing(fold_thr)) fold_thr <- 0
  if(missing(pvalue_thr)) pvalue_thr <- 0.05
 
  message("--- Concocting with DRUID ---")
  message("")
  
  druid_data <- compendium_selection()
  message("") 

  if (druid_data != 6) {
    message(paste("Running DRUID on the", names(cauldron::druid_potion)[druid_data], "dataset"))
    res <- run_druid(dge_matrix = dge_matrix, 
                     druid_direction = druid_direction, 
                     fold_thr = fold_thr, 
                     pvalue_thr = pvalue_thr, 
                     entrez = entrez, 
                     num_random = num_random, 
                     selection = druid_data)
  } else { 
    message(":: Running DRUID on all data sets ::")
    message("!!Warning: depending on processor speed, this could take a while. Go get a coffee or something.")
    message("")
    
    res <- vector(mode = "list", length = 5)
    for (i in seq(1, 5)) {
      message(paste("Running DRUID on the", names(cauldron::druid_potion)[i], "dataset"))
      res[[i]] <- run_druid(dge_matrix = dge_matrix, 
                            druid_direction = druid_direction, 
                            fold_thr = fold_thr, 
                            pvalue_thr = pvalue_thr, 
                            entrez = entrez, 
                            num_random = num_random, 
                            selection = i)
      message("")
    }
    res <- dplyr::bind_rows(res)
  }

  message("DONE.")
  return(res)
  
}