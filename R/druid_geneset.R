#' DRUID geneset
#'
#' This function generates the query gene set that is to be used as input into DRUID.
#'
#' @param dge_matrix This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.
#' @param desired_effect Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".
#' @param fold_thr Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)
#' @param pvalue_thr Threshold for the p-value of the fold change to be considered. Defaults to 0.05.
#' @param entrez EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.
#' @param gene_space EntrezIDs from TF.IDF matrix (column names of this matrix).
#' @return A gene set to be used as input in DRUID as a Nx4 matrix.
druid_geneset <- function(dge_matrix, 
                          desired_effect = c("pos", "neg"), 
                          fold_thr, 
                          pvalue_thr, 
                          entrez, 
                          gene_space)
{
  if(ncol(dge_matrix) != 2) stop("Differential expression data needs to be Nx2 matrix.")
  if(missing(desired_effect)) desired_effect <- "neg"
  if(missing(fold_thr)) fold_thr <- 0
  if(missing(pvalue_thr)) pvalue_thr <- 0.05
  if(missing(entrez)) stop("Need EntrezID of genes in data.")
  
  query_vector <- integer(length = length(gene_space))
  
  if(desired_effect == "neg")
  {
    dge_matrix[, 1] <- -1 * dge_matrix[, 1]
  }

  a2 <- which(dge_matrix[, 1] > 0)
  a3 <- which(dge_matrix[, 1] < 0)
  gs_dir <- character(length = nrow(dge_matrix))
  gs_dir[a2] <- "up"
  gs_dir[a3] <- "down"
  gs_eff <- paste(entrez, gs_dir) # <--- entrez with direction
  
  x1 <- which(abs(dge_matrix[, 1]) > fold_thr & dge_matrix[, 2] < pvalue_thr)
  x2 <- gs_eff[x1]
  
  # query vector
  query_vector[which(gene_space %in% x2)] <- 1
  qsize <- sum(query_vector)
  if(qsize == 0) stop("No genes in gene set.")
  
  return(query_vector)
}