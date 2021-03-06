% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/druid_geneset.R
\name{druid_geneset}
\alias{druid_geneset}
\title{DRUID geneset}
\usage{
druid_geneset(dge_matrix, desired_effect = c("pos", "neg"), fold_thr,
  pvalue_thr, entrez, gene_space)
}
\arguments{
\item{dge_matrix}{This is a 2 column matrix for gene expression changes, where column 1 is the gene fold change and column 2 is the corresponding p-value for the fold change. NOTE: Use log2 of the fold changes as output, for example, from `limma` or `DESeq2`.}

\item{desired_effect}{Desired effect for DRUID to run on: "pos" mimics query phenotype, "neg" reverts query phenotype. Defaults to "neg".}

\item{fold_thr}{Threshold for the fold change to be considered. Defaults to 0 (i.e., log2(1), where fold change is not used as filter)}

\item{pvalue_thr}{Threshold for the p-value of the fold change to be considered. Defaults to 0.05.}

\item{entrez}{EntrezIDs for genes in differentially expressed set. Must be same order as the input matrix.}

\item{gene_space}{EntrezIDs from TF.IDF matrix (column names of this matrix).}
}
\value{
A gene set to be used as input in DRUID as a Nx4 matrix.
}
\description{
This function generates the query gene set that is to be used as input into DRUID.
}
