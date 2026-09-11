#' Connectivity Map data for DRUID
#'
#' Built from CMAP 2.0: each row is a drug profile; columns are gene-direction
#' features (Entrez ID plus up/down suffix).
#'
#' @format A list with three components:
#' \describe{
#'   \item{tfidf}{Sparse combined TF-IDF matrix (6100 x 19672)}
#'   \item{cpm}{Numeric vector of row squared-norms of \code{tfidf} (length 6100)}
#'   \item{drugs}{Data frame of drug metadata: \code{id}, \code{name},
#'     \code{concentration}, \code{duration}, \code{cell_line}, \code{vehicle},
#'     \code{vendor}}
#' }
#' @docType data
#' @name cmap_druid
#' @aliases cmap_druid
#' @keywords datasets
NULL
