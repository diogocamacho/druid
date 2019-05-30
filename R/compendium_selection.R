#' Compendium selection
#'
#' This function will make the selection of which data source to use to run DRUID against. All the datasets were generated with the companion package `cauldron`.
#' 
#' @return An integer.
compendium_selection <- function() {
  message("Available data sets:")
  message("[1] CMAP")
  message("[2] LINCS")
  message("[3] Small molecules")
  message("[4] Natural products")
  message("[5] Comparative Toxicogenomic Database")
  message("[6] Run all")
  message("")
  selection <- readline(promp = "Selection: ")
  selection <- as.integer(selection)
  if(selection %in% seq(1, 6)) {
    return(selection)
  } else {
    message("Not an available data set. Stopping.")
  }
}