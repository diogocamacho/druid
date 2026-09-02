#' Normalize DRUID dataset selection to an integer id and name
#'
#' @param selection Integer 1-6, or character name (cmap, lincs, ...), or "all"
#' @return list(id=, name=)
normalize_dataset_selection <- function(selection) {
  name_map <- c(
    "cmap" = 1L,
    "lincs" = 2L,
    "small_molecules" = 3L,
    "small-molecules" = 3L,
    "natural_products" = 4L,
    "natural-products" = 4L,
    "ctd" = 5L,
    "all" = 6L
  )
  if (is.character(selection)) {
    key <- tolower(gsub("\\s+", "_", selection))
    if (!key %in% names(name_map)) {
      stop("Unknown dataset '", selection, "'. Use cmap, lincs, small_molecules, natural_products, ctd, or all.")
    }
    id <- unname(name_map[[key]])
  } else {
    id <- as.integer(selection)
    if (length(id) != 1 || is.na(id) || id < 1 || id > 6) {
      stop("Dataset selection must be an integer in 1..6 or a known dataset name.")
    }
  }
  labels <- c("cmap", "lincs", "small_molecules", "natural_products", "ctd", "all")
  list(id = id, name = labels[[id]])
}

#' Load one drug-compendium slice (tfidf, cpm, drugs)
#'
#' Uses \code{cauldron::druid_potion} when available. For CMAP only, falls back
#' to the bundled \code{cmap_druid} data object.
#'
#' @param selection Integer 1-5 or dataset name (not "all")
#' @return list(name=, tfidf=, cpm=, drugs=)
get_compendium_slice <- function(selection) {
  sel <- normalize_dataset_selection(selection)
  if (sel$id == 6L) stop("get_compendium_slice() does not accept 'all'; pass a single dataset.")

  if (requireNamespace("cauldron", quietly = TRUE)) {
    potion <- cauldron::druid_potion
    nm <- names(potion)[sel$id]
    return(list(
      name = nm,
      tfidf = potion[[sel$id]]$tfidf,
      cpm = potion[[sel$id]]$cpm,
      drugs = potion[[sel$id]]$drugs
    ))
  }

  if (sel$id == 1L) {
    dat <- NULL
    if (exists("cmap_druid", inherits = TRUE)) {
      dat <- get("cmap_druid", inherits = TRUE)
    } else {
      # package data / development checkout
      pkg_data <- system.file("data", "cmap_druid.RData", package = "DRUID")
      local_data <- file.path("data", "cmap_druid.RData")
      path <- if (nzchar(pkg_data)) pkg_data else if (file.exists(local_data)) local_data else ""
      if (!nzchar(path)) stop("Bundled cmap_druid data not found.")
      e <- new.env(parent = emptyenv())
      load(path, envir = e)
      dat <- e$cmap_druid
    }
    return(list(name = "cmap", tfidf = dat$tfidf, cpm = dat$cpm, drugs = dat$drugs))
  }

  stop(
    "Package 'cauldron' is required for dataset '", sel$name,
    "'. Install cauldron, or use dataset = 'cmap' with bundled data."
  )
}

#' Prepare corpus matrices for a DRUID run
#'
#' Recovers binary support from the stored (usually combined) TF-IDF matrix and
#' rebuilds weights for the requested mode.
#'
#' @param selection Dataset id or name
#' @param tfidf_mode One of \code{geom_tf}, \code{combined}, \code{binary}, \code{drug_tfidf}
#' @return list(name=, tfidf=, cpm=, drugs=, binary=, tfidf_mode=)
prepare_druid_corpus <- function(selection,
                                 tfidf_mode = c("geom_tf", "combined", "binary", "drug_tfidf")) {
  tfidf_mode <- match.arg(tfidf_mode)
  slice <- get_compendium_slice(selection)
  B <- binary_from_weighted(slice$tfidf)

  W <- switch(
    tfidf_mode,
    geom_tf = geom_mean_tfidf(B),
    combined = combined_tfidf(B),
    binary = B,
    drug_tfidf = drug_tfidf(B)
  )
  dimnames(W) <- dimnames(B)
  cpm <- as.numeric(Matrix::rowSums(W * W))

  list(
    name = slice$name,
    tfidf = W,
    cpm = cpm,
    drugs = slice$drugs,
    binary = B,
    tfidf_mode = tfidf_mode
  )
}
