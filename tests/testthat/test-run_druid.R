test_that("run_druid works on bundled cmap with geom_tf", {
  skip_if_not_installed("Matrix")
  cmap_druid <- NULL
  if (exists("cmap_druid", envir = asNamespace("DRUID"), inherits = FALSE)) {
    cmap_druid <- get("cmap_druid", envir = asNamespace("DRUID"))
  } else {
    path <- system.file("data", "cmap_druid.RData", package = "DRUID")
    if (!nzchar(path)) path <- file.path("data", "cmap_druid.RData")
    if (file.exists(path)) {
      e <- new.env(parent = emptyenv())
      load(path, envir = e)
      cmap_druid <- e$cmap_druid
    }
  }
  skip_if(is.null(cmap_druid), "cmap_druid data unavailable")

  B <- DRUID:::binary_from_weighted(cmap_druid$tfidf)
  feat <- colnames(B)[which(as.numeric(B[10, ]) != 0)][1:20]
  entrez <- sub(" (up|down)$", "", feat)
  dir_up <- grepl(" up$", feat)
  dge <- cbind(ifelse(dir_up, 1, -1), 0.001)

  res <- run_druid(
    dge_matrix = dge,
    entrez = entrez,
    selection = "cmap",
    tfidf_mode = "geom_tf",
    num_random = 20,
    n_cores = 1,
    min_matches = 3
  )

  expect_s3_class(res, "tbl")
  expect_true(nrow(res) >= 1)
  expect_equal(unique(res$tfidf_mode), "geom_tf")
  expect_true(all(is.finite(res$druid_score)))
})

test_that("run_druid validates inputs", {
  expect_error(
    run_druid(matrix(1, 2, 1), entrez = c("1", "2"), selection = "cmap"),
    "Nx2"
  )
  expect_error(
    run_druid(cbind(1, 0.01), entrez = "1", selection = "bogus"),
    "Unknown dataset"
  )
})
