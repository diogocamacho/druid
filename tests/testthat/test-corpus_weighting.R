test_that("combined_tfidf matches closed form on tiny binary matrix", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  B <- Matrix::rsparsematrix(5, 8, density = 0.4)
  B@x <- rep(1, length(B@x))
  rownames(B) <- paste0("d", seq_len(nrow(B)))
  colnames(B) <- paste0("g", seq_len(ncol(B)))

  W <- DRUID:::combined_tfidf(B)
  D <- nrow(B)
  Tm <- ncol(B)
  s <- as.numeric(Matrix::rowSums(B))
  df <- as.numeric(Matrix::colSums(B))

  sm <- as.data.frame(Matrix::summary(B))
  d <- sm$i[1]
  g <- sm$j[1]
  expected <- (1 / (s[d] * df[g])) * log(D / df[g]) * log(Tm / s[d])
  expect_equal(as.numeric(W[d, g]), expected, tolerance = 1e-10)
})

test_that("geom_mean is sqrt of product for dual views", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  B <- Matrix::rsparsematrix(4, 6, density = 0.5)
  B@x <- rep(1, length(B@x))
  rownames(B) <- paste0("d", seq_len(nrow(B)))
  colnames(B) <- paste0("g", seq_len(ncol(B)))

  G <- DRUID:::geom_mean_tfidf(B)
  Wd <- DRUID:::drug_tfidf(B)
  Wg <- DRUID:::gene_tfidf(B)
  expect_equal(G@x, sqrt(Wd@x * Wg@x), tolerance = 1e-12)
})
