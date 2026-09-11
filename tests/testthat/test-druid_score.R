test_that("druid_score handles zero max similarity", {
  sim <- c(0, 0, 0)
  pr <- c(0.1, 0.2, 0.3)
  out <- DRUID:::druid_score(sim, pr, num_random = 1000)
  expect_true(all(out == 1))
  expect_true(all(is.finite(out)))
})

test_that("druid_score handles p=0 floor", {
  sim <- c(0.5, 0.25)
  pr <- c(0, 0.01)
  out <- DRUID:::druid_score(sim, pr, num_random = 1000)
  expect_true(all(is.finite(out)))
  expect_true(out[1] > out[2])
})
