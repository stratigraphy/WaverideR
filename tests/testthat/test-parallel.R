test_that("wr_lapply is deterministic and independent of cores", {
  f <- function(i) stats::rnorm(5)
  a <- WaverideR:::wr_lapply(1:6, f, n_cores = 1, seed = 1, combine = "cbind")
  b <- WaverideR:::wr_lapply(1:6, f, n_cores = 2, run_multicore = TRUE, seed = 1, combine = "cbind")
  expect_equal(a, b)
  expect_equal(dim(a), c(5, 6))
})
test_that("wr_lapply restores the master RNG state", {
  set.seed(10); x1 <- stats::runif(1)
  set.seed(10)
  WaverideR:::wr_lapply(1:3, function(i) stats::rnorm(1), n_cores = 1, seed = 5)
  expect_equal(stats::runif(1), x1)
})
test_that("worker errors are reported", {
  expect_error(WaverideR:::wr_lapply(1:2, function(i) stop("boom"), n_cores = 1), "boom")
})
