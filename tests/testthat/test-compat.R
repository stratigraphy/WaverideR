test_that("base replacements behave", {
  expect_equal(WaverideR:::wr_closest(c(1, 5, 9), 6), 5)
  m <- matrix(c(1, 2, 3, 4, 5, 7), 2)
  expect_equal(WaverideR:::wr_row_sds(m), apply(m, 1, stats::sd))
  expect_equal(WaverideR:::wr_row_mins(m), apply(m, 1, min))
  set.seed(1)
  x <- WaverideR:::wr_rtruncnorm(200, a = 0, b = 1, mean = 0.5, sd = 0.5)
  expect_true(all(x >= 0 & x <= 1))
  expect_length(WaverideR:::wr_coloured_noise(50, 0, 1, 0.5), 50)
})
test_that("palettes work without extra packages", {
  expect_length(suppressWarnings(WaverideR:::wr_brewer_pal(5, "Spectral")), 11)
  expect_length(WaverideR:::wr_viridis(10, "inferno"), 10)
})
