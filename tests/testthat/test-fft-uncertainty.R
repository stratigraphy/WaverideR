test_that("fft frequency uncertainty reproduces the worked examples", {
  a <- fft_frequency_uncertainty(T_length = 2, f = 100)
  expect_equal(a$df_practical, 0.5)
  expect_equal(a$df_gabor, 0.5 / (4 * pi), tolerance = 1e-8)
  expect_equal(a$rel_practical, 0.005)
  b <- fft_frequency_uncertainty(T_length = 0.5, f = 1000)
  expect_equal(c(b$df_practical, round(b$df_gabor, 3), b$rel_practical), c(2, 0.159, 0.002))
  expect_equal(fft_frequency_uncertainty(N = 1000, fs = 10)$df_practical, 0.01)
  expect_equal(a$ratio, 4 * pi)
  expect_error(fft_frequency_uncertainty())
})

test_that("fft_uncertainty gives ordered time models", {
  tr <- data.frame(depth = seq(0, 100, by = 1), period = 10)
  u <- fft_uncertainty(tr, period_of_tracked_cycle = 100)
  g <- fft_uncertainty(tr, period_of_tracked_cycle = 100, type = "gabor")
  expect_equal(ncol(u), 8)
  expect_true(all(u[, "time min x_times"] <= u[, "time mean"] + 1e-8))
  expect_true(all(u[, "time plus x_times"] >= u[, "time mean"] - 1e-8))
  w_p <- u[nrow(u), "time plus x_times"] - u[nrow(u), "time min x_times"]
  w_g <- g[nrow(g), "time plus x_times"] - g[nrow(g), "time min x_times"]
  expect_gt(w_p, w_g)
})
