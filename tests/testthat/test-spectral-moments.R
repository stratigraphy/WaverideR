test_that("spectral_moments keeps the window spectra and plot_spectral_moments draws them", {
  d <- mag[mag[, 1] <= 20, 1:2]
  sm <- spectral_moments(d, window = 4, step = 0.5)
  expect_true(!is.null(sm$spectra))
  expect_equal(ncol(sm$spectra$power), nrow(sm$result))
  expect_equal(length(sm$spectra$frequency), nrow(sm$spectra$power))
  expect_true(all(is.finite(sm$result[, "mean_frequency"])))
  f <- tempfile(fileext = ".png"); png(f)
  expect_silent(plot_spectral_moments(sm))
  expect_silent(plot_spectral_moments(sm, yaxis = "period", scale = "quantile"))
  expect_silent(plot_spectral_moments(sm, image = FALSE))
  dev.off()
})

test_that("the anchored bandwidth trend recovers a varying sedimentation rate (Sinnesael et al. 2018)", {
  set.seed(3); t <- seq(0, 3000, 0.5)
  x <- 0.8 * sin(2 * pi * t / 41) + sin(2 * pi * t / 23.7) + 0.5 * sin(2 * pi * t / 405) + 0.3 * rnorm(length(t))
  rate <- 0.02 + 0.01 * sin(2 * pi * t / 1500)
  z <- cumsum(c(0, rate[-1] * diff(t))); zg <- seq(0, max(z), 0.02)
  d <- data.frame(zg, approx(z, x, zg)$y); truth <- approx(z, rate, zg)$y
  z0 <- 5; r0 <- approx(zg, truth, z0)$y            # the obliquity line identified at 5 m: 41 kyr at 1 / (41 r0) cycles/m
  sm <- spectral_moments(d, window = 5, step = 0.1, anchor = c(z0, 41, 1 / (41 * r0)), freq_range = c(0.1, 8), noise_floor = TRUE)
  r <- sm$result; tt <- approx(zg, truth, r[, "depth"])$y
  expect_true(all(c("moment_trend", "sedrate_anchored") %in% colnames(r)))
  expect_gt(cor(r[, "sedrate_anchored"], tt), 0.9)
  expect_gt(mean(abs(r[, "sedrate_anchored"] / tt - 1) < 0.2), 0.7)
  sm2 <- spectral_moments(d, window = 5, step = 0.1, anchor = list(depth = z0, rate = r0), moment = "mean_frequency", trend = "lowess")
  expect_gt(cor(sm2$result[, "sedrate_anchored"], tt), 0.8)
  f <- tempfile(fileext = ".png"); png(f); expect_silent(plot_spectral_moments(sm, freq_max = 8)); dev.off()
})
