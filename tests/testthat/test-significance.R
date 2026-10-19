noise_wt <- function(phi = 0.5, n = 512, seed = 1, omega = 8) {
  set.seed(seed)
  z <- seq(0, by = 0.1, length.out = n)
  y <- as.numeric(stats::arima.sim(list(ar = phi), n))
  analyse_wavelet(data.frame(z, y), dj = 1 / 10, lowerPeriod = 0.4, upperPeriod = 16, verbose = FALSE, omega_nr = omega)
}

test_that("analytic (Torrence and Compo) levels give about 5 percent false positives", {
  fr <- vapply(1:12, function(k) {
    s <- wavelet_significance(noise_wt(0.5, seed = k))
    mean(s$significant[s$valid])
  }, 1)
  expect_gt(mean(fr), 0.02)
  expect_lt(mean(fr), 0.08)
  s <- wavelet_significance(noise_wt(0.5))
  expect_s3_class(s, "wr_significance")
  expect_equal(length(s$crit), length(s$period))
  expect_true(all(s$crit > 0))
  expect_error(wavelet_significance(noise_wt(), background = "powerlaw"), "supports")
})

test_that("surrogate levels agree roughly with the analytic ones and do not depend on cores", {
  w <- noise_wt(0.5)
  a <- wavelet_significance(w)
  s <- wavelet_significance(w, method = "surrogate", n_sim = 24, seed = 5, candidates = 4)
  expect_gt(stats::median(a$crit / s$crit), 0.7)
  expect_lt(stats::median(a$crit / s$crit), 1.5)
  s2 <- wavelet_significance(w, method = "surrogate", n_sim = 24, seed = 5, run_multicore = TRUE, n_cores = 2)
  expect_equal(s$crit, s2$crit)
  expect_equal(nrow(s$candidates), 1)
  for (b in c("white", "ar1_corrected", "powerlaw")) {
    sb <- wavelet_significance(w, method = "surrogate", background = b, n_sim = 6, seed = 1, notch_periods = 4)
    expect_equal(length(sb$crit), length(w$Period))
  }
})

test_that("significance works for superlets and can be plotted", {
  set.seed(2)
  z <- seq(0, 40, by = 0.1)
  y <- sin(2 * pi * z / 4) + as.numeric(stats::arima.sim(list(ar = 0.5), length(z)))
  sl <- analyse_superlet(data.frame(z, y), 15, 0.4, 25, c1 = 3, o = c(1, 3), order_scaling = "log2")
  expect_error(wavelet_significance(sl, method = "analytic"), "only valid for wavelets")
  s <- wavelet_significance(sl, n_sim = 12, seed = 1)
  expect_equal(dim(s$significant), dim(sl$Power))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_error(plot_superlet(sl, dev_new = FALSE, add_signif = s, add_avg = TRUE), NA)
  expect_error(plot_superlet(sl, dev_new = FALSE, add_signif = s, plot_horizontal = FALSE), NA)
  w <- noise_wt()
  expect_error(plot_wavelet(w, dev_new = FALSE, add_signif = wavelet_significance(w), add_avg = TRUE), NA)
  expect_error(plot_wavelet(w, dev_new = FALSE, add_signif = wavelet_significance(w), plot_horizontal = FALSE), NA)
})

test_that("superlet cone of influence is smooth and the transform can run on several cores", {
  z <- seq(0, 40, by = 0.1)
  d <- data.frame(z, sin(2 * pi * z / 4))
  a <- analyse_superlet(d, 15, 0.4, 30, c1 = 3, o = c(1, 5), order_scaling = "log2")
  expect_lt(max(abs(log(a$coi[-1] / a$coi[-length(a$coi)]))), log(1.2))      # no steps
  b <- analyse_superlet(d, 15, 0.4, 30, c1 = 3, o = c(1, 5), order_scaling = "log2", run_multicore = TRUE, n_cores = 2)
  expect_equal(a$Wave, b$Wave)
})
