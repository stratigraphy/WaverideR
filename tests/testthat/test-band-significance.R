test_that("exact band-averaged significance reproduces Monte Carlo moments for several omega", {
  skip_on_cran()
  set.seed(11)
  for (om in c(6, 10)) {
    N <- 1024
    pr <- c(8, 32)
    vals <- NULL
    for (r in 1:10) {
      w <- analyse_wavelet(cbind(1:N, rnorm(N)), dj = 0.125, lowerPeriod = 4, upperPeriod = 64,
          verbose = FALSE, omega_nr = om)
      b <- wavelet_band_significance(w, pr, background = "white")
      vals <- c(vals, b$power[b$valid])
    }
    expect_equal(mean(vals), b$background_mean, tolerance = 0.05)
    expect_equal(2 * mean(vals)^2 / var(vals), b$dof, tolerance = 0.2)
  }
})

test_that("false positive rate is near nominal for omega_nr 4, 6 and 10 (exact) and 6 (torrence)", {
  skip_on_cran()
  set.seed(12)
  cases <- list(list(4, "exact"), list(6, "exact"), list(10, "exact"), list(6, "torrence"), list(10, "torrence"))
  for (cs in cases) {
    hit <- 0
    tot <- 0
    for (r in 1:60) {
      N <- 1024
      y <- as.numeric(arima.sim(list(ar = 0.4), N))
      w <- analyse_wavelet(cbind(1:N, y), dj = 0.125, lowerPeriod = 4, upperPeriod = 128, verbose = FALSE,
          omega_nr = cs[[1]])
      b <- wavelet_band_significance(w, c(8, 24), method = cs[[2]])
      hit <- hit + sum(b$significant)
      tot <- tot + sum(b$valid)
    }
    expect_true(abs(hit / tot - 0.05) < 0.015, info = paste(cs[[1]], cs[[2]], hit / tot))
  }
})

test_that("a strong cycle is detected inside the band and not outside it", {
  set.seed(13)
  z <- seq(0, 100, by = 0.1)
  y <- sin(2 * pi * z / 4) + 0.5 * as.numeric(arima.sim(list(ar = 0.3), length(z)))
  w <- analyse_wavelet(cbind(z, y), dj = 1 / 8, lowerPeriod = 0.4, upperPeriod = 30, verbose = FALSE, omega_nr = 8)
  inside <- wavelet_band_significance(w, c(3, 5))
  outside <- wavelet_band_significance(w, c(0.6, 1.2))
  expect_gt(mean(inside$significant[inside$valid]), 0.9)
  expect_lt(mean(outside$significant[outside$valid]), 0.15)
  expect_s3_class(inside, "wr_band_significance")
  expect_output(print(inside), "scale-averaged")
  expect_true(is.finite(inside$dof) && inside$dof > 2)
  expect_error(wavelet_band_significance(w, c(5, 3)))
})
