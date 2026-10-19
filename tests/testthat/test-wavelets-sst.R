test_that("all mother wavelets find the cycle and Morlet is unchanged", {
  d <- .wr_wavelet("morlet", omega_nr = 8)
  N <- 256; dt <- 2
  ok <- 1:floor(N / 2) * 2 * pi / (N * dt)
  om <- c(0, ok, -ok[floor((N - 1) / 2):1])
  s <- 3.7
  leg <- pi^(1/4) * sqrt(2 * s / dt) * exp(-((s * om - 8)^2 / 2) * (om > 0)) * (om > 0)
  expect_equal(d$daughter(s, om, dt), leg, tolerance = 1e-12)
  expect_equal(.wr_wavelet("paul", 4)$psi0, 1.079, tolerance = 1e-3)
  expect_equal(.wr_wavelet("dog", 2)$psi0, 0.867, tolerance = 1e-3)
  t <- seq(0, 1000, by = 2)
  x <- sin(2 * pi * t / 50)
  for (w in c("morlet", "paul", "dog", "morse", "bump")) {
    r <- analyse_wavelet(data.frame(t, x), lowerPeriod = 10, upperPeriod = 400, dj = 1 / 20, wavelet = w)
    expect_equal(r$wavelet, w)
    pk <- r$Period[which.max(rowMeans(r$Power[, 100:400]))]
    expect_true(abs(pk - 50) / 50 < 0.2, info = w)
  }
  expect_error(analyse_wavelet(data.frame(t, x), wavelet = "nope"))
})

test_that("synchrosqueezing sharpens ridges and keeps the reconstruction", {
  set.seed(1)
  t <- seq(0, 1000, by = 2)
  x <- sin(2 * pi * t / 45) + sin(2 * pi * t / 60) + rnorm(length(t), sd = 0.2)
  a <- list(lowerPeriod = 20, upperPeriod = 150, dj = 1 / 40)
  cwt <- do.call(analyse_wavelet, c(list(data.frame(t, x)), a))
  sst <- do.call(analyse_wavelet_sst, c(list(data.frame(t, x)), a))
  expect_true(isTRUE(sst$synchrosqueezed))
  expect_s3_class(sst, "analyze.wavelet")
  avg <- rowMeans(sst$Power[, 100:400])
  expect_true(abs(sst$Period[which.max(avg)] - 45) < 3 || abs(sst$Period[which.max(avg)] - 60) < 3)
  rec <- function(o) colSums(Re(o$Wave) / sqrt(o$Scale)) * o$dj * sqrt(o$dt) / o$psi0
  expect_gt(cor(rec(cwt), rec(sst)), 0.99)
  # the two cycles are separated by a deeper valley than in the scalogram
  valley <- function(o) { v <- rowMeans(o$Power[, 100:400]); i <- which.min(abs(o$Period - 52)); max(v) / v[i] }
  expect_gt(valley(sst), valley(cwt))
  expect_error(analyse_wavelet_sst(data.frame(t, x), wavelet = "dog"))
  expect_identical(analyze_wavelet_sst, analyse_wavelet_sst)
  p <- analyse_wavelet(data.frame(t, x), lowerPeriod = 20, upperPeriod = 150, dj = 1 / 20, wavelet = "paul")
  expect_silent(suppressMessages(ws <- wavelet_significance(p)))
  expect_error(wavelet_significance(sst))
})
