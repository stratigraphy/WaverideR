test_that("all mother wavelets are defined and calibrated", {
  for (w in c(setdiff(.wr_wavelet_names, "custom"), "mexh", "gabor")) {
    d <- .wr_wavelet(w)
    expect_true(is.finite(d$peak_omega) && d$peak_omega > 0)
    expect_true(is.finite(d$psi0) && d$psi0 > 0)
  }
  expect_equal(.wr_wavelet("mexh")$name, "dog")
  expect_equal(.wr_wavelet("gabor", omega_nr = 6)$order, 6)
  # equal relative bandwidth: the mapped order reproduces the Morlet quality factor
  for (w in c("paul", "dog", "morse", "bump", "cmor", "cgau", "fbsp", "loggabor")) {
    o <- .wr_order_for_omega(w, 8)
    expect_equal(.wr_quality(w, o), 8 * sqrt(2), tolerance = 1e-3)
  }
})

test_that("wavelet agnostic Awavelet, superlets, cross wavelet and coherence", {
  x <- seq(0, 600, by = 1)
  d <- data.frame(x, sin(2 * pi * x / 40) + 0.3 * sin(2 * pi * x / 100))
  for (w in c("morlet", "paul", "cmor", "bump", "shannon", "meyer")) {
    aw <- analyse_Awavelet(d, dj = 1 / 20, lowerPeriod = 15, upperPeriod = 200, wavelet = w)
    expect_equal(aw$wavelet, w)
    expect_lt(abs(aw$Period[which.max(aw$Power.avg)] - 40), 8)
    sl <- analyse_superlet(d, upperPeriod = 200, lowerPeriod = 15, Nf = 40, c1 = 3, o = c(1, 2), wavelet = w)
    expect_lt(abs(sl$Period[which.max(sl$Power.avg)] - 40), 6)
    xs <- analyse_Xsuperlet(d, d, upperPeriod = 200, lowerPeriod = 15, Nf = 30, o = c(1, 2), wavelet = w)
    expect_equal(xs$wavelet, w)
    xw <- analyse_Xwavelet(d, d, dj = 1 / 20, lowerPeriod = 15, upperPeriod = 200, wavelet = w)
    expect_equal(xw$wavelet, w)
    expect_true(all(abs(xw$Phase) < 1e-8))
  }
  co <- analyse_wavelet_coherence(d, d, dj = 1 / 10, lowerPeriod = 15, upperPeriod = 200, wavelet = "paul",
                                  abs_window_t = 50, abs_window_s = 0.5)
  expect_equal(co$wavelet, "paul")
  # legacy Morlet paths unchanged by the new arguments
  a1 <- analyse_Awavelet(d, dj = 1 / 20, lowerPeriod = 15, upperPeriod = 200)
  a2 <- analyse_Awavelet(d, dj = 1 / 20, lowerPeriod = 15, upperPeriod = 200, wavelet = "morlet")
  expect_equal(a1$Wave, a2$Wave)
})

test_that("synchrosqueezing accepts the new analytic wavelets", {
  t <- 1:500
  x <- sin(2 * pi * t / 50)
  for (w in c("cmor", "cgau", "loggabor", "fbsp")) {
    s <- analyse_wavelet_sst(data.frame(t, x), lowerPeriod = 20, upperPeriod = 120, dj = 1 / 20, wavelet = w)
    expect_lt(abs(s$Period[which.max(s$Power.avg)] - 50), 5)
  }
})
