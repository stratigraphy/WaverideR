test_that("extended wavelet families are available and find the cycle", {
  expect_equal(.wr_resolve("db4")$name, "db")
  expect_equal(.wr_resolve("db4")$order, 4)
  expect_equal(.wr_resolve("cmor1.5-1")$order, c(1.5, 1))
  expect_equal(.wr_resolve("gaus2")$name, "dog")
  expect_equal(.wr_resolve("cmhat")$order, 2)
  expect_error(.wr_resolve("nonsense"))
  # Daubechies filter of order 2 (known coefficients)
  h <- .wr_daub_filter(2) * sqrt(2)
  expect_equal(sort(h), sort(c(0.4829629, 0.8365163, 0.2241439, -0.1294095)), tolerance = 1e-6)
  # Symlet 8 (known coefficients)
  h <- .wr_daub_filter(8, "sym") * sqrt(2)
  expect_equal(h[9], 0.7771857517, tolerance = 1e-6)
  t <- 1:600
  d <- data.frame(t, sin(2 * pi * t / 50))
  for (w in c("morlet0", "morl", "hhat", "poisson", "gammatone", "hann", "haar", "db4", "sym6", "battle")) {
    r <- analyse_wavelet(d, lowerPeriod = 15, upperPeriod = 200, dj = 1 / 20, wavelet = w)
    expect_lt(abs(log2(r$Period[which.max(r$Power.avg)] / 50)), 0.15)
  }
  f <- function(w) ifelse(w > 0, exp(-(w - 7)^2 / 2), 0)
  r <- analyse_wavelet(d, lowerPeriod = 15, upperPeriod = 200, dj = 1 / 20, wavelet = f)
  expect_equal(r$wavelet, "custom")
  expect_lt(abs(r$Period[which.max(r$Power.avg)] / 50 - 1), 0.05)
  expect_error(analyse_wavelet_sst(d, wavelet = "db4"))
})

test_that("second-order synchrosqueezing sharpens a chirp", {
  z <- 0:1499
  f <- 1 / 60 + (1 / 20 - 1 / 60) * z / 1500
  d <- data.frame(z, cos(2 * pi * cumsum(f)))
  conc <- function(w) {
    j <- which(z > 200 & z < 1300)
    mean(sapply(j, function(i) sum(w$Power[abs(w$Period * f[i] - 1) <= 0.03, i]) / sum(w$Power[, i])))
  }
  a <- analyse_wavelet_sst(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = 20)
  b <- analyse_wavelet_sst(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = 20, sst_order = 2)
  expect_gt(conc(b), conc(a) + 0.1)
  e <- analyse_wavelet_sst(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = 8, type = "set")
  expect_equal(e$type, "set")
  ad <- analyse_wavelet_sst(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = c(6, 20))
  expect_true(all(ad$omega_selected %in% c(6, 20)))
  expect_gt(conc(ad), conc(a))
})

test_that("superlet synchrosqueezing separates adjacent close cycles in noise", {
  set.seed(3)
  z <- 0:1199
  Pt <- ifelse(z < 600, 40, 46)
  d <- data.frame(z, sin(2 * pi * cumsum(1 / Pt)) + rnorm(length(z), sd = 1))
  s <- analyse_superlet_sst(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 40, o = 4)
  expect_true(isTRUE(s$superlet))
  r1 <- s$Period[apply(s$Power[, z > 150 & z < 500], 2, which.max)]
  r2 <- s$Period[apply(s$Power[, z > 700 & z < 1050], 2, which.max)]
  expect_lt(abs(median(r1) / 40 - 1), 0.04)
  expect_lt(abs(median(r2) / 46 - 1), 0.04)
  expect_identical(analyze_superlet_sst, analyse_superlet_sst)
})

test_that("multisynchrosqueezing, reassignment, ConceFT and significance work", {
  set.seed(5)
  z <- 0:799
  d <- data.frame(z, sin(2 * pi * z / 40) + rnorm(length(z), sd = 0.8))
  fr <- function(w) sum(w$Power[abs(w$Period / 40 - 1) <= 0.03, 100:700]) / sum(w$Power[, 100:700])
  a <- analyse_wavelet_sst(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 40, omega_nr = 12)
  m <- analyse_wavelet_sst(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 40, omega_nr = 12, n_iter = 4)
  expect_gte(fr(m), fr(a) - 0.01)
  r <- analyse_wavelet_sst(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 40, omega_nr = 12, type = "reassign")
  expect_gt(fr(r), 0.5)
  cf <- analyse_conceft(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 40, n_comb = 4)
  expect_true(isTRUE(cf$conceft))
  expect_gt(fr(cf), 0.5)
  sig <- sst_significance(m, n_sim = 5)
  expect_true(is.logical(sig$significant))
  expect_gt(sum(sig$Power[abs(sig$Period / 40 - 1) <= 0.03, ]), 0.5 * sum(m$Power[abs(m$Period / 40 - 1) <= 0.03, ]))
  expect_lt(sum(sig$Power), sum(m$Power) + 1e-9)
  expect_silent(suppressMessages({pdf(NULL); plot_wavelet(sig); dev.off()}))
})
