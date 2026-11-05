test_that("chirplet synchrosqueezing is sharper than the superlet-SST on a chirping record", {
  set.seed(1)
  z <- seq(0, 600, by = 0.5)
  sr <- 1 + z / 600
  t <- cumsum(c(0, diff(z) / sr[-1]))
  x <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 41) + rnorm(length(z), sd = 0.5)
  cs <- analyse_chirplet_sst(data.frame(z, x), lowerPeriod = 10, upperPeriod = 200, dj = 1 / 50)
  ss <- analyse_superlet_sst(data.frame(z, x), lowerPeriod = 10, upperPeriod = 200, dj = 1 / 50)
  expect_s3_class(cs, "analyze.wavelet")
  expect_equal(dim(cs$Power), dim(ss$Power))
  expect_equal(dim(cs$Chirp), dim(cs$Power))
  conc <- function(m) {
    tot <- 0
    for (p in c(23.7, 41)) {
      for (i in seq_along(z)) {
        b <- abs(log(m$Period / (p * sr[i]))) < log(1.03)
        tot <- tot + sum(m$Power[b, i])
      }
    }
    tot / sum(m$Power)
  }
  expect_gt(conc(cs), conc(ss))
  expect_gt(conc(cs), 0.95)
  # the chirp map has the sign of the change of the period (period grows with depth here)
  i <- which(z == 300)
  j <- which.min(abs(cs$Period - 23.7 * sr[i]))
  expect_gt(cs$Chirp[j, i], 0)
})

test_that("the iterative tracker keeps or improves the age model and returns its iterations", {
  set.seed(3)
  t <- seq(0, 1500, by = 0.5)
  ecc <- 0.4 * cos(2 * pi * t / 405) + 0.25 * cos(2 * pi * t / 125)
  x_t <- (1 + ecc) * (sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19)) + 0.7 * sin(2 * pi * t / 41)
  sr <- ifelse(t < 700, 1.5, 2.6)
  zt <- c(0, cumsum(sr[-1] * diff(t)))
  z <- seq(0, max(zt), by = 1)
  x <- stats::approx(zt, x_t, z)$y + rnorm(length(z), sd = 1)
  res <- track_astro_iterative(data.frame(z, x), sedrate = c(0.8, 5), n_iter = 2, n_samples = 20, verbose = FALSE)
  expect_s3_class(res, "wr_track")
  its <- attr(res, "iterations")
  expect_true(nrow(its) >= 2)
  expect_true(!is.null(attr(res, "warped")$scalogram))
  srz <- stats::approx(zt, sr, z)$y
  j <- which(z > 150 & z < max(z) - 150)
  expect_lt(stats::median(abs(log2(res$sedrate[j] / srz[j]))), 0.1)
  expect_true(all(res$time_upper >= res$time_lower))
})

test_that("L1 path smoothing keeps jumps and removes wiggles", {
  v <- c(rep(0, 40), rep(1, 40)) + sin(seq_len(80)) * 0.05
  u <- .wr_l1_smooth(v, rep(1, 80), lambda = 1, order = 1)
  expect_lt(sd(u[1:40]), sd(v[1:40]))
  expect_gt(mean(u[41:80]) - mean(u[1:40]), 0.9)
  m <- analyse_superlet_sst(data.frame(seq(0, 500, by = 1), rnorm(501)), lowerPeriod = 10, upperPeriod = 100, dj = 1 / 20)
  expect_silent(suppressMessages(track_astro_hmm(m, astro_periods = c(41, 23.7), sedrate = c(0.5, 2), n_samples = 5, smooth_path = "tv")))
})

test_that("detect_transients finds a spike", {
  set.seed(2)
  z <- seq(0, 1000, by = 1)
  x <- sin(2 * pi * z / 30) + 0.5 * sin(2 * pi * z / 70) + rnorm(length(z), sd = 0.3)
  x[500:502] <- x[500:502] + 5
  tr <- detect_transients(data.frame(z, x), lowerPeriod = 10, upperPeriod = 200)
  expect_true(any(abs(tr$transients$depth - 500) <= 3))
  expect_equal(tr$transients$depth[1], 500, tolerance = 3)
})

test_that("plot_astro_evidence draws the evidence and the posterior", {
  set.seed(5)
  t <- seq(0, 1200, by = 0.5)
  x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
  sr <- 2 + 0.6 * sin(2 * pi * t / 500)
  zt <- c(0, cumsum(sr[-1] * diff(t)))
  z <- seq(0, max(zt), by = 1)
  x <- stats::approx(zt, x_t, z)$y + rnorm(length(z), sd = 1)
  w <- analyse_wavelet(data.frame(z, x), lowerPeriod = 10, upperPeriod = 400, dj = 1 / 40, omega_nr = 10, verbose = FALSE)
  pdf(NULL)
  r <- suppressMessages(plot_astro_evidence(w, sedrate = c(0.5, 8), reference = data.frame(z, stats::approx(zt, sr, z)$y)))
  expect_equal(nrow(r$map), length(r$sedrate))
  tr <- suppressMessages(track_astro_hmm(w, sedrate = c(0.5, 8), n_samples = 5))
  expect_true(!is.null(attr(tr, "posterior")))
  r2 <- plot_astro_evidence(track = tr, type = "posterior")
  expect_equal(dim(r2$map), dim(attr(tr, "posterior")$post))
  expect_silent(plot_wavelet(w, whiten = TRUE))
  dev.off()
})
