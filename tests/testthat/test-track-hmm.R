test_that("the HMM engine recovers a path through a gap", {
  set.seed(1)
  ns <- 60
  nt <- 200
  truep <- round(30 + 10 * sin(seq(0, 3, length.out = nt)))
  logE <- matrix(rnorm(ns * nt), ns, nt)
  for (t in 1:nt) {
    logE[truep[t] + (-1:1), t] <- logE[truep[t] + (-1:1), t] + 2
  }
  logE[, 80:120] <- rnorm(ns * 41)
  h <- .wr_hmm(logE, 2, 1e-3, n_samples = 20)
  expect_lt(mean(abs(h$path - truep)), 3)
  expect_equal(dim(h$samples), c(nt, 20))
  expect_equal(colSums(h$post), rep(1, nt), tolerance = 1e-8)
})

test_that("track_ridge_hmm gives a continuous track and age model", {
  set.seed(2)
  z <- 0:1199
  Pt <- 40 + 8 * sin(2 * pi * z / 1000)
  amp <- ifelse(z > 500 & z < 650, 0.05, 1)          # the cycle almost disappears: a gap
  x <- amp * sin(2 * pi * cumsum(1 / Pt)) + rnorm(length(z), sd = 0.4)
  s <- analyse_wavelet_sst(data.frame(z, x), lowerPeriod = 15, upperPeriod = 120, dj = 1 / 40, n_iter = 4)
  tr <- track_ridge_hmm(s, period = c(25, 70), period_of_tracked_cycle = 20, n_samples = 50)
  expect_s3_class(tr, "wr_track")
  expect_equal(nrow(tr), length(z))
  expect_true(all(is.finite(tr$period)))
  j <- 100:1100
  expect_lt(sqrt(mean(log2(tr$period[j] / Pt[j])^2)), 0.05)
  true_t <- c(0, cumsum(20 / Pt[-1]))
  expect_lt(abs(tail(tr$time, 1) / tail(true_t, 1) - 1), 0.03)
  expect_true(all(tr$time_lower <= tr$time_upper))
  expect_silent({
    pdf(NULL)
    plot_track(tr, s)
    dev.off()
  })
})

test_that("track_astro_hmm tracks the sedimentation rate with all cycles", {
  set.seed(3)
  t <- seq(0, 1200, by = 0.5)
  ecc <- 0.4 * cos(2 * pi * t / 405) + 0.25 * cos(2 * pi * t / 125)
  x_t <- (1 + ecc) * (sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19)) + 0.8 * sin(2 * pi * t / 41)
  sr <- 2 + 0.5 * sin(2 * pi * t / 600)
  zt <- c(0, cumsum(sr[-1] * diff(t)))
  z <- seq(0, max(zt), by = 1)
  x <- stats::approx(zt, x_t, z)$y + rnorm(length(z), sd = 0.5)
  w <- analyse_wavelet(data.frame(z, x), lowerPeriod = 15, upperPeriod = 1000, dj = 1 / 40, omega_nr = 12)
  tr <- track_astro_hmm(w, astro_periods = c(405, 125, 41, 23.7, 19), sedrate = c(0.8, 6), n_samples = 50)
  srz <- stats::approx(zt, sr, z)$y
  j <- which(z > 200 & z < max(z) - 200)
  expect_lt(stats::median(abs(log2(tr$sedrate[j] / srz[j]))), 0.05)
  expect_true(all(c("period_41", "period_23.7") %in% names(tr)))
  expect_lt(abs(tail(tr$time, 1) / 1200 - 1), 0.1)
})

test_that("known jump depths raise the jump probability locally", {
  zc <- seq(0, 1000, by = 10)
  r <- .wr_jump_rho(zc, 1e-4, c(500), 0.05)
  expect_equal(max(r), 0.05)
  expect_equal(sum(r > 1e-4), 3)
  expect_true(all(r[zc < 400] == 1e-4))
})

test_that("deep-time periods and the automatic workflow", {
  p0 <- astro_periods_tracking(0)
  expect_equal(p0$period[p0$term == "k+s3"], 41, tolerance = 0.01)
  expect_equal(p0$period[p0$term == "k+g5"], 23.7, tolerance = 0.01)
  p4 <- astro_periods_tracking(420)
  expect_lt(p4$period[p4$term == "k+s3"], 35)            # shorter obliquity in the Silurian
  expect_equal(p4$period[p4$term == "g2-g5"], p0$period[p0$term == "g2-g5"])
  set.seed(4)
  t <- seq(0, 900, by = 0.5)
  x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41) + 0.4 * sin(2 * pi * t / 405)
  sr <- ifelse(t < 450, 1.5, 2.5)
  zt <- c(0, cumsum(sr[-1] * diff(t)))
  z <- seq(0, max(zt), by = 1)
  x <- stats::approx(zt, x_t, z)$y + rnorm(length(z), sd = 0.5)
  res <- track_astro_auto(data.frame(z, x), sedrate = c(0.8, 4), n_samples = 20, method = "hmm")
  expect_s3_class(res$track, "wr_track")
  expect_true(all(c("sedrate_long", "flagged") %in% names(res$check)))
  srz <- stats::approx(zt, sr, z)$y
  j <- which(z > 150 & z < max(z) - 150)
  expect_lt(stats::median(abs(log2(res$track$sedrate[j] / srz[j]))), 0.05)
  expect_silent({
    pdf(NULL)
    plot_track(res$track, res$scalogram, show = c(41, 23.7), flags = res$flags)
    dev.off()
  })
  tl <- track_astro_hmm(res$scalogram, age_ma = 0, sedrate = c(0.8, 4), long_weight = 1, n_samples = 10)
  expect_true(all(is.finite(tl$sedrate)))
})

test_that("the constrained Viterbi path meets the duration constraint", {
  set.seed(1)
  ns <- 40
  nt <- 300
  logE <- matrix(rnorm(ns * nt, sd = 0.5), ns, nt)
  logE[25:35, ] <- logE[25:35, ] + 1
  zc <- seq(0, 3, length.out = nt)
  sgrid <- exp(seq(log(0.002), log(0.015), length.out = ns))
  Tm <- matrix(0.01 / sgrid, ns, nt)
  m <- zc >= 1 & zc <= 2
  p0 <- .wr_viterbi(logE, 1, 1e-4)
  free <- sum(Tm[cbind(p0, seq_len(nt))][m])
  cons <- data.frame(from = 1, to = 2, duration = 3 * free, tolerance = 0.025)
  p <- .wr_viterbi_constrained(logE, 1, 1e-4, Tm, zc, cons, bin_frac = 0.0005)
  expect_equal(sum(Tm[cbind(p, seq_len(nt))][m]), 3 * free, tolerance = 0.03)
  hc <- suppressWarnings(.wr_hmm_constrained(logE, 1, 1e-4, 5, 1 / sgrid, zc, cons[, 1:3]))
  expect_lt(abs(hc$constraints$achieved / hc$constraints$duration - 1), 0.025)
  expect_equal(dim(hc$h$post), c(ns, nt))
})

test_that("eccentricity-only period sets weight the 405-kyr cycle fully", {
  ap <- astro_periods_tracking(420, include = "eccentricity")
  expect_equal(ap$weight[ap$term == "g2-g5"], 1)
  expect_equal(astro_periods_tracking(0)$weight[1], 0.5)
})
