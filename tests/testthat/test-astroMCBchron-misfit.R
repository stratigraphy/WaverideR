## the uncertainty of an anchored age model may never have a local minimum between two dated horizons, and the stored
## curves of misfit = "window" (default) keep time everywhere
nmin <- function(v) {
  sum(diff(sign(round(diff(v), 10))) == 2)
}
interior_minima <- function(fit, dates) {
  s <- fit$summary[order(fit$summary$position), ]
  zz <- sort(unique(dates$position))
  n <- numeric(length(zz) - 1)
  for (i in seq_len(length(zz) - 1)) {
    k <- s$position > zz[i] & s$position < zz[i + 1]
    n[i] <- nmin(s$sd2_kyr[k]) + nmin((s$median - s$lower)[k]) + nmin((s$upper - s$median)[k])
  }
  sum(n)
}
sim_record <- function() {
  set.seed(1)
  rate <- function(z) {
    0.02 * (1 + 0.5 * sin(2 * pi * z / 40))
  }
  zz <- seq(0, 100, 0.05)
  age_z <- 400000 + c(0, cumsum(diff(zz) / rate(zz[-1]))) + ifelse(zz > 60, 300, 0)
  z_dates <- c(10, 35, 55, 70, 95)
  dates <- data.frame(id = paste0("ash", 1:5), position = z_dates,
                      age = approx(zz, age_z, z_dates)$y / 1000 + rnorm(5, 0, 0.08), sd = 0.08, thickness = 0.5)
  trk <- function(z) {
    period <- matrix(NA_real_, length(z), 20)
    for (i in 1:20) {
      period[, i] <- 405 * rate(z) * exp(rnorm(1, 0, 0.02))
    }
    list(position = z, period = period)
  }
  list(dates = dates, tracks = list(trk(seq(0, 60, 0.25)), trk(seq(60, 100, 0.25))),
       hi = data.frame(position = 60, mean = 280, sd = 60))
}

test_that("misfit = 'localised' and 'window' give no interior uncertainty minimum between dates (with a hiatus)", {
  skip_on_cran()
  r <- sim_record()
  fits <- list()
  for (m in c("localised", "window", "stretch")) {
    fits[[m]] <- astroMCBchron(r$dates, r$tracks, hiatuses = r$hi, n_iter = 4000, loo = FALSE, n_tracks = 200,
                               predict_positions = seq(0, 100, 0.5), seed = 2, verbose = FALSE, misfit = m)
  }
  expect_equal(interior_minima(fits$localised, r$dates), 0)
  expect_equal(interior_minima(fits$window, r$dates), 0)
  ## same posterior and mean age model as the rigid stretch
  expect_equal(fits$localised$summary$mean, fits$stretch$summary$mean)
  expect_equal(fits$window$summary$mean, fits$stretch$summary$mean)
  ## between two dates the uncertainty is never below the lower dated end
  for (m in c("localised", "window")) {
    s <- fits[[m]]$summary[order(fits[[m]]$summary$position), ]
    zd <- sort(r$dates$position)
    for (i in 1:4) {
      k <- s$position >= zd[i] & s$position <= zd[i + 1]
      ends <- s$sd2_kyr[k][c(1, sum(k))]
      expect_true(all(s$sd2_kyr[k] >= min(ends) - 1e-6))
    }
  }
})

test_that("misfit = 'window' keeps time in every step and stays close to 'localised'", {
  skip_on_cran()
  r <- sim_record()
  fl <- astroMCBchron(r$dates, r$tracks, hiatuses = r$hi, n_iter = 4000, loo = FALSE, n_tracks = 200,
                      predict_positions = seq(0, 100, 0.5), seed = 2, verbose = FALSE, misfit = "localised")
  fw <- astroMCBchron(r$dates, r$tracks, hiatuses = r$hi, n_iter = 4000, loo = FALSE, n_tracks = 200,
                      predict_positions = seq(0, 100, 0.5), seed = 2, verbose = FALSE)
  expect_identical(fw$settings$misfit, "window")
  dt <- diff(fw$model)
  m <- rowMeans(dt)
  hs <- which(fw$predictPositions[-1] > 60 & fw$predictPositions[-length(fw$predictPositions)] <= 60)
  empty <- 0
  for (j in seq_len(ncol(dt))) {
    empty <- empty + sum(dt[-hs, j] < 0.05 * m[-hs])
  }
  expect_equal(empty, 0)
  zd <- range(r$dates$position)
  k <- fw$predictPositions >= zd[1] & fw$predictPositions <= zd[2]
  expect_lt(max(abs(fw$summary$sd2_kyr[k] - fl$summary$sd2_kyr[k]) / fl$summary$sd2_kyr[k]), 0.15)
  ## analytic band: symmetric about the mean, 1.96 sd wide; the curve percentiles are returned as well
  s <- fw$summary
  expect_equal((s$upper - s$mean) * 1000, 1.959964 * s$sd2_kyr / 2, tolerance = 1e-8)
  expect_true(all(c("lower_mc", "upper_mc", "guard_raise_kyr") %in% names(s)))
  expect_true(all(s$guard_raise_kyr >= 0))
})
