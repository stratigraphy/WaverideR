## track_period_wavelet() needs mouse input; the clicks are scripted by replacing the two input wrappers.
make_wt <- function() {
  z <- seq(500, 1500, by = 1)          # deep record: depths > 1024 m overflowed 2^depth in the old code
  x <- sin(2 * pi * z / (40 + 10 * (z - 500) / 1000)) + 0.3 * sin(2 * pi * z / 9)
  analyze_wavelet(data = cbind(z, x), dj = 1/20, lowerPeriod = 4, upperPeriod = 120, verbose = FALSE, omega_nr = 6)
}

run_track <- function(wt, clicks_depth, clicks_period, horizontal, track_peaks, fun = "track_period_wavelet") {
  i <- 0
  check_panel <- function() {          # the clicks must land on the scalogram panel
    usr <- graphics::par("usr")
    dr <- if (horizontal) { usr[1:2] } else { usr[3:4] }
    expect_true(min(dr) <= min(wt$x) + 1 && max(dr) >= max(wt$x) - 1)
  }
  locator <- function() {
    i <<- i + 1
    if (i == 1) {
      check_panel()
    }
    if (i > length(clicks_depth)) {
      return(NULL)
    }
    lp <- log2(clicks_period[i])
    if (horizontal) {
      list(x = clicks_depth[i], y = lp)
    } else {
      list(x = lp, y = clicks_depth[i])
    }
  }
  identify <- function(x, y) {
    i <<- i + 1
    if (i == 1) {
      check_panel()
    }
    if (i > length(clicks_depth)) {
      return(integer(0))
    }
    lp <- log2(clicks_period[i])
    px <- if (horizontal) { clicks_depth[i] } else { lp }
    py <- if (horizontal) { lp } else { clicks_depth[i] }
    usr <- graphics::par("usr")
    which.min(((x - px) / diff(usr[1:2]))^2 + ((y - py) / diff(usr[3:4]))^2)
  }
  local_mocked_bindings(.tpw_locator = locator, .tpw_identify = identify, wr_open_device = function(...) invisible(TRUE))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (fun == "track_period_wavelet") {
    track_period_wavelet(wavelet = wt, astro_cycle = 405, track_peaks = track_peaks, plot_horizontal = horizontal)
  } else {
    track_period(scalogram = wt, astro_cycle = 405, track_peaks = track_peaks, plot_horizontal = horizontal)
  }
}

test_that("track_peaks = FALSE returns depth, period, sedrate in both layouts", {
  wt <- make_wt()
  d <- c(seq(560, 1440, by = 80), 1100.4)
  p <- 40 + 25 * (d - 500) / 1000
  for (h in c(TRUE, FALSE)) {
    out <- run_track(wt, d, p, horizontal = h, track_peaks = FALSE)
    expect_named(out, c("depth", "period", "sedrate"))
    expect_equal(out$depth, sort(round(d)))
    expect_equal(out$period, p[order(d)], tolerance = 1e-8)
    expect_equal(out$sedrate, out$period / 405 * 100)
    expect_true(all(is.finite(as.matrix(out))))
    cmp <- suppressWarnings(loess_auto(time_series = out[, 1:2], genplot = FALSE, print_span = FALSE))
    expect_true(all(is.finite(curve2time(tracked_cycle_curve = cmp, tracked_cycle_period = 405, genplot = FALSE)[, 2])))
  }
})

test_that("clicking a pick again removes it", {
  wt <- make_wt()
  out <- run_track(wt, c(700, 1000, 700, 1300), c(42, 45, 42, 48), horizontal = TRUE, track_peaks = FALSE)
  expect_equal(out$depth, c(1000, 1300))
})

test_that("track_peaks = TRUE snaps to the power maxima in both layouts", {
  wt <- make_wt()
  for (h in c(TRUE, FALSE)) {
    out <- run_track(wt, c(600, 1000, 1400), c(48, 58, 68), horizontal = h, track_peaks = TRUE)
    expect_equal(nrow(out), 3)
    zs <- c(600, 1000, 1400)
    ridge <- numeric(3)
    b <- wt$Period > 20
    for (m in 1:3) {
      j <- which.min(abs(wt$x - zs[m]))
      ridge[m] <- wt$Period[b][which.max(wt$Power[b, j])]
    }
    expect_true(all(abs(out$depth - c(600, 1000, 1400)) <= 30))
    expect_true(all(abs(out$period / ridge - 1) < 0.1))
  }
})

test_that("track_period() (wavelet and superlet) returns depth, period, sedrate in both layouts", {
  wt <- make_wt()
  z <- wt$x
  sl <- analyse_superlet(data.frame(z, wt$y), upperPeriod = 120, lowerPeriod = 4, Nf = 40, c1 = 3, o = c(1, 3),
      order_scaling = "log2", verbose = FALSE)
  d <- c(seq(560, 1440, by = 80))
  p <- 40 + 25 * (d - 500) / 1000
  for (obj in list(wt, sl)) {
    for (h in c(TRUE, FALSE)) {
      out <- run_track(obj, d, p, horizontal = h, track_peaks = FALSE, fun = "track_period")
      expect_equal(out$depth, d)
      expect_equal(out$period, p, tolerance = 1e-8)
      out2 <- run_track(obj, c(700, 1300), c(50, 64), horizontal = h, track_peaks = TRUE, fun = "track_period")
      expect_equal(nrow(out2), 2)
      expect_true(all(is.finite(as.matrix(out2))))
    }
  }
})
