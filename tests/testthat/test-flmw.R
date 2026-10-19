test_that("peak finder returns alternating maxima", {
  p <- WaverideR:::.flmw_peaks(1:20, sin(seq(0, 6 * pi, length.out = 20)))
  expect_true(length(p$period) >= 2)
  expect_true(all(diff(p$period) > 0))
})

test_that("classic metrics find a perfect match", {
  pk <- c(4.06, 1.1, 0.41)
  m <- WaverideR:::.flmw_metrics(pk, c(406, 110, 41), 1)
  expect_equal(m[1, 3], 3)
  expect_lt(m[1, 1], 1e-8)
  expect_gt(m[1, 2], 1 - 1e-8)
})

test_that("flmw methods and flmw_path recover a constant sedimentation rate", {
  set.seed(1)
  t <- seq(0, 600, by = 1)
  s <- cos(2 * pi * t / 41) + cos(2 * pi * t / 23 + 1) + cos(2 * pi * t / 100 + 2) + rnorm(length(t), 0, 0.5)
  d <- data.frame(depth = t * 0.03, y = s)            # 3 cm/kyr
  w <- analyse_wavelet(d, dj = 1/10, lowerPeriod = 0.1, upperPeriod = 8, verbose = FALSE, omega_nr = 8)
  for (m in c("classic", "weighted", "template")) {
    r <- flmw(w, 1, 6, 0.5, c(100, 41, 23), method = m)
    expect_equal(dim(r$score_mat), c(length(seq(1, 6, 0.5)), length(r$depth)))
    expect_equal(r$method, m)
  }
  r <- flmw(w, 1, 6, 0.5, c(100, 41, 23), method = "template")
  p <- flmw_path(r)
  expect_equal(median(p$path), 3)
  expect_equal(dim(p$interval), c(length(r$depth), 4))
  expect_equal(colSums(p$posterior), rep(1, length(r$depth)), tolerance = 1e-8)
  p2 <- flmw_path(r, sedrate_range = c(4, 6))
  expect_true(all(p2$path >= 4))
})
