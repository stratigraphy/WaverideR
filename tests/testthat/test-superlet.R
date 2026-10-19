sine <- function(offset = 0, per = 4, a = 5, lag = 0) {
  z <- seq(0, 100, by = 0.05)
  data.frame(z = z, y = offset + a * sin(2 * pi * (z - lag) / per))
}

test_that("superlet finds the period of a record with a large mean", {
  d <- sine(offset = 100)
  s <- analyse_superlet(d, upperPeriod = 40, lowerPeriod = 0.5, Nf = 40, c1 = 3, o = c(1, 5), order_scaling = "log2")
  mid <- which(s$x > 30 & s$x < 70)
  expect_equal(s$Period[which.max(colMeans(s$Power[mid, ]))], 4, tolerance = 0.08)   # coarse period grid (Nf = 40)
  expect_equal(range(s$order), c(1, 5))
  expect_length(s$coi, 40)
})

test_that("superlet amplitude equals the direct geometric mean of Morlet responses", {
  d <- sine()
  s <- analyse_superlet(d, 40, 0.5, 40, c1 = 3, o = c(3, 3), detrend = "none", standardize = FALSE)
  i <- 20; j <- which(s$x > 40)[1]
  Pper <- s$Period[i]; Fc <- (1 / Pper) * s$dt
  cx <- function(Nc) {
    sd <- (Nc / 2) / Fc / 2.5; wl <- 2 * floor(6 * sd / 2) + 1; off <- floor(wl / 2)
    t <- seq_len(wl) - 1 - off; g <- exp(-t^2 / (2 * sd^2)) / (sd * sqrt(2 * pi))
    g * exp(2i * pi * Fc * t) / sum(g)
  }
  amps <- vapply(1:3, function(k) {
    w <- cx(3 * k); h <- (length(w) - 1) / 2
    Mod(sum(s$y[(j - h):(j + h)] * Conj(w)))
  }, 1)
  expect_equal(Mod(s$Wave[i, j]), exp(mean(log(amps))), tolerance = 1e-6)
})

test_that("cross superlet phase has the usual sign", {
  xs <- analyse_Xsuperlet(sine(100), sine(20, a = 3, lag = 0.5), upperPeriod = 40, lowerPeriod = 0.5, Nf = 40, c1 = 3, o = c(1, 5))
  i4 <- which.min(abs(xs$Period - 4))
  mid <- which(xs$x1 > 30 & xs$x1 < 70)
  expect_equal(stats::median(xs$Phase[mid, i4]), 2 * pi * 0.5 / 4, tolerance = 0.02)
})

test_that("order scaling options give increasing orders", {
  d <- sine()
  for (sc in c("none", "log2", "linear", "sqrt", "quadratic", "power")) {
    s <- analyse_superlet(d, 40, 0.5, 20, c1 = 3, o = c(1, 6), order_scaling = sc, order_alpha = 2)
    expect_equal(range(s$order), c(1, 6))
    expect_true(all(diff(s$order) <= 1e-12))           # longest period first, so orders fall
  }
  s1 <- analyse_superlet(d, 40, 0.5, 20, c1 = 3)
  expect_true(all(s1$order == 1))
})

test_that("flmw accepts a superlet object", {
  d <- sine(per = 1.2); d$y <- d$y + sin(2 * pi * d$z / 0.5)
  s <- analyse_superlet(d, 8, 0.2, 40, c1 = 3, o = c(1, 3))
  r <- flmw(s, 1, 5, 1, c(100, 41), method = "template")
  expect_equal(dim(r$score_mat), c(5, length(s$x)))
})

test_that("superlet objects carry a cone of influence that the plots can draw", {
  d <- sine()
  s <- analyse_superlet(d, 40, 0.5, 30, c1 = 3, o = c(1, 4))
  expect_equal(length(s$coi.1), length(s$coi.2))
  expect_equal(sum(is.na(s$coi.1)), 0)                  # one polygon
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_error(plot_superlet(s, dev_new = FALSE), NA)
  expect_error(plot_superlet(s, dev_new = FALSE, plot_horizontal = FALSE), NA)
  xs <- analyse_Xsuperlet(d, d, upperPeriod = 40, lowerPeriod = 0.5, Nf = 30, c1 = 3, o = c(1, 4))
  expect_error(plot_Xsuperlet(xs, dev_new = FALSE), NA)
})

test_that("additive superlets follow c_i = c_1 + i - 1 (Moca et al. 2021)", {
  set.seed(4)
  d <- data.frame(x = seq(0, 100, by = 0.5), y = rnorm(201))
  # order 1 has c_1 cycles for both the multiplicative and the additive rule
  a <- analyse_superlet(d, 16, 1, 12, c1 = 3, o = c(1, 1), mult = TRUE)
  b <- analyse_superlet(d, 16, 1, 12, c1 = 3, o = c(1, 1), mult = FALSE)
  expect_equal(a$Power, b$Power)
  expect_equal(a$coi, b$coi)
  # higher orders widen the wavelets more slowly in the additive case (narrower band, less sharpening)
  m <- analyse_superlet(d, 16, 1, 12, c1 = 3, o = c(3, 3), mult = TRUE)
  ad <- analyse_superlet(d, 16, 1, 12, c1 = 3, o = c(3, 3), mult = FALSE)
  expect_true(all(ad$coi <= m$coi + 1e-12))
})

test_that("the cone of influence of a superlet is one polygon whose sides do not overlap", {
  set.seed(5)
  d <- data.frame(x = seq(0, 60, by = 0.5), y = rnorm(121))
  s <- analyse_superlet(d, 40, 1, 20, c1 = 3, o = c(1, 3))
  expect_false(anyNA(s$coi.1))
  half <- mean(range(d$x))
  nl <- length(s$coi.1)
  left <- s$coi.1[2:(1 + length(s$coi))]
  right <- s$coi.1[(2 + length(s$coi)):(1 + 2 * length(s$coi))]
  expect_true(all(left <= half + 1e-9) && all(right >= half - 1e-9))
})
