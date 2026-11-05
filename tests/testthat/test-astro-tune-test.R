test_that("rate tracking recovers a linearly changing sedimentation rate", {
  skip_on_cran()
  set.seed(11)
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(NULL)
  dx <- 0.15
  z <- seq(0, 200, by = dx)
  rate <- 9 - 3 * z / max(z)
  tt <- c(0, cumsum(diff(z) * 100 / ((rate[-1] + rate[-length(z)]) / 2)))
  x <- 0
  for (p in c(tg$E[1:3], tg$O[1], tg$P[1:2])) {
    x <- x + cos(2 * pi * tt / p + stats::runif(1, 0, 2 * pi))
  }
  x <- x + 0.5 * stats::rnorm(length(z))
  r <- astro_tune_test(cbind(z, x), sedmin = 3, sedmax = 20, numsed = 30, n_sur = 3, n_paths = 10)
  expect_s3_class(r, "astro_tune_test")
  # a few pure sinusoids (two precession terms only): within 15% (ETP records: ~1-3%)
  expect_lt(abs(r$summary$duration_kyr / max(tt) - 1), 0.15)
  expect_true(all(r$summary$p_loglik > 0 & r$summary$p_loglik <= 1))
  expect_equal(nrow(r$tests), 1)
})

test_that("multivariate surrogates keep each spectrum and the cross-correlation", {
  set.seed(12)
  ns <- asNamespace("WaverideR")
  a <- as.numeric(stats::filter(rnorm(512), 0.8, method = "recursive"))
  Y <- cbind(a, a + rnorm(512, sd = 0.3))
  S <- ns$wr_att_surrogate(Y)
  expect_equal(Mod(stats::fft(S[, 1] - mean(S[, 1]))), Mod(stats::fft(Y[, 1] - mean(Y[, 1]))),
      tolerance = 1e-8)
  expect_equal(cor(S[, 1], S[, 2]), cor(Y[, 1], Y[, 2]), tolerance = 0.05)
})

test_that("plot and print of astro_tune_test leave par unchanged", {
  skip_on_cran()
  set.seed(13)
  z <- seq(0, 60, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(z)), 0.7, method = "recursive"))
  r <- astro_tune_test(cbind(z, x), sedmin = 1, sedmax = 5, numsed = 15, n_sur = 3, n_paths = 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p0 <- graphics::par(no.readonly = TRUE)
  plot(r)
  expect_identical(p0, graphics::par(no.readonly = TRUE))
  expect_output(print(r), "mean rate")
})
