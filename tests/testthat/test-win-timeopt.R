test_that("win_timeOpt splits the timeOpt output into the correct rows", {
  skip_on_cran()
  set.seed(3)
  depth <- seq(0, 8, by = 0.1)
  x <- cumsum(rnorm(length(depth)))
  dat <- cbind(depth, x)
  res <- win_timeOpt(dat, window_size = 4, sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      normalize = FALSE, bias_correction = "none")
  w <- 41
  rows <- res$x >= res$x[w] - 2 - 1e-9 & res$x <= res$x[w] + 2 + 1e-9
  ref <- astrochron::timeOpt(dat[rows, ], sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      roll = 10^6, targetE = c(405.7, 130.7, 123.8, 98.9, 94.9), output = 1, genplot = FALSE,
      verbose = FALSE, limit = FALSE)
  expect_equal(unname(res$r_2_envelope[, w]), unname(ref[, 2]))
  expect_equal(unname(res$r_2_power[, w]), unname(ref[, 3]))
  expect_equal(unname(res$r_2_opt[, w]), unname(ref[, 4]))
  expect_null(res$z_opt)
})

test_that("bias correction keeps the raw values and returns z-scores", {
  skip_on_cran()
  set.seed(4)
  depth <- seq(0, 8, by = 0.1)
  x <- cumsum(rnorm(length(depth)))
  dat <- cbind(depth, x)
  set.seed(1)
  r0 <- win_timeOpt(dat, window_size = 4, sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      bias_correction = "none", step = 4)
  set.seed(1)
  r1 <- win_timeOpt(dat, window_size = 4, sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      bias_correction = "null", n_null = 8, step = 4, min_cycles = 0)
  set.seed(1)
  r2 <- win_timeOpt(dat, window_size = 4, sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      bias_correction = "null", n_null = 8, step = 4, min_cycles = 0)
  expect_identical(r0$r_2_opt, r1$r_2_opt)
  expect_identical(r1$z_opt, r2$z_opt)
  expect_equal(dim(r1$z_opt), dim(r1$r_2_opt))
  expect_true(all(is.finite(r1$z_opt[r1$testable])))
  expect_equal(r1$z_opt, (r1$z_envelope + r1$z_power) / sqrt(2))
  expect_equal(ncol(r1$noise_param), length(r1$x))
  expect_true(all(r1$p_fw[r1$testable] > 0 & r1$p_fw[r1$testable] <= 1))
  ## reusing the noise tables gives the same result without new simulations
  r3 <- win_timeOpt(dat, window_size = 4, sedmin = 0.5, sedmax = 2, numsed = 6, fit = 2,
      bias_correction = "null", n_null = 8, step = 4, min_cycles = 0, null_cache = r1$null_cache)
  expect_identical(r1$z_opt, r3$z_opt)
})

test_that("time windows have the requested duration and match a direct timeOpt run", {
  skip_on_cran()
  set.seed(5)
  dx <- 0.05
  depth <- seq(0, 40, by = dx)
  x <- cumsum(rnorm(length(depth)))
  dat <- cbind(depth, x)
  res <- win_timeOpt(dat, sedmin = 0.5, sedmax = 2, numsed = 4, fit = 1, window_mode = "time",
      window_duration = 1000, step = 100, normalize = FALSE, bias_correction = "none")
  k <- which(res$x == 20)
  for (j in seq_along(res$sedrate)) {
    m <- res$npts[j, k]
    expect_equal((m - 1) * dx * 100 / res$sedrate[j], 1000, tolerance = 0.01)
    c0 <- which(abs(depth - 20) < 1e-9)
    rows <- (c0 - floor((m - 1) / 2)):(c0 - floor((m - 1) / 2) + m - 1)
    ref <- astrochron::timeOpt(dat[rows, ], sedmin = res$sedrate[j], sedmax = res$sedrate[j],
        numsed = 1, fit = 1, roll = 10^6, targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
        targetP = c(20.9, 19.9, 17.1, 17.2), output = 1, genplot = FALSE, verbose = FALSE,
        limit = FALSE)
    expect_equal(res$r_2_opt[j, k], ref[1, 4])
  }
  ## windows that do not fit in the record are untestable (code 3)
  expect_true(all(res$untestable_code[, 1] == 3L))
  expect_true(all(is.na(res$r_2_opt[, 1])))
})

test_that("the testable mask flags too few samples per carrier cycle and too few cycles", {
  set.seed(6)
  depth <- seq(0, 10, by = 0.05)
  dat <- cbind(depth, rnorm(length(depth)))
  res <- win_timeOpt(dat, window_size = 10, sedmin = 0.5, sedmax = 4, numsed = 5, fit = 1,
      bias_correction = "none", step = 50)
  few <- min(c(20.9, 19.9, 17.1, 17.2)) * res$sedrate / 100 < 3 * 0.05
  for (j in seq_along(res$sedrate)) {
    if (few[j]) {
      expect_true(all(res$untestable_code[j, ] == 1L))
    }
  }
  ## at 4 cm/kyr a 10 m window spans 250 kyr < 2 x 405 kyr
  expect_true(all(res$untestable_code[5, ] == 2L))
  expect_true(all(!is.finite(res$best_sedrate) | res$best_sedrate < 4))
})

test_that("the noise model fit recognises AR(1) and power-law noise", {
  set.seed(7)
  n <- 2000
  fa <- WaverideR:::wr_noise_fit(WaverideR:::wr_ar1(n, 0.8), seq_len(n))
  fp <- WaverideR:::wr_noise_fit(WaverideR:::wr_powerlaw(n, 1), seq_len(n))
  expect_equal(fa$rho, 0.8, tolerance = 0.05)
  expect_gt(fa$ll_ar1, fa$ll_pl)
  expect_equal(fp$beta, 1, tolerance = 0.15)
  expect_gt(fp$ll_pl, fp$ll_ar1)
})

test_that("win_timeOpt_multi compares the fits and the plots leave par unchanged", {
  skip_on_cran()
  set.seed(8)
  depth <- seq(0, 20, by = 0.05)
  dat <- cbind(depth, cumsum(rnorm(length(depth))))
  mr <- win_timeOpt_multi(dat, fits = c(1, 2, 3), window_size = 20, sedmin = 1, sedmax = 3,
      numsed = 5, n_null = 6, step = 100, min_cycles = 0)
  expect_s3_class(mr, "win_timeOpt_multi.result")
  expect_true(all(c("best_fit1", "best_fit2", "best_fit3", "agree_13", "ambiguous_13",
      "ambiguous_23") %in% names(mr$summary)))
  expect_type(mr$summary$ambiguous_13, "logical")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p0 <- graphics::par(no.readonly = TRUE)
  plot_win_timeOpt_multi(mr)
  plot_win_timeOpt(mr$fit1)
  plot_win_timeOpt(mr$fit2, corrected = FALSE)
  expect_identical(p0, graphics::par(no.readonly = TRUE))
})
