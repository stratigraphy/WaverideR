test_that("the TimeOpt part of astro_sedrate_test reproduces astrochron::timeOpt", {
  skip_on_cran()
  set.seed(1)
  depth <- seq(0, 60, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(depth)), 0.8, method = "recursive"))
  dat <- cbind(depth, x)
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(NULL)
  y <- residuals(lm(x ~ depth))
  y <- (y - mean(y)) / sd(y)
  for (fit in 1:2) {
    cfg <- ns$wr_ast_cfg(tg, fit, FALSE, TRUE, NULL, 10^7, "BIC", 3, FALSE)
    o <- ns$wr_ast_eval(y, 0.1, c(2, 4), cfg)
    for (j in 1:2) {
      sr <- c(2, 4)[j]
      ref <- astrochron::timeOpt(dat, sedmin = sr, sedmax = sr, numsed = 1, fit = fit,
          flow = cfg$flow, fhigh = cfg$fhigh, roll = 10^7, targetE = tg$E, targetP = tg$P,
          output = 1, genplot = FALSE, verbose = FALSE, limit = FALSE)
      expect_equal(o$r2_env[j], ref[1, 2], tolerance = 1e-8)
      expect_equal(o$r2_pwr[j], ref[1, 3], tolerance = 1e-8)
    }
  }
})

test_that("age-dependent targets use the AstroGeo22 precession frequency", {
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(418)
  a <- astrogeo22(418)
  expect_equal(tg$p_arcsec, a$p_arcsec)
  expect_equal(max(tg$E), 405.6912, tolerance = 1e-3)
  expect_true(all(tg$P < 21) && all(tg$P > 15))
  tg0 <- ns$wr_ast_targets(NULL)
  expect_equal(round(tg0$P, 1), c(23.7, 22.4, 19.1, 19.0))
})

test_that("unresolved target periods are merged and too short or long ones dropped", {
  ns <- asNamespace("WaverideR")
  p <- c(405.7, 130.7, 123.8, 98.9, 94.9)
  expect_length(ns$wr_ast_resolvable(p, dur = 3000, dt = 2, min_samples = 3), 5)
  r <- ns$wr_ast_resolvable(p, dur = 1000, dt = 2, min_samples = 3)
  expect_length(r, 3)
  expect_length(ns$wr_ast_resolvable(p, dur = 300, dt = 2, min_samples = 3), 1)
})

test_that("noise models are recognised", {
  set.seed(3)
  ns <- asNamespace("WaverideR")
  fa <- ns$wr_ast_noise_fit(ns$wr_ar1(3000, 0.8))
  expect_true(fa$model[which.min(fa$AIC)] %in% c("ar1", "ar2"))
  fp <- ns$wr_ast_noise_fit(ns$wr_ast_noise_series(3000, "powerlaw", 1))
  expect_true(fp$model[which.min(fp$AIC)] %in% c("powerlaw", "bending"))
})

test_that("the band mask removes the astronomical bands and the local fit keeps enough frequencies", {
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(NULL)
  n <- 2000
  dt <- 2
  keep <- ns$wr_ast_band_mask(n, dt, tg, pad = 2)
  fk <- (1:floor((n - 1) / 2)) / (n * dt)
  expect_false(any(keep[abs(fk - 1 / 405) < 1e-4]))
  expect_false(any(keep[fk > 1 / 24 & fk < 1 / 19]))
  expect_true(mean(keep) > 0.5)
  kl <- ns$wr_ast_band_mask(n, dt, tg, pad = 2, fmax = 4 / 95)
  expect_true(sum(kl) >= 40)
  expect_true(all(fk[kl] <= max(4 / 95, fk[kl][40])))
})

test_that("noise fit uncertainty is returned and draws stay valid", {
  set.seed(5)
  ns <- asNamespace("WaverideR")
  x <- ns$wr_ast_noise_series(1500, "powerlaw", 1)
  f <- ns$wr_ast_noise_fit(x, rep(c(TRUE, FALSE), length.out = 749))
  expect_true(is.matrix(f$theta_cov[["powerlaw"]]))
  p <- ns$wr_ast_draw_par("ar2", f$param[["ar2"]], f$theta_cov[["ar2"]], c(3, -3))
  expect_true(abs(p[2]) < 1 && abs(p[1]) < 2)
})

test_that("minimal tuning recovers a variable-rate 405-kyr cycle", {
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(NULL)
  cfg <- ns$wr_ast_cfg(tg, 2, FALSE, TRUE, NULL, 10^7, "BIC", 3, FALSE)
  cfg <- ns$wr_ast_add_tune(cfg, tg, "minimal", 0.3, FALSE)
  expect_false(any(abs(cfg$spec_tuned - tg$E[1]) < 1e-6))
  dx <- 0.2
  z <- seq(0, 250, by = dx)
  rate <- 6 + 3 * z / max(z)
  tt <- c(0, cumsum(diff(z) * 100 / ((rate[-1] + rate[-length(z)]) / 2)))
  y <- cos(2 * pi * tt / tg$E[1])
  tu <- ns$wr_ast_min_tune(y, dx, mean(rate), cfg)
  expect_false(is.null(tu))
  dur_true <- max(tt)
  expect_lt(abs(max(tu$model$time) - dur_true) / dur_true, 0.1)
  expect_true(all(diff(tu$model$time) > 0))
})

test_that("masked and local noise fits run and report per-rate fits", {
  skip_on_cran()
  set.seed(6)
  depth <- seq(0, 40, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive"))
  for (nf in c("masked", "local")) {
    r <- astro_sedrate_test(cbind(depth, x), sedmin = 1, sedmax = 4, numsed = 8, fit = 1,
        null_fit = nf, n_null = 6, sensitivity = FALSE, target_sensitivity = FALSE, coco = FALSE)
    expect_equal(r$settings$null_fit, nf)
    expect_true(all(c("sedrate", "model", "AIC", "kept") %in% names(r$noise_fits_rate)))
    ok <- r$scores$untestable == 0
    expect_true(all(r$scores$p_opt[ok] > 0 & r$scores$p_opt[ok] <= 1))
  }
})

test_that("astro_sedrate_test returns scores, probabilities and leaves par unchanged", {
  skip_on_cran()
  set.seed(4)
  depth <- seq(0, 40, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive"))
  r <- astro_sedrate_test(cbind(depth, x), sedmin = 1, sedmax = 4, numsed = 12, fit = 1,
      n_null = 8, sensitivity = TRUE, target_sensitivity = FALSE)
  expect_s3_class(r, "astro_sedrate_test")
  expect_equal(nrow(r$scores), 12)
  ok <- r$scores$untestable == 0
  expect_true(all(r$scores$p_comb[ok] > 0 & r$scores$p_comb[ok] <= 1))
  expect_true(all(is.na(r$scores$z_opt[!ok])))
  expect_equal(sum(r$scores$weight_lik, na.rm = TRUE), 1, tolerance = 1e-8)
  expect_equal(length(unique(r$sensitivity$noise_model)), 4)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p0 <- graphics::par(no.readonly = TRUE)
  plot(r)
  expect_identical(p0, graphics::par(no.readonly = TRUE))
  expect_output(print(r), "best rates")
})

test_that("fit-2 band guard, obliquity-modulation terms and edge trim are set as documented", {
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(418)
  expect_equal(round(tg$Omod), c(173, 151))
  c0 <- ns$wr_ast_cfg(tg, 2, FALSE, TRUE, NULL, 10^7, "BIC", 3, FALSE)
  c1 <- ns$wr_ast_cfg(tg, 2, FALSE, TRUE, NULL, 10^7, "BIC", 3, FALSE, TRUE, TRUE, "auto")
  expect_equal(c0$trim, 0)
  expect_lt(1 / c0$flow, 152)
  expect_gt(1 / c0$flow, 150)
  expect_lt(1 / c1$flow, 141)
  expect_gt(1 / c1$flow, 139)
  expect_true(all(round(tg$Omod, 3) %in% round(c1$spec, 3)))
  expect_equal(c1$trim, 1 / (c1$fhigh - c1$flow))
})

test_that("the significance spectrum is calibrated per rate and resolves small p", {
  set.seed(9)
  ns <- asNamespace("WaverideR")
  nr <- 30
  nsim <- 200
  ## noise whose spread grows with the rate index (like r2 under noise)
  V <- matrix(rnorm(nr * nsim), nr, nsim) * rep(seq(0.5, 3, length.out = nr), nsim)
  obs <- rnorm(nr) * seq(0.5, 3, length.out = nr)
  obs[8] <- 5 * (0.5 + 7 * 2.5 / 29)
  sg <- ns$wr_ast_sigspec(obs, V, rep(TRUE, nr), 2)
  expect_equal(which.min(sg$p), 8)
  expect_lt(sg$fw[8], 0.05)
  expect_lt(min(sg$p), 1 / nsim)
  expect_true(all(sg$fw > 0 & sg$fw <= 1))
  expect_gt(sg$levels[["p01"]], sg$levels[["p05"]])
  ## pure noise: the family-wise p of each noise series is roughly uniform
  fwk <- numeric(40)
  for (k in 1:40) {
    s2 <- ns$wr_ast_sigspec(V[, k], V[, -k], rep(TRUE, nr), 2)
    fwk[k] <- min(s2$fw)
  }
  expect_lt(mean(fwk < 0.05), 0.2)
})

test_that("astro_sedrate_test returns the significance spectrum and its levels", {
  skip_on_cran()
  set.seed(6)
  depth <- seq(0, 40, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive"))
  r <- astro_sedrate_test(cbind(depth, x), sedmin = 1, sedmax = 4, numsed = 8, fit = 1,
      n_null = 10, sensitivity = FALSE, target_sensitivity = FALSE, coco = FALSE)
  ok <- r$scores$untestable == 0
  expect_true(all(c("S_opt", "p_S") %in% names(r$scores)))
  expect_true(all(r$scores$p_S[ok] > 0 & r$scores$p_S[ok] <= 1))
  expect_true(all(is.na(r$scores$S_opt[!ok])))
  expect_equal(dim(r$significance_levels), c(3L, 2L))
  expect_equal(r$best$score[2], "TimeOpt significance S_opt")
})

test_that("tuned-band peak helpers pick the untuned bands and find a spectral peak", {
  ns <- asNamespace("WaverideR")
  tg <- ns$wr_ast_targets(NULL)
  f2 <- ns$wr_ast_peak_families(tg, 2)
  f1 <- ns$wr_ast_peak_families(tg, 1)
  expect_length(f2, 5)
  expect_length(f1, 3)
  expect_false(any(abs(f2 - tg$E[1]) < 1))
  expect_true(any(abs(f1 - 40.3) < 0.5))
  set.seed(2)
  tt <- 0:3999
  noise <- as.numeric(stats::filter(rnorm(4000), 0.7, method = "recursive"))
  v1 <- noise + 1.5 * cos(2 * pi * tt / 40.3) + 1.5 * cos(2 * pi * tt / 23)
  expect_gt(ns$wr_ast_comb(v1, 1, c(40.3, 23), 0.05), ns$wr_ast_comb(noise, 1, c(40.3, 23), 0.05) + 1)
})

test_that("astro_sedrate_test returns tuned-band peak and joint scores", {
  skip_on_cran()
  set.seed(6)
  depth <- seq(0, 60, by = 0.1)
  x <- as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive"))
  r <- astro_sedrate_test(cbind(depth, x), sedmin = 1, sedmax = 3, numsed = 6, fit = 1,
      n_null = 10, sensitivity = FALSE, target_sensitivity = FALSE, coco = FALSE)
  expect_true(all(c("S_peaks", "p_peaks", "S_joint", "p_joint") %in% names(r$scores)))
  ok <- is.finite(r$scores$S_joint)
  expect_true(any(ok))
  expect_true(all(r$scores$p_joint[ok] > 0 & r$scores$p_joint[ok] <= 1))
  expect_equal(rownames(r$significance_levels), c("TimeOpt", "peaks", "joint"))
  r0 <- astro_sedrate_test(cbind(depth, x), sedmin = 1, sedmax = 3, numsed = 6, fit = 1,
      n_null = 5, sensitivity = FALSE, target_sensitivity = FALSE, coco = FALSE, tuned_peaks = FALSE)
  expect_true(all(is.na(r0$scores$S_peaks)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(r))
})

test_that("astro_sedrate_multi combines series of one section and sections", {
  skip_on_cran()
  set.seed(8)
  depth <- seq(0, 60, by = 0.1)
  a <- cbind(depth, as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive")))
  b <- cbind(depth, a[, 2] + rnorm(length(depth)))
  cc <- cbind(depth, as.numeric(stats::filter(rnorm(length(depth)), 0.7, method = "recursive")))
  r <- astro_sedrate_multi(list(a = a, b = b, c = cc), section = c("s1", "s1", "s2"),
      sedmin = c(s1 = 1, s2 = 1.2), sedmax = 3, numsed = 6, fit = 1, n_null = 10, coco = FALSE)
  expect_s3_class(r, "astro_sedrate_multi")
  expect_equal(sort(names(r$sections)), c("s1", "s2"))
  expect_true(all(r$best$p_fw > 0 & r$best$p_fw <= 1))
  expect_true(is.finite(r$p_all))
  expect_true(all(r$sections$s1$n_components[is.finite(r$sections$s1$S_multi)] <= 4))
  r2 <- astro_sedrate_multi(list(a = a, b = b), sedmin = 1, sedmax = 3, numsed = 6, fit = 1,
      n_null = 10, coco = FALSE, combine = "stouffer", scores = "TimeOpt")
  expect_equal(nrow(r2$best), 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(r))
  expect_output(print(r), "astro_sedrate_multi")
})
