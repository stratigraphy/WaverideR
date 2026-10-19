# a short record with a smoothly oscillating sedimentation rate, moderate noise
mk <- function(seed = 1, noise = 0.8) {
  set.seed(seed)
  t <- seq(0, 2400, by = 1)
  x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41) + 0.4 * cos(2 * pi * t / 100)
  sr <- 2 + 0.7 * sin(2 * pi * t / 600)
  zt <- c(0, cumsum(sr[-1] * diff(t)))
  z <- seq(0, max(zt), by = 1)
  x <- stats::approx(zt, x_t, z)$y
  list(d = data.frame(z, x + rnorm(length(z), sd = noise * sd(x))), z = z, sr = stats::approx(zt, sr, z)$y)
}

test_that("chirp-matched EHA recovers the rate gradient and the rate with the template", {
  a <- mk()
  ap <- data.frame(period = c(100, 41, 23.7, 19), weight = c(0.3, 1, 1, 0.6))
  ce <- analyse_chirp_eha(a$d, window = 600, step = 10, chirp = seq(-0.004, 0.004, by = 0.001), nw = 2, lowerPeriod = 10,
                          upperPeriod = 300, dj = 1 / 40, criterion = "template", astro_periods = ap, sedrate = c(0.8, 5))
  expect_s3_class(ce, "analyze.wavelet")
  expect_equal(dim(ce$template_score_all)[3], length(ce$axis.1))
  st <- stats::approx(a$z, a$sr, ce$axis.1)$y
  gt <- c(NA, diff(log(st)) / diff(ce$axis.1))
  expect_gt(cor(ce$chirp, gt, use = "complete"), 0.3)
  tk <- track_chirp_hmm(ce, n_samples = 0)
  expect_s3_class(tk, "wr_track")
  expect_gt(cor(log(tk$sedrate), log(st)), 0.8)
  expect_gt(mean(abs(log(tk$sedrate / st)) < log(1.15)), 0.7)
  expect_true(!is.null(attr(tk, "posterior")))
  # steered map
  cs <- analyse_chirp_eha(a$d, window = 600, step = 10, chirp = seq(-0.004, 0.004, by = 0.001), nw = 2, lowerPeriod = 10,
                          upperPeriod = 300, dj = 1 / 40, steer = tk)
  expect_equal(dim(cs$Power_norm), dim(cs$Power))
})

test_that("chirp-matched wavelet evidence and the gradient HMM on a plain scalogram run", {
  a <- mk(2, 0.5)
  ap <- data.frame(period = c(100, 41, 23.7, 19), weight = c(0.3, 1, 1, 0.6))
  cw <- analyse_chirp_wavelet(a$d, lowerPeriod = 10, upperPeriod = 300, dj = 1 / 40, omega_nr = 8,
                              chirp = seq(-0.004, 0.004, by = 0.002), astro_periods = ap, sedrate = c(0.8, 5), step = 10)
  expect_equal(length(cw$positions), dim(cw$template_score_all)[3])
  tk <- track_chirp_hmm(cw, n_samples = 0)
  st <- stats::approx(a$z, a$sr, tk$depth)$y
  expect_gt(cor(log(tk$sedrate), log(st)), 0.5)
  w <- analyse_wavelet(a$d, lowerPeriod = 10, upperPeriod = 300, dj = 1 / 40, omega_nr = 8, verbose = FALSE)
  tk2 <- suppressMessages(track_chirp_hmm(w, astro_periods = ap, sedrate = c(0.8, 5), chirp = seq(-0.004, 0.004, by = 0.002), n_samples = 0))
  expect_s3_class(tk2, "wr_track")
  cws <- analyse_chirp_wavelet(a$d, lowerPeriod = 10, upperPeriod = 300, dj = 1 / 40, omega_nr = 8,
                               chirp = seq(-0.004, 0.004, by = 0.002), astro_periods = ap, sedrate = c(0.8, 5), step = 10,
                               steer = tk, squeeze = TRUE)
  expect_equal(dim(cws$Power_sst), dim(cws$Power))
})

test_that("automatic windows, jumps, combined evidence and the surrogate significance run", {
  a <- mk(3, 0.6)
  ap <- data.frame(period = c(100, 41, 23.7, 19), weight = c(0.3, 1, 1, 0.6))
  ce <- suppressMessages(analyse_chirp_eha(a$d, window = NULL, chirp = "auto", nw = 2, lowerPeriod = 10, upperPeriod = 300,
                                           dj = 1 / 40, criterion = "template", astro_periods = ap, sedrate = c(0.8, 5)))
  expect_length(ce$windows, 2)
  expect_equal(length(ce$positions), dim(ce$template_score_all)[3])
  w <- analyse_wavelet(a$d, lowerPeriod = 10, upperPeriod = 300, dj = 1 / 40, omega_nr = 8, verbose = FALSE)
  tk <- suppressMessages(track_chirp_hmm(ce, jump_prob = 1e-3, scalogram = w, scalogram_args = list(astro_periods = ap), n_samples = 0))
  st <- stats::approx(a$z, a$sr, tk$depth)$y
  expect_gt(cor(log(tk$sedrate), log(st)), 0.7)
  cs <- analyse_chirp_eha(a$d, window = 600, step = 20, chirp = seq(-0.004, 0.004, by = 0.001), nw = 2, lowerPeriod = 10,
                          upperPeriod = 300, dj = 1 / 40, steer = tk)
  cs <- chirp_eha_significance(cs, a$d, n_sim = 3)
  expect_equal(dim(cs$Power_sig), dim(cs$Power_norm))
  expect_true(all(cs$Power_pval >= 0 & cs$Power_pval <= 1))
})

test_that("track_astro_auto(method = 'both') returns both tracks and their agreement; plot and rate_gram work", {
  a <- mk(4, 0.6)
  ap <- data.frame(period = c(100, 41, 23.7, 19), weight = c(0.3, 1, 1, 0.6))
  res <- suppressMessages(suppressWarnings(track_astro_auto(a$d, sedrate = c(0.8, 5), astro_periods = ap, scalogram = "cwt",
                                                            omega_nr = 8, dj = 1 / 40, n_samples = 5)))
  expect_true(all(c("track", "track_hmm", "chirp_eha", "agreement") %in% names(res)))
  expect_true("disagreement" %in% names(res$check))
  st <- stats::approx(a$z, a$sr, res$track$depth)$y
  expect_gt(cor(log(res$track$sedrate), log(st)), 0.7)
  pdf(NULL)
  rg <- rate_gram(res$chirp_eha, k = 5)
  expect_equal(length(rg$best), length(rg$depth))
  expect_silent(plot_chirp_eha(res$chirp_eha, style = "plain"))
  dev.off()
})

test_that("analyse_dechirped maps a transform to the time of a track and back", {
  a <- mk(5, 0.5)
  ap <- data.frame(period = c(100, 41, 23.7, 19), weight = c(0.3, 1, 1, 0.6))
  w <- analyse_wavelet(a$d, lowerPeriod = 10, upperPeriod = 300, dj = 1 / 40, omega_nr = 8, verbose = FALSE)
  tk <- suppressMessages(track_astro_hmm(w, astro_periods = ap, sedrate = c(0.8, 5), n_samples = 5))
  dc <- analyse_dechirped(a$d, tk, "wavelet", lowerPeriod = 12, upperPeriod = 150, dj = 1 / 40, omega_nr = 8, n_sim = 2)
  expect_s3_class(dc, "wr_dechirped")
  expect_equal(ncol(dc$Power_norm), length(dc$depth$axis.1))
  expect_equal(dim(dc$Power_sig), dim(dc$Power_norm))
  # in the tuned scalogram the 23.7-kyr line is at 23.7 kyr
  j <- which.min(abs(dc$tuned$Period - 23.7)); band <- abs(log2(dc$tuned$Period / 23.7)) < 0.3
  expect_equal(dc$tuned$Period[band][which.max(rowMeans(dc$tuned$Power[band, ]))], 23.7, tolerance = 0.1)
  pdf(NULL); plot_dechirped(dc, "tuned", show = 23.7); plot_dechirped(dc, "outline", show = 23.7); dev.off()
})
