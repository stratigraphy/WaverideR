## tests for the date-model options of astroMCBchron (version 2): outlier mixture, laboratory groups, depth groups,
## systematic envelopes and single-zircon likelihoods (small chains; checks behaviour, not exact values)
test_that("astroMCBchron version 2 date models behave as expected", {
  skip_on_cran()
  set.seed(1)
  rate <- function(z) 0.02 * (1 + 0.4 * sin(2 * pi * z / 40))
  zz <- seq(0, 100, 0.05)
  age_z <- 400000 + c(0, cumsum(diff(zz) / rate(zz[-1]))) + ifelse(zz > 60, 300, 0)
  true_age <- function(z) approx(zz, age_z, z)$y / 1000
  trk <- function(z) {
    P <- matrix(NA_real_, length(z), 50)
    for (i in 1:50) {
      P[, i] <- 405 * rate(z) * exp(rnorm(1, 0, 0.02))
    }
    list(position = z, period = P)
  }
  TR <- list(trk(seq(0, 60, 0.5)), trk(seq(60, 100, 0.5)))
  HI <- data.frame(position = 60, mean = 300, sd = 40)
  zd <- c(8, 22, 37, 52, 68, 83, 96)
  d <- data.frame(id = paste0("ash", 1:7), position = zd, age = true_age(zd) + rnorm(7, 0, 0.04), sd = 0.04, thickness = 0.005)
  d$age[4] <- d$age[4] + 0.2                                                       # one date 5 sd too old
  run <- function(x, ...) astroMCBchron(x, TR, hiatuses = HI, n_iter = 3000, thin = 10, n_tracks = 50, loo = FALSE,
      predict_positions = seq(0, 100, 1), verbose = FALSE, seed = 2, ...)
  f <- run(d, outlier_type = "old")
  expect_true(f$dates$p_outlier[4] > 0.5)
  expect_true(all(f$dates$p_outlier[-4] < 0.5))
  expect_true(all(f$summary$sd2_total_kyr > f$summary$sd2_tracer_kyr & f$summary$sd2_tracer_kyr > f$summary$sd2_kyr))
  f2 <- run(d, date_model = "two_stage")
  expect_equal(f2$settings$date_model, "two_stage")
  ## laboratory group with a +150 kyr bias
  d3 <- d
  d3$age[4] <- d3$age[4] - 0.2
  d3$group <- c("A", "B", "A", "A", "B", "A", "A")
  d3$tracer_rel <- ifelse(d3$group == "B", 5e-4, 1.5e-4)
  d3$age[d3$group == "B"] <- d3$age[d3$group == "B"] + 0.15
  f3 <- run(d3)
  expect_true(median(f3$group_offsets[, "B"]) > 50)
  ## dates on a depth scale 1.5 m too deep
  d4 <- d
  d4$age[4] <- d4$age[4] - 0.2
  d4$position[3:5] <- d4$position[3:5] + 1.5
  d4$depth_group <- c(NA, NA, "box", "box", "box", NA, NA)
  f4 <- run(d4, depth_offset_sd = 2)
  expect_true(median(f4$depth_offsets) < 0)
  ## single-zircon likelihood
  g <- data.frame(id = rep(d$id[1:2], each = 6), age = rep(true_age(zd[1:2]), each = 6) + runif(12, 0, 0.1) + rnorm(12, 0, 0.05), sd = 0.05)
  gl <- upb_grain_loglik(g)
  expect_equal(nrow(gl$summary), 2)
  expect_true(all(is.finite(gl$loglik[[1]](gl$summary$te_mean[1] + c(-0.1, 0, 0.1)))))
})
