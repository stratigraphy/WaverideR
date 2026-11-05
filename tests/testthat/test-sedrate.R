test_that("model_red_noise_wt and sum_power_sedrate run and are consistent", {
  skip_on_cran()
  w <- analyse_wavelet(data = mag[1:150, ], dj = 1/8, lowerPeriod = 0.2, upperPeriod = 16,
      verbose = FALSE, omega_nr = 8)
  set.seed(2)
  rn <- model_red_noise_wt(w, n_simulations = 4)
  expect_equal(dim(rn), c(length(w$Period), 4))
  expect_true(all(is.finite(rn)) && all(rn >= 0))
  sr <- sum_power_sedrate(red_noise = rn, wavelet = w, percentile = 0.75, sedrate_low = 0.5,
      sedrate_high = 6, spacing = 1, cycles = c(2376, 1600, 1180, 696, 406, 110))
  expect_equal(dim(sr$power_max_mat), c(6, 150))
  expect_true(all(sr$nr_compenents_ma >= 0 & sr$nr_compenents_ma <= 6))
  expect_true(all(sr$power_max_mat[sr$nr_compenents_ma == 0] == 0))
})
