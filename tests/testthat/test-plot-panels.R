# the scalogram plots share one set of panels (R/plot_wavelet_panels.R); every layout must draw without error
draw <- function(expr) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 600, height = 500)
  on.exit({ grDevices::dev.off(); unlink(f) })
  force(expr)
  TRUE
}

set.seed(1)
z <- seq(0, 30, by = 0.25)
d <- data.frame(z, sin(2 * pi * z / 4) + as.numeric(arima.sim(list(ar = 0.4), length(z))))
wt <- analyse_wavelet(d, dj = 1 / 6, lowerPeriod = 0.8, upperPeriod = 16, verbose = FALSE, omega_nr = 6)

test_that("plot_wavelet draws all 16 layouts", {
  for (h in c(TRUE, FALSE)) for (dd in c(TRUE, FALSE)) for (av in c(TRUE, FALSE)) for (mt in c(TRUE, FALSE)) {
    expect_true(draw(plot_wavelet(wt, dev_new = FALSE, plot_horizontal = h, add_data = dd, add_avg = av,
                                  add_MTM = mt, add_MTM_peaks = mt, add_abline_h = 3, add_abline_v = 10,
                                  add_lines = cbind(z, rep(2, length(z))), keep_editable = TRUE)),
                info = paste(h, dd, av, mt))
  }
})

test_that("plot_wavelet draws p-values and significance contours", {
  wtp <- analyse_wavelet(d, dj = 1 / 6, lowerPeriod = 0.8, upperPeriod = 16, verbose = FALSE, omega_nr = 6,
                         pval = TRUE, n_simulations = 3)
  expect_equal(dim(wtp$Power.pval), dim(wtp$Power))
  sg <- wavelet_significance(wt)
  for (h in c(TRUE, FALSE)) {
    expect_true(draw(plot_wavelet(wtp, dev_new = FALSE, plot_horizontal = h, add_data = TRUE, add_avg = TRUE,
                                  add_pval = TRUE, add_signif = sg, keep_editable = TRUE)))
  }
})

test_that("all palette families work", {
  for (cb in list(c("grDevices", "rainbow"), c("grDevices", "heat.colors"), c("grDevices", "Viridis"),
                  c("viridis", "viridis"), c("viridis", "magma"), c("colorRamps", "matlab.like"),
                  c("RColorBrewer", "Spectral"))) {
    cols <- WaverideR:::.wr_key_cols(cb[2], cb[1], 20)
    expect_length(cols, 20)
    expect_true(draw(plot_wavelet(wt, dev_new = FALSE, palette_name = cb[2], color_brewer = cb[1],
                                  keep_editable = TRUE)), info = paste(cb, collapse = "/"))
  }
  expect_error(WaverideR:::.wr_key_cols("x", "nonsense", 5))
})

test_that("superlet, Awavelet, eha and cross spectrum plots share the panels", {
  sl <- analyse_superlet(d, upperPeriod = 16, lowerPeriod = 0.8, Nf = 20, c1 = 3, o = c(1, 3))
  aw <- analyse_Awavelet(d, dj = 1 / 6, lowerPeriod = 0.8, upperPeriod = 16, verbose = FALSE, omega_min = 6,
                         omega_max = 8, n_simulations = 2)
  d2 <- data.frame(z, sin(2 * pi * z / 4 + 1) + rnorm(length(z)))
  xw <- analyse_Xwavelet(d, d2, dj = 1 / 6, lowerPeriod = 0.8, upperPeriod = 16, verbose = FALSE, omega_nr = 6)
  for (h in c(TRUE, FALSE)) for (dd in c(TRUE, FALSE)) for (av in c(TRUE, FALSE)) {
    expect_true(draw(plot_superlet(sl, dev_new = FALSE, plot_horizontal = h, add_data = dd, add_avg = av,
                                   keep_editable = TRUE)))
    expect_true(draw(plot_Awavelet(aw, dev_new = FALSE, plot_horizontal = h, add_data = dd, add_avg = av,
                                   keep_editable = TRUE)))
    expect_true(draw(plot_Xwavelet(xw, dev_new = FALSE, plot_horizontal = h, add_data = dd, add_avg = av,
                                   n_arrows = 15, keep_editable = TRUE)))
  }
})
