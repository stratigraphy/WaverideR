test_that("the common colour breaks are increasing and pick log for sparse maps", {
  M <- matrix(0, 50, 50)
  M[25, ] <- 10^seq(-6, 2, length.out = 50)
  M[10, 1:20] <- 1e-3
  br <- .wr_colour_breaks(M, "auto")
  expect_identical(attr(br, "scale"), "log")
  expect_true(all(diff(br) > 0))
  expect_lte(br[1], min(M))
  expect_gte(br[length(br)], max(M))
  D <- matrix(runif(2500, 1, 2), 50)
  expect_identical(attr(.wr_colour_breaks(D, "auto"), "scale"), "quantile")
  for (sc in c("log", "quantile", "sqrt", "linear")) {
    expect_true(all(diff(.wr_colour_breaks(D, sc)) > 0))
  }
  N <- D
  N[1:5, ] <- NA
  expect_false(anyNA(.wr_fill_na(N)))
})

test_that("extract_amplitude is calibrated in both modes", {
  t <- seq(0, 2000, 1)
  a <- extract_amplitude(cbind(t, 3 * cos(2 * pi * t / 100)))
  expect_equal(stats::median(a[200:1800, 2]), 3, tolerance = 0.02)
  env <- 2 + sin(2 * pi * t / 400)
  x <- env * cos(2 * pi * t / 41) + 0.5 * cos(2 * pi * t / 100)
  a2 <- extract_amplitude(cbind(t, x), cycle = 41, band = c(0.8, 1.25))
  expect_gt(cor(a2[200:1800, 2], env[200:1800]), 0.98)
  expect_equal(stats::median(a2[200:1800, 2] / env[200:1800]), 1, tolerance = 0.05)
  expect_message(extract_amplitude(cbind(t, x), pts = 3), "no longer used")
})

test_that("the superlet cone of influence is monotonic in the period", {
  d <- mag[mag[, 1] <= 20, 1:2]
  s <- analyse_superlet(d, lowerPeriod = 0.2, upperPeriod = 8, Nf = 40, c1 = 3, o = c(1, 10), order_scaling = "log2", verbose = FALSE)
  expect_true(all(diff(s$coi[order(s$Period)]) >= 0))
})

test_that("plot_mesa_evolutive draws on a log axis with the common scale", {
  d <- mag[mag[, 1] <= 20, 1:2]
  me <- mesa_evolutive(d, window = 4, step = 0.5, fmax = 5)
  f <- tempfile(fileext = ".png")
  png(f)
  expect_silent(plot_mesa_evolutive(me))
  expect_silent(plot_mesa_evolutive(me, xaxis = "frequency", normalise = "period"))
  expect_silent(plot_mesa_evolutive(me, log_axis = FALSE, scale = "linear"))
  dev.off()
})
