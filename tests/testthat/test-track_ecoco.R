make_ecoco <- function() {
  sr <- seq(1, 10, by = 0.25)
  depth <- seq(0, 20, by = 0.5)
  true_sr <- 4 + 2 * depth / 20
  rho <- 0.8 * exp(-outer(sr, true_sr, "-")^2 / 0.5)
  list(sr = sr, depth = depth, rho = rho, pcoco = 2 * rho,
      p_global = matrix(0.01, length(sr), length(depth)),
      norbit = matrix(9, length(sr), length(depth)))
}

test_that("track_ecoco snaps scripted points to the ridge", {
  e <- make_ecoco()
  trk <- track_ecoco(e, points = data.frame(depth = c(2, 10, 18), sr = c(4.3, 5.1, 5.6)))
  expect_equal(nrow(trk), 3)
  expect_equal(trk$depth, c(2, 10, 18))
  expect_equal(trk$sedrate, c(4.25, 5, 5.75), tolerance = 0.3)
  expect_true(all(c("rho", "p_global", "norbit", "pcoco") %in% names(trk)))
})

test_that("track_ecoco interpolates and respects criteria", {
  e <- make_ecoco()
  trk <- track_ecoco(e, points = cbind(c(2, 18), c(4.25, 5.75)), interpolate = TRUE)
  expect_equal(nrow(trk), length(e$depth))
  expect_true(all(diff(trk$sedrate) >= 0))
  expect_error(track_ecoco(e, rho_min = 5, points = cbind(1, 1)), "no selectable")
  expect_error(track_ecoco(e), "mouse clicks")
})

test_that("track_ecoco can plot", {
  e <- make_ecoco()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(track_ecoco(e, points = cbind(10, 5), genplot = TRUE, plot_horizontal = FALSE), NA)
})
