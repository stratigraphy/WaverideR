make_records <- function() {
  t <- seq(0, 400, by = 1)
  astro <- data.frame(t, sin(2 * pi * t / 20) + 0.4 * sin(2 * pi * t / 41))
  d <- seq(0, 40, by = 0.1)
  # proxy: the same signal, sedimentation rate 0.1 m/kyr with a slowdown
  tt <- ifelse(d < 20, d / 0.1, 200 + (d - 20) / 0.05)
  proxy <- data.frame(d, sin(2 * pi * tt / 20) + 0.4 * sin(2 * pi * tt / 41))
  list(astro = astro[astro[, 1] <= 400, ], proxy = proxy)
}

test_that("clicking ties, unties and selects peaks in any order", {
  r <- make_records()
  pk_p <- as.data.frame(max_detect(r$proxy, pts = 3)[, 1:2]); names(pk_p) <- c("x", "y")
  pk_a <- as.data.frame(max_detect(r$astro, pts = 3)[, 1:2]); names(pk_a) <- c("x", "y")
  st <- list(ties = data.frame(p = integer(0), a = integer(0)), pending = NULL, history = list())
  st <- WaverideR:::.aai_click(st, "proxy", pk_p$x[3], pk_p, pk_a)
  expect_equal(st$pending$panel, "proxy")
  st <- WaverideR:::.aai_click(st, "astro", pk_a$x[5], pk_p, pk_a)
  expect_equal(nrow(st$ties), 1)
  # astronomical peak first, then the proxy peak
  st <- WaverideR:::.aai_click(st, "astro", pk_a$x[8], pk_p, pk_a)
  st <- WaverideR:::.aai_click(st, "proxy", pk_p$x[6], pk_p, pk_a)
  expect_equal(nrow(st$ties), 2)
  # clicking a tied peak (either record) unties it; ties can be removed in random order
  st <- WaverideR:::.aai_click(st, "astro", pk_a$x[5], pk_p, pk_a)
  expect_equal(st$ties$a, 8)
  st <- WaverideR:::.aai_click(st, "proxy", pk_p$x[6], pk_p, pk_a)
  expect_equal(nrow(st$ties), 0)
})

test_that("tie points in the wrong order are flagged", {
  pk_p <- data.frame(x = 1:5, y = 0)
  pk_a <- data.frame(x = c(10, 20, 15, 40, 50), y = 0)
  ties <- data.frame(p = 1:5, a = 1:5)
  fl <- WaverideR:::.aai_flag(ties, pk_p, pk_a, TRUE)
  expect_equal(fl, c(FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(WaverideR:::.aai_flag(ties, pk_p, pk_a, FALSE)[1], TRUE)
})

test_that("the window can be driven with scripted clicks and returns the tie points", {
  r <- make_records()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  pk_p <- max_detect(r$proxy, pts = 3)[, 1:2]
  pk_a <- max_detect(r$astro, pts = 3)[, 1:2]
  at <- function(geo, name, x) {
    g <- geo[[name]]
    f <- (x - g[5]) / (g[6] - g[5])
    if (g[9] == 1) f <- 1 - f
    list(x = g[1] + f * (g[2] - g[1]), y = mean(g[3:4]))
  }
  btn <- function(geo, nm) {
    b <- geo$buttons[[nm]]
    list(x = mean(b[1:2]), y = mean(b[3:4]))
  }
  # tie the first 4 proxy peaks to astronomical peaks, remove the second tie, tune, finish
  steps <- list(
    function(g) at(g, "proxy", pk_p[2, 1]), function(g) at(g, "astro", pk_a[2, 1]),
    function(g) at(g, "astro", pk_a[3, 1]), function(g) at(g, "proxy", pk_p[3, 1]),
    function(g) at(g, "proxy", pk_p[4, 1]), function(g) at(g, "astro", pk_a[4, 1]),
    function(g) at(g, "proxy", pk_p[3, 1]),
    function(g) btn(g, "Zoom -"), function(g) btn(g, "P >"), function(g) btn(g, "Fit"),
    function(g) btn(g, "Done"))
  k <- 0
  loc <- function() {
    k <<- k + 1
    if (k > length(steps)) return(NULL)
    steps[[k]](WaverideR:::.aai_env$geo)
  }
  out <- astro_anchor_interactive(r$astro, r$proxy, locator_fun = loc)
  expect_equal(names(out), c("data_x", "data_y", "insolation_x", "insolation_y"))
  expect_equal(nrow(out), 2)
  expect_equal(out$data_x, pk_p[c(2, 4), 1])
  expect_equal(out$insolation_x, pk_a[c(2, 4), 1])
})

test_that("tie points can be given, are snapped and tuning works", {
  r <- make_records()
  ties <- data.frame(data_x = c(2.1, 10.0, 25), insolation_x = c(21, 101, 301))
  out <- astro_anchor_interactive(r$astro, r$proxy, ties = ties)
  expect_equal(nrow(out), 3)
  expect_true(all(abs(out$data_x - ties$data_x) < 1.5))
  expect_error(astro_anchor_interactive(r$astro, r$proxy), "interactive")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  pk_a <- as.data.frame(max_detect(r$astro, pts = 3)[, 1:2]); names(pk_a) <- c("x", "y")
  full <- astro_anchor_interactive(r$astro, r$proxy, ties = data.frame(
    data_x = c(2, 6, 10, 14, 18, 24, 30, 36), insolation_x = c(20, 60, 100, 140, 180, 280, 380, 400)))
  astro_xy <- setNames(r$astro, c("x", "y"))
  res <- WaverideR:::.aai_inspect(r$proxy, astro_xy, pk_a, full, TRUE, NULL, 0.25, 100, FALSE)
  expect_true(is.finite(res$correlation))
  expect_true(all(res$sedrate > 0))
  bad <- full[c(1, 3, 2, 4:8), ]
  bad$insolation_x[c(2, 3)] <- bad$insolation_x[c(3, 2)]
  bad <- bad[order(bad$data_x), ]
  bad$insolation_x[2] <- bad$insolation_x[3] + 5
  expect_error(WaverideR:::.aai_inspect(r$proxy, astro_xy, pk_a, bad, TRUE, NULL, 0.25, 100, FALSE), "wrong order")
})

test_that("sync moves the window of the other record to the partners of the visible tie points", {
  pk_p <- data.frame(x = c(1, 2, 3, 40, 41), y = 0)
  pk_a <- data.frame(x = c(100, 120, 140, 600, 620), y = 0)
  proxy <- data.frame(x = 0:50, y = 0)
  astro <- data.frame(x = seq(0, 700, by = 10), y = 0)
  st <- list(ties = data.frame(p = c(4, 5), a = c(4, 5)), pending = NULL, history = list(),
             view_p = c(38, 43), view_a = c(0, 200))
  st2 <- WaverideR:::.aai_sync(st, pk_p, pk_a, proxy, astro)
  expect_true(all(c(600, 620) >= st2$view_a[1] & c(600, 620) <= st2$view_a[2]))
  expect_equal(diff(st2$view_a), 200)
})
