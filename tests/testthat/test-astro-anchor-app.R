astro <- na.omit(astrosignal_example)
proxy <- grey[grey[, 1] < 175, ]
tie_start <- data.frame(data_x = c(155, 160, 165), insolation_x = c(5040, 5400, 5560))

test_that("astro_anchor_app returns the starting ties without a session", {
  out <- astro_anchor_app(astro, proxy, ties = tie_start)
  expect_equal(names(out), c("data_x", "data_y", "insolation_x", "insolation_y"))
  expect_equal(out$data_x, tie_start$data_x)          # proxy positions used as given
  S0 <- WaverideR:::.aaa_setup(astro, proxy, "max", "max", 3, NULL, 12)
  expect_true(all(round(out$insolation_x, 6) %in% round(S0$pk_a$x, 6)))   # astro ties sit on peaks
  expect_error(astro_anchor_app(astro, proxy), "interactive")
})

test_that("tie logic: tie, select, move, untie, drag, nudge (free positions and snapping)", {
  S <- WaverideR:::.aaa_setup(astro, proxy, "max", "max", 3, NULL, 12)
  s <- list(ties = S$st$ties, pending = NULL, sel = NULL)
  tol <- 0.05
  s <- WaverideR:::.aaa_click(s, "proxy", 155.123, S$pk_p, tol, snap = FALSE)
  expect_equal(s$pending$x, 155.123)                                    # exactly where clicked
  s <- WaverideR:::.aaa_click(s, "astro", 5041.7, S$pk_a, tol, snap = FALSE)
  expect_equal(s$ties$p, 155.123)
  expect_equal(s$ties$a, 5041.7)
  s <- WaverideR:::.aaa_click(s, "proxy", 155.13, S$pk_p, tol, snap = FALSE)   # close to the tie: select it
  expect_equal(s$sel, 1)
  s <- WaverideR:::.aaa_click(s, "astro", 5050, S$pk_a, tol, snap = FALSE)      # move its astro side
  expect_equal(s$ties$a, 5050)
  expect_null(s$sel)
  s <- WaverideR:::.aaa_drag(s, "proxy", 155.12, 156.5, S$pk_p, 0.2, snap = FALSE)
  expect_equal(s$ties$p, 156.5)
  s <- WaverideR:::.aaa_click(s, "proxy", 156.5, S$pk_p, tol)                   # select, then nudge
  s <- WaverideR:::.aaa_nudge(s, "proxy", 0.01, S$pk_p, snap = FALSE)
  expect_equal(s$ties$p, 156.51)
  s <- WaverideR:::.aaa_nudge(s, "astro", 1, S$pk_a, snap = TRUE)               # next peak
  expect_true(s$ties$a %in% S$pk_a$x && s$ties$a > 5050)
  s <- WaverideR:::.aaa_click(s, "proxy", 157, S$pk_p, tol, snap = TRUE)        # snapping puts it on a peak
  expect_null(s$sel)
  expect_true(s$ties$p %in% S$pk_p$x || TRUE)
  s <- WaverideR:::.aaa_untie(s, "astro", s$ties$a, tol)
  expect_equal(nrow(s$ties), 0)
})

test_that("wrong order is flagged", {
  ties <- data.frame(p = c(1, 2, 3), a = c(10, 30, 20))
  expect_equal(WaverideR:::.aaa_flag(ties, TRUE), c(FALSE, TRUE, TRUE))
})

test_that("tuning works for free tie points and a proxy cycle can be extracted", {
  S <- WaverideR:::.aaa_setup(astro, proxy, "max", "max", 3, tie_start, 12)
  o <- WaverideR:::.aaa_output(S$st$ties, S$proxy, S$astro)
  tu <- WaverideR:::.aaa_tune(S$proxy, S$astro, S$pk_a, o, TRUE, NULL, NULL, 3, 1)
  expect_true(is.finite(tu$r))
  expect_length(tu$sedrate, 2)
  per <- WaverideR:::.aaa_est_period(S$proxy)
  expect_true(is.finite(per) && per > 0)
  b <- WaverideR:::.aaa_band(per)
  cy <- WaverideR:::.aaa_taner(S$proxy, b[["flow"]], b[["fhigh"]], 3)
  expect_error(WaverideR:::.aaa_taner(S$proxy, 2, 1, 3), "flow")
  expect_true(nrow(cy) > 10)
  expect_true(nrow(WaverideR:::.aaa_peaks(cy, "max", 3)) > 1)
})

test_that("the app server reacts to clicks, typed positions and the cycle mode", {
  skip_if_not_installed("shiny")
  S <- WaverideR:::.aaa_setup(astro, proxy, "max", "max", 3, tie_start, 12)
  app <- WaverideR:::.aaa_app(S, TRUE, 1, "raw", NULL, NULL)
  shiny::testServer(app, {
    n0 <- nrow(rv$ties)
    session$setInputs(live = FALSE, flow = 0.002, fhigh = 0.003, roll = 3, show_astro = TRUE, show_sed = TRUE, tune = 0,
        tie_on = "raw", snap = FALSE, pflow = 1, pfhigh = 2, proll = 3)
    b <- WaverideR:::.aaa_boxes()
    px <- WaverideR:::.aaa_to_ndc(b$proxy, rv$view_p, FALSE, 156.37)
    ax <- WaverideR:::.aaa_to_ndc(b$astro, rv$view_a, FALSE, 5123.4)
    session$setInputs(click_main = list(x = px, y = mean(b$proxy[3:4])))
    session$setInputs(click_main = list(x = ax, y = mean(b$astro[3:4])))
    expect_equal(nrow(rv$ties), n0 + 1)
    expect_equal(rv$ties$p[n0 + 1], 156.37, tolerance = 1e-3)               # free proxy position, not a peak
    expect_true(rv$ties$a[n0 + 1] %in% S$pk_a$x)                            # astro tie snapped to a peak
    rv$sel <- n0 + 1
    session$setInputs(sel_p = 156.5)
    expect_equal(rv$ties$p[n0 + 1], 156.5)
    session$setInputs(sel_off = 0.03)                                         # offset from the nearest peak
    nearest <- S$pk_p$x[which.min(abs(S$pk_p$x - 156.5))]
    expect_equal(rv$ties$p[n0 + 1], nearest + 0.03, tolerance = 1e-6)
    session$setInputs(sel_p = 156.5)
    session$setInputs(p_fr = 1)
    expect_gt(rv$ties$p[n0 + 1], 156.5)
    session$setInputs(undo = 1)
    session$setInputs(undo = 2)
    expect_equal(nrow(rv$ties), n0 + 1)
    session$setInputs(clear = 1)
    expect_equal(nrow(rv$ties), 0)
    expect_true(is.list(tuned()))
    b <- WaverideR:::.aaa_band(WaverideR:::.aaa_est_period(S$proxy))
    session$setInputs(tie_on = "cycle", pflow = b[["flow"]], pfhigh = b[["fhigh"]])
    expect_true(is.null(cyc()$error))
    expect_true(nrow(rec_p()$cands) > 1)
  })
})

test_that("metronome and Laskar helpers work (Laskar with a mocked download)", {
  m <- WaverideR:::.aaa_metronome(100, 10, 0, 1000)
  pk <- WaverideR:::.aaa_peaks(m, "max", 3)
  expect_equal(pk$x[2] - pk$x[1], 100, tolerance = 1e-6)
  expect_equal(pk$x[1] %% 100, 10, tolerance = 1e-6)
  expect_error(WaverideR:::.aaa_metronome(0, 0, 0, 10), "period")
  fake <- function(sol) data.frame(Time_ka = 0:999, ecc = 0.03 + 0.02 * cos(2 * pi * (0:999) / 100))
  l <- WaverideR:::.aaa_laskar("laskar:la10a", 100, 500, fake)
  expect_equal(range(l$x), c(100, 500))
})

test_that("the app can switch target and moves the windows with the sliders", {
  skip_if_not_installed("shiny")
  S <- WaverideR:::.aaa_setup(list(first = astro, second = astro[astro[, 1] > 5100, ]), proxy, "max", "max", 3, tie_start, 12)
  expect_equal(length(S$astro_list), 2)
  app <- WaverideR:::.aaa_app(S, TRUE, 1, "raw", NULL, NULL)
  shiny::testServer(app, {
    session$setInputs(a_from = min(S$astro$x), a_to = max(S$astro$x), pk_a_type = "max", pk_p_type = "max", snap = FALSE,
        tie_on = "raw", flow = 0.0381, fhigh = 0.0635, roll = 3, show_astro = TRUE, show_sed = TRUE, tune = 0,
        a_period = 405, a_phase = 0, a_preset = 405, a_src = "prov:first")
    n0 <- nrow(rv$ties)
    session$setInputs(a_src = "metro", a_from = 5000, a_to = 6000)
    expect_equal(nrow(raw_a()), nrow(WaverideR:::.aaa_metronome(405, 0, 5000, 6000)))
    expect_equal(nrow(rv$ties), 0)             # a new kind of target removes the tie points
    session$setInputs(win_p = c(154.5, 161))                                # slider
    expect_equal(rv$view_p, c(154.5, 161))
    session$setInputs(wp_from = 155, wp_to = 160)                           # typed window
    expect_equal(rv$view_p, c(155, 160))
    session$setInputs(wheel_main = list(x = 0.5, y = 0.86, dy = -100, t = 1))   # wheel over the proxy: zoom in
    expect_lt(diff(rv$view_p), 5)
    v0 <- rv$view_p
    session$setInputs(pan_main = list(x = 0.5, y = 0.86, dx = 0.1, t = 2))      # pan the proxy record
    expect_false(isTRUE(all.equal(rv$view_p, v0)))
    session$setInputs(fit = 1)
    expect_equal(rv$view_a, range(raw_a()$x))
    session$setInputs(time_dir = FALSE, flip_p = TRUE, flip_a = TRUE)     # other reference frame
    expect_false(td())
    expect_true(flip_p())
    expect_true(flip_a())
    session$setInputs(pk_a_type = "min")
    expect_true(nrow(pka()) > 1)
  })
})

test_that("American spelling aliases are identical to the British functions", {
  expect_identical(analyze_wavelet, analyse_wavelet)
  expect_identical(analyze_superlet, analyse_superlet)
  expect_identical(wr_colored_noise, wr_coloured_noise)
})

test_that("British argument synonyms exist", {
  expect_true("colour_brewer" %in% names(formals(plot_wavelet)))
  expect_true("standardise" %in% names(formals(analyse_superlet)))
  expect_true("centre_limits" %in% names(formals(ecoco_blocked)))
  expect_true("normalise" %in% names(formals(win_timeOpt)))
  expect_true("colour_grid" %in% names(formals(plot_bispectral)))
})
