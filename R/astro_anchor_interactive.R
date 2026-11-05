#' @title Interactively tie a proxy record to an astronomical solution (random order, tie and untie, inspect)
#'
#' @description A more flexible alternative to \code{\link{astro_anchor}}. The proxy record (top) and the astronomical
#' solution (bottom) are shown in one window with their peaks (or troughs). Tie points can be set, removed and changed in
#' any order: click a peak of the proxy record and then a peak of the astronomical solution to tie them (or first the
#' astronomical peak and then the proxy peak); click a tied peak again to untie it. The windows can be zoomed and moved
#' with the buttons above the plot, so that long records can be tied piece by piece. Press \emph{Tune and inspect} at any
#' moment: a second window shows the record tuned to the current tie points together with a cycle extracted from it,
#' the astronomical solution and the sedimentation rate between the tie points, so that the result can be judged
#' while tie points are still being added or removed.
#'
#' @details
#' Buttons (click them with the left mouse button): \emph{Tune and inspect} (open or update the inspection window),
#' \emph{Undo} (remove the last change), \emph{Clear} (remove all tie points), \emph{Zoom +}/\emph{Zoom -} (both
#' windows), \emph{P <}, \emph{P >} (move the proxy window), \emph{A <}, \emph{A >} (move the astronomical window),
#' \emph{Fit} (show the whole records), \emph{Sync} (move the astronomical window to the partners of the tie points in the proxy window, or the other way round if there are none), \emph{Wavelet} (switch the wavelet scalogram of the tuned record in the inspection
#' window on or off) and \emph{Done}. A right click also finishes. Tie points that are in the wrong order (the tuned
#' time would run backwards) are drawn in orange, and the tuning refuses to run until they are fixed.
#'
#' The result has the same four columns as the output of \code{\link{astro_anchor}} and can be used directly by
#' \code{\link{plot_astro_anchor}} and \code{\link{anchor2time}}.
#'
#' @param astro_solution Astronomical solution to tie to (matrix or data frame with two columns: time/age (e.g.
#'   kyr) in column 1, value in column 2; rows with \code{NA} are removed and rows are sorted on column 1).
#'   Required, no default (\code{Default=NULL} gives an error).
#' @param proxy_signal Proxy record to be tied (matrix or data frame with two columns: depth/time in column 1,
#'   proxy value in column 2; rows with \code{NA} are removed). Use \code{\link{extract_signal}} first if a single
#'   cycle of the record should be tied. Required, no default (\code{Default=NULL} gives an error).
#' @param proxy_min_or_max Tie the maxima or the minima of the proxy record (character): \code{"max"} or
#'   \code{"min"}. \code{Default="max"}.
#' @param astrosolution_min_or_max Tie the maxima or the minima of the astronomical solution (character):
#'   \code{"max"} or \code{"min"}. \code{Default="max"}.
#' @param pts Number of points to the left and right of a value that the peak detection (\code{\link{max_detect}},
#'   \code{\link{min_detect}}) compares it with (integer, >= 1); larger values give fewer, more certain peaks.
#'   \code{Default=3}.
#' @param ties Existing tie points to start from (data frame with the columns \code{data_x} (proxy depth/time) and
#'   \code{insolation_x} (astronomical age), otherwise the first two columns are used; e.g. the output of
#'   \code{\link{astro_anchor}}). They are snapped to the nearest peaks. \code{NULL} starts without tie points.
#'   \code{Default=NULL}.
#' @param time_dir Direction of time in the proxy record (logical): \code{TRUE} if time increases with increasing
#'   depth/time values, \code{FALSE} if it decreases (the proxy axis is then reversed). \code{Default=TRUE}.
#' @param peaks_in_view Approximate number of peaks shown in each window at the start (numeric, > 0).
#'   \code{Default=12}.
#' @param inspect_period Period of the cycle extracted from the tuned record in the inspection window, in the time
#'   unit of the astronomical solution (numeric, > 0); \code{NULL} uses the median spacing of the peaks of the
#'   astronomical solution. \code{Default=NULL}.
#' @param inspect_band Relative half-width of the extracted band (numeric, between 0 and 1): \code{0.25} extracts
#'   0.75 to 1.25 times \code{inspect_period}. \code{Default=0.25}.
#' @param sedrate_factor Factor applied to the sedimentation rate (depth units per time unit) shown in the inspection
#'   window (numeric, > 0), for example 100 to go from m/kyr to cm/kyr. \code{Default=1}.
#' @param inspect_wavelet Also show the wavelet scalogram of the tuned record in the inspection window at the start
#'   (logical, TRUE or FALSE); it can be switched on and off with the \emph{Wavelet} button. \code{Default=FALSE}.
#' @param locator_fun Function without arguments that returns the next mouse click in normalised device
#'   coordinates as a list with \code{x} and \code{y}, like \code{\link[graphics]{locator}(1)}; returning
#'   \code{NULL} finishes. \code{NULL} uses \code{graphics::locator(1)}. Another function (for example one that
#'   replays a list of clicks) allows scripted use and testing. When no interactive device is available and
#'   \code{ties} is given, the snapped tie points are returned without opening a window. \code{Default=NULL}.
#' @return A data frame with the tie points sorted by depth/time: \code{data_x}, \code{data_y} (proxy),
#'   \code{insolation_x}, \code{insolation_y} (astronomical solution).
#' @examples
#' \donttest{
#' # Tie points can be given as a starting point and are snapped to the nearest peaks;
#' # without an interactive graphics device (as here) they are returned as they are.
#' astrosignal <- na.omit(astrosignal_example)
#' tie_start <- data.frame(data_x = c(155, 160, 165), insolation_x = c(5040, 5400, 5560))
#' proxy <- grey[grey[, 1] < 175, ]
#' ties <- astro_anchor_interactive(astro_solution = astrosignal,
#'                                  proxy_signal = proxy,
#'                                  ties = tie_start)
#' ties
#' # in an interactive session, leave out 'ties' to tie the peaks with the mouse
#' }
#' @seealso \code{\link{astro_anchor}}, \code{\link{plot_astro_anchor}}, \code{\link{anchor2time}}
#' @export
#' @importFrom graphics locator rect text segments lines points plot.new plot.window axis box title abline polygon legend
#' @importFrom graphics grconvertX grconvertY
#' @importFrom grDevices dev.cur dev.list dev.set dev.off
#' @importFrom stats approx median cor sd
astro_anchor_interactive <- function(astro_solution = NULL,
    proxy_signal = NULL,
    proxy_min_or_max = "max",
    astrosolution_min_or_max = "max",
    pts = 3,
    ties = NULL,
    time_dir = TRUE,
    peaks_in_view = 12,
    inspect_period = NULL,
    inspect_band = 0.25,
    sedrate_factor = 1,
    inspect_wavelet = FALSE,
    locator_fun = NULL) {
  astro <- as.data.frame(na.omit(as.data.frame(astro_solution)))[, 1:2]
  proxy <- as.data.frame(na.omit(as.data.frame(proxy_signal)))[, 1:2]
  astro <- astro[order(astro[, 1]), ]
  proxy <- proxy[order(proxy[, 1]), ]
  names(astro) <- c("x", "y")
  names(proxy) <- c("x", "y")

  pk_a <- if (astrosolution_min_or_max == "max") { max_detect(astro, pts = pts) } else { min_detect(astro, pts = pts) }
  pk_p <- if (proxy_min_or_max == "max") { max_detect(proxy, pts = pts) } else { min_detect(proxy, pts = pts) }
  pk_a <- as.data.frame(pk_a[, 1:2])
  pk_p <- as.data.frame(pk_p[, 1:2])
  names(pk_a) <- names(pk_p) <- c("x", "y")
  if (nrow(pk_a) < 2 || nrow(pk_p) < 2) {
    stop("fewer than two peaks found; change 'pts' or the min/max setting")
  }

  st <- list(ties = data.frame(p = integer(0), a = integer(0)), pending = NULL, history = list(),
      view_p = .aai_first_view(pk_p$x, proxy$x, peaks_in_view),
      view_a = .aai_first_view(pk_a$x, astro$x, peaks_in_view),
      wavelet = isTRUE(inspect_wavelet))

  if (!is.null(ties)) {
    ties <- as.data.frame(ties)
    tx <- if ("data_x" %in% names(ties)) { ties$data_x } else { ties[, 1] }
    ta <- if ("insolation_x" %in% names(ties)) { ties$insolation_x } else { ties[, 2] }
    ok <- !is.na(tx) & !is.na(ta)
    tx_ok <- tx[ok]
    ta_ok <- ta[ok]
    p_i <- integer(length(tx_ok))
    for (i in seq_along(tx_ok)) {
      p_i[i] <- which.min(abs(pk_p$x - tx_ok[i]))
    }
    a_i <- integer(length(ta_ok))
    for (i in seq_along(ta_ok)) {
      a_i[i] <- which.min(abs(pk_a$x - ta_ok[i]))
    }
    keep <- !duplicated(p_i) & !duplicated(a_i)
    st$ties <- data.frame(p = p_i[keep], a = a_i[keep])
    if (nrow(st$ties) > 0) {      # start with the windows around the first tie point
      f <- which.min(pk_p$x[st$ties$p])
      st$view_p <- .aai_fix_view(pk_p$x[st$ties$p[f]] + c(-0.5, 0.5) * diff(st$view_p), range(proxy$x), diff(st$view_p))
      st$view_a <- .aai_fix_view(pk_a$x[st$ties$a[f]] + c(-0.5, 0.5) * diff(st$view_a), range(astro$x), diff(st$view_a))
    }
  }

  have_locator <- !is.null(locator_fun)
  if (!have_locator) {
    if (!interactive()) {
      if (is.null(ties)) {
        stop("astro_anchor_interactive() needs mouse input and can only be used in an interactive R session ",
            "(or with 'ties' to snap existing tie points).", call. = FALSE)
      }
      return(.aai_output(st, pk_p, pk_a))
    }
    wr_open_device(width = 14, height = 8, interactive_input = TRUE)
    locator_fun <- function() graphics::locator(1)
  }
  main_dev <- dev.cur()
  if (main_dev == 1) {
    grDevices::pdf(NULL)
    main_dev <- dev.cur()
  }
  op_main <- graphics::par(no.readonly = TRUE)
  on.exit(try({
    if (main_dev %in% dev.list()) {
      dev.set(main_dev)
      graphics::par(op_main)
    }
  }, silent = TRUE), add = TRUE)
  inspect_dev <- NULL

  repeat {
    dev.set(main_dev)
    geo <- .aai_draw(st, proxy, astro, pk_p, pk_a, time_dir)
    .aai_env$geo <- geo
    cl <- locator_fun()
    if (is.null(cl) || length(cl$x) == 0) {
      break
    }
    act <- .aai_hit(geo, cl$x[1], cl$y[1])
    if (is.null(act)) {
      next
    }
    if (act$type == "button") {
      b <- act$name
      if (b == "Done") {
        break
      }
      if (b == "Undo") {
        if (length(st$history) > 0) {
          last <- st$history[[length(st$history)]]
          st$ties <- last$ties
          st$pending <- last$pending
          st$history <- st$history[-length(st$history)]
        }
      } else if (b == "Clear") {
        st <- .aai_push(st)
        st$ties <- st$ties[0, ]
        st$pending <- NULL
      } else if (b == "Zoom +") {
        st$view_p <- .aai_zoom(st$view_p, 0.5, range(proxy$x))
        st$view_a <- .aai_zoom(st$view_a, 0.5, range(astro$x))
      } else if (b == "Zoom -") {
        st$view_p <- .aai_zoom(st$view_p, 2, range(proxy$x))
        st$view_a <- .aai_zoom(st$view_a, 2, range(astro$x))
      } else if (b == "Fit") {
        st$view_p <- range(proxy$x)
        st$view_a <- range(astro$x)
      } else if (b == "Sync") {
        st <- .aai_sync(st, pk_p, pk_a, proxy, astro)
      } else if (b %in% c("P <", "P >")) {
        st$view_p <- .aai_pan(st$view_p, if (b == "P <") { -0.5 } else { 0.5 }, range(proxy$x))
      } else if (b %in% c("A <", "A >")) {
        st$view_a <- .aai_pan(st$view_a, if (b == "A <") { -0.5 } else { 0.5 }, range(astro$x))
      } else if (b == "Wavelet") {
        st$wavelet <- !st$wavelet
      } else if (b == "Tune and inspect") {
        out <- .aai_output(st, pk_p, pk_a)
        inspect_dev <- .aai_inspect_device(inspect_dev)
        res <- tryCatch(
            .aai_inspect(proxy, astro, pk_a, out, time_dir, inspect_period, inspect_band, sedrate_factor, st$wavelet),
            error = function(e) {
              message("Tuning is not possible: ", conditionMessage(e))
              NULL
            })
        dev.set(main_dev)
      }
    } else {
      st <- .aai_click(st, act$panel, act$x, pk_p, pk_a, push = TRUE)
    }
  }
  if (!is.null(inspect_dev) && inspect_dev %in% dev.list()) {
    # leave the inspection window open for the user in an interactive session only
    if (!interactive() || have_locator) {
      try(dev.off(inspect_dev), silent = TRUE)
    }
  }
  .aai_output(st, pk_p, pk_a)
}

# Internal state and geometry

.aai_env <- new.env(parent = emptyenv())

.aai_first_view <- function(peak_x, all_x, n) {
  r <- range(all_x)
  if (length(peak_x) < 2) {
    return(r)
  }
  w <- stats::median(diff(peak_x)) * n
  c(r[1], min(r[2], r[1] + w))
}

.aai_zoom <- function(v, f, lim) {
  mid <- mean(v)
  half <- diff(v) / 2 * f
  half <- min(half, diff(lim) / 2)
  .aai_fix_view(c(mid - half, mid + half), lim, 2 * half)
}

.aai_pan <- function(v, frac, lim) {
  w <- diff(v)
  new <- v + frac * w
  .aai_fix_view(new, lim, w)
}

# keep a window of width w inside lim
.aai_fix_view <- function(v, lim, w) {
  w <- min(w, diff(lim))
  if (v[1] < lim[1]) {
    v <- c(lim[1], lim[1] + w)
  }
  if (v[2] > lim[2]) {
    v <- c(lim[2] - w, lim[2])
  }
  v
}

# move the window of one record to the partners of the tie points visible in the other (proxy first)
.aai_sync <- function(st, pk_p, pk_a, proxy, astro) {
  if (nrow(st$ties) == 0) {
    return(st)
  }
  xp <- pk_p$x[st$ties$p]
  xa <- pk_a$x[st$ties$a]
  in_p <- xp >= st$view_p[1] & xp <= st$view_p[2]
  in_a <- xa >= st$view_a[1] & xa <= st$view_a[2]
  if (any(in_p)) {
    w <- diff(st$view_a)
    st$view_a <- .aai_fix_view(mean(range(xa[in_p])) + c(-0.5, 0.5) * w, range(astro$x), w)
  } else if (any(in_a)) {
    w <- diff(st$view_p)
    st$view_p <- .aai_fix_view(mean(range(xp[in_a])) + c(-0.5, 0.5) * w, range(proxy$x), w)
  }
  st
}

.aai_push <- function(st) {
  st$history[[length(st$history) + 1]] <- list(ties = st$ties, pending = st$pending)
  st
}

# click on a peak of one record: select, tie or untie
.aai_click <- function(st, panel, x, pk_p, pk_a, push = TRUE) {
  pk <- if (panel == "proxy") { pk_p } else { pk_a }
  idx <- which.min(abs(pk$x - x))
  col <- if (panel == "proxy") { "p" } else { "a" }
  other <- if (panel == "proxy") { "a" } else { "p" }
  if (push) {
    st <- .aai_push(st)
  }
  hit <- which(st$ties[[col]] == idx)
  if (length(hit) > 0) {                       # tied peak: untie
    st$ties <- st$ties[-hit, , drop = FALSE]
    st$pending <- NULL
    return(st)
  }
  if (!is.null(st$pending)) {
    if (st$pending$panel == panel) {           # same record again: change or drop the selection
      st$pending <- if (st$pending$idx == idx) { NULL } else { list(panel = panel, idx = idx) }
    } else {                                   # other record: tie
      p_i <- if (panel == "proxy") { idx } else { st$pending$idx }
      a_i <- if (panel == "astro") { idx } else { st$pending$idx }
      st$ties <- rbind(st$ties, data.frame(p = p_i, a = a_i))
      st$pending <- NULL
    }
  } else {
    st$pending <- list(panel = panel, idx = idx)
  }
  st
}

# tie points are in the wrong order when the tuned time would run backwards
.aai_flag <- function(ties, pk_p, pk_a, time_dir) {
  n <- nrow(ties)
  flag <- rep(FALSE, n)
  if (n < 2) {
    return(flag)
  }
  o <- order(pk_p$x[ties$p])
  tp <- pk_p$x[ties$p][o]
  ta <- pk_a$x[ties$a][o]
  d <- diff(ta)
  bad <- if (time_dir) { d <= 0 } else { d >= 0 }
  fo <- rep(FALSE, n)
  fo[-1] <- fo[-1] | bad
  fo[-n] <- fo[-n] | bad
  flag[o] <- fo
  flag
}

.aai_output <- function(st, pk_p, pk_a) {
  if (nrow(st$ties) == 0) {
    return(data.frame(data_x = numeric(0), data_y = numeric(0), insolation_x = numeric(0), insolation_y = numeric(0)))
  }
  out <- data.frame(data_x = pk_p$x[st$ties$p], data_y = pk_p$y[st$ties$p],
      insolation_x = pk_a$x[st$ties$a], insolation_y = pk_a$y[st$ties$a])
  out <- out[order(out$data_x), ]
  rownames(out) <- NULL
  out
}

# boxes (normalised device coordinates) of the two panels and the buttons
.aai_boxes <- function() {
  list(proxy = c(0.07, 0.97, 0.52, 0.89), astro = c(0.07, 0.97, 0.08, 0.37))
}

.aai_buttons <- function() {
  labs <- c("Tune and inspect", "Undo", "Clear", "Zoom +", "Zoom -", "P <", "P >", "A <", "A >", "Fit", "Sync", "Wavelet", "Done")
  w <- c(0.13, 0.055, 0.055, 0.06, 0.06, 0.045, 0.045, 0.045, 0.045, 0.045, 0.055, 0.07, 0.055)
  gap <- 0.008
  x0 <- 0.07
  b <- list()
  for (i in seq_along(labs)) {
    b[[labs[i]]] <- c(x0, x0 + w[i], 0.935, 0.985)
    x0 <- x0 + w[i] + gap
  }
  b
}

.aai_draw <- function(st, proxy, astro, pk_p, pk_a, time_dir) {
  # no par() restore here: locator() must see the coordinate system of the overlay drawn last
  bx <- .aai_boxes()
  flag <- .aai_flag(st$ties, pk_p, pk_a, time_dir)
  panel <- function(name, dat, pk, view, tied, pend, reversed, ylab, flagged, new) {
    b <- bx[[name]]
    graphics::par(fig = c(0, 1, 0, 1), plt = b[c(1, 2, 3, 4)], new = new, mgp = c(2.2, 0.7, 0), cex.axis = 0.8)
    inside <- dat[, 1] >= view[1] & dat[, 1] <= view[2]
    yr <- range(dat[inside, 2], na.rm = TRUE)
    if (!is.finite(diff(yr)) || diff(yr) == 0) {
      yr <- range(dat[, 2])
    }
    yr <- yr + c(-0.05, 0.05) * diff(yr)
    xl <- if (reversed) { rev(view) } else { view }
    graphics::plot(dat[, 1], dat[, 2], type = "l", xlim = xl, ylim = yr, xaxs = "i", yaxs = "i", xlab = "", ylab = ylab)
    graphics::box()
    graphics::points(pk$x, pk$y, pch = 1, col = "grey40")
    if (length(tied) > 0) {
      graphics::points(pk$x[tied], pk$y[tied], pch = 19, cex = 1.4, col = ifelse(flagged, "orange", "green3"))
    }
    if (!is.null(pend)) {
      graphics::points(pk$x[pend], pk$y[pend], pch = 19, cex = 1.6, col = "red")
    }
    c(b, view, yr, if (reversed) { 1 } else { 0 })
  }
  g_p <- panel("proxy", proxy, pk_p, st$view_p, st$ties$p, if (!is.null(st$pending) && st$pending$panel == "proxy") { st$pending$idx } else { NULL },
      !time_dir, "proxy", flag, FALSE)
  g_a <- panel("astro", astro, pk_a, st$view_a, st$ties$a, if (!is.null(st$pending) && st$pending$panel == "astro") { st$pending$idx } else { NULL },
      FALSE, "astronomical solution", flag, TRUE)
  # overlay in normalised device coordinates: tie lines, buttons, status
  graphics::par(fig = c(0, 1, 0, 1), plt = c(0, 1, 0, 1), new = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  to_ndc <- function(g, x) {
    f <- (x - g[5]) / (g[6] - g[5])
    if (g[9] == 1) {
      f <- 1 - f
    }
    g[1] + f * (g[2] - g[1])
  }
  y_ndc <- function(g, y) g[3] + (y - g[7]) / (g[8] - g[7]) * (g[4] - g[3])
  if (nrow(st$ties) > 0) {
    stub <- 0.045
    for (k in seq_len(nrow(st$ties))) {
      xp <- pk_p$x[st$ties$p[k]]
      xa <- pk_a$x[st$ties$a[k]]
      col <- if (flag[k]) { "orange" } else { "green4" }
      p_in <- xp >= st$view_p[1] && xp <= st$view_p[2]
      a_in <- xa >= st$view_a[1] && xa <= st$view_a[2]
      xpn <- if (p_in) { to_ndc(g_p, xp) } else { NA }
      xan <- if (a_in) { to_ndc(g_a, xa) } else { NA }
      if (p_in) {
        graphics::segments(xpn, y_ndc(g_p, pk_p$y[st$ties$p[k]]), xpn, bx$proxy[3], lwd = 1, lty = 3, col = col)
      }
      if (a_in) {
        graphics::segments(xan, y_ndc(g_a, pk_a$y[st$ties$a[k]]), xan, bx$astro[4], lwd = 1, lty = 3, col = col)
      }
      if (p_in && a_in) {
        graphics::segments(xpn, bx$proxy[3], xan, bx$astro[4], lwd = 2, col = col)
      } else if (p_in) {           # partner outside the window of the other record: short stub with its position
        graphics::segments(xpn, bx$proxy[3], xpn, bx$proxy[3] - stub, lwd = 2, col = col)
        graphics::text(xpn, bx$proxy[3] - stub - 0.012, paste0("> ", signif(xa, 5)), cex = 0.6, col = col)
      } else if (a_in) {
        graphics::segments(xan, bx$astro[4], xan, bx$astro[4] + stub, lwd = 2, col = col)
        graphics::text(xan, bx$astro[4] + stub + 0.012, paste0("> ", signif(xp, 5)), cex = 0.6, col = col)
      }
    }
  }
  btn <- .aai_buttons()
  for (nm in names(btn)) {
    r <- btn[[nm]]
    on <- nm == "Wavelet" && isTRUE(st$wavelet)
    graphics::rect(r[1], r[3], r[2], r[4], col = if (on) { "lightgreen" } else { "grey90" }, border = "grey30")
    graphics::text(mean(r[1:2]), mean(r[3:4]), nm, cex = 0.75)
  }
  msg <- sprintf("tie points: %d%s%s   (click a peak of each record to tie, click a tied peak to untie)", nrow(st$ties),
      if (any(flag)) { sprintf(", %d in the wrong order (orange)", sum(flag)) } else { "" },
      if (!is.null(st$pending)) { sprintf("   selected: %s peak (red)", st$pending$panel) } else { "" })
  graphics::text(0.07, 0.915, msg, adj = c(0, 0.5), cex = 0.75)
  graphics::text(0.07, 0.395, "astronomical solution", adj = c(0, 0.5), cex = 0.7, col = "grey30")
  list(proxy = g_p, astro = g_a, buttons = btn)
}

# which button or panel was clicked (x, y in normalised device coordinates)
.aai_hit <- function(geo, x, y) {
  for (nm in names(geo$buttons)) {
    r <- geo$buttons[[nm]]
    if (x >= r[1] && x <= r[2] && y >= r[3] && y <= r[4]) {
      return(list(type = "button", name = nm))
    }
  }
  for (nm in c("proxy", "astro")) {
    g <- geo[[nm]]
    pad <- 0.03
    if (x >= g[1] - pad && x <= g[2] + pad && y >= g[3] - pad && y <= g[4] + pad) {
      f <- (min(max(x, g[1]), g[2]) - g[1]) / (g[2] - g[1])
      if (g[9] == 1) {
        f <- 1 - f
      }
      return(list(type = "panel", panel = nm, x = g[5] + f * (g[6] - g[5])))
    }
  }
  NULL
}

# Tuning and inspection

.aai_inspect_device <- function(dev) {
  if (!is.null(dev) && dev %in% dev.list()) {
    dev.set(dev)
    return(dev)
  }
  if (interactive()) {
    try(wr_open_device(width = 12, height = 9), silent = TRUE)
  } else {
    grDevices::pdf(NULL)
  }
  dev.cur()
}

# tune the record to the tie points and draw the inspection plots on the current device
.aai_inspect <- function(proxy, astro, pk_a, ties, time_dir, period, band, sedrate_factor, show_wavelet) {
  if (nrow(ties) < 2) {
    stop("at least two tie points are needed")
  }
  t_ord <- ties[order(ties$data_x), ]
  bad <- if (time_dir) { diff(t_ord$insolation_x) <= 0 } else { diff(t_ord$insolation_x) >= 0 }
  if (any(bad)) {
    stop("tie points are in the wrong order (shown in orange); untie or change them first")
  }
  tuned <- anchor2time(anchor_points = cbind(t_ord$data_x, t_ord$insolation_x), data = proxy)
  tuned <- as.data.frame(tuned)
  names(tuned) <- c("time", "value")
  tuned <- tuned[order(tuned$time), ]
  tuned <- tuned[is.finite(tuned$time) & is.finite(tuned$value), ]
  dt <- stats::median(diff(tuned$time))
  if (!is.finite(dt) || dt <= 0) {
    stop("the tuned time axis is not increasing")
  }
  n_grid <- min(20000, ceiling(diff(range(tuned$time)) / dt) + 1)
  grid <- seq(min(tuned$time), max(tuned$time), length.out = n_grid)
  reg <- data.frame(time = grid, value = stats::approx(tuned$time, tuned$value, grid)$y)
  dtg <- grid[2] - grid[1]
  if (is.null(period)) {
    period <- stats::median(diff(pk_a$x))
  }
  lo <- max(2 * dtg, period * 0.35)
  up <- min(diff(range(grid)) / 2, period * 3)
  if (up <= lo * 1.5) {
    stop("the tuned record is too short for a cycle of period ", signif(period, 3))
  }
  w <- analyse_wavelet(reg, dj = 1 / 40, lowerPeriod = lo, upperPeriod = up, verbose = FALSE, omega_nr = 6)
  ext <- extract_signal_stable(w, cycle = period, period_up = 1 + band, period_down = 1 - band, add_mean = FALSE)
  z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(matrix(1:3, ncol = 1), heights = c(1.1, 1.1, 0.8))
  graphics::par(mar = c(1.2, 4.2, 2, 1), mgp = c(2.4, 0.7, 0))
  xr <- range(grid)
  graphics::plot(reg$time, z(reg$value), type = "l", col = "grey35", xlim = xr, xaxt = "n", ylab = "z-score",
      main = sprintf("Tuned record (%d tie points) and the extracted %s cycle", nrow(t_ord), signif(period, 3)))
  graphics::lines(ext[, 1], z(ext[, 2]) * 1, col = "red", lwd = 2)
  graphics::abline(v = t_ord$insolation_x, col = "green4", lty = 3)
  graphics::legend("topright", legend = c("tuned proxy", "extracted cycle", "tie points"), col = c("grey35", "red", "green4"),
      lty = c(1, 1, 3), lwd = c(1, 2, 1), bty = "n", cex = 0.8)
  ai <- astro[astro$x >= xr[1] & astro$x <= xr[2], ]
  a_ext <- stats::approx(astro$x, astro$y, ext[, 1])$y
  r <- suppressWarnings(stats::cor(ext[, 2], a_ext, use = "complete.obs"))
  graphics::par(mar = c(1.2, 4.2, 1, 1))
  graphics::plot(ext[, 1], z(ext[, 2]), type = "l", col = "red", lwd = 2, xlim = xr, xaxt = "n", ylab = "z-score",
      ylim = range(c(z(ext[, 2]), z(ai$y)), na.rm = TRUE),
      main = sprintf("Extracted cycle (red) against the astronomical solution (blue), r = %s", format(round(r, 2), nsmall = 2)))
  graphics::lines(ai$x, z(ai$y), col = "blue", lwd = 1.5)
  graphics::abline(v = t_ord$insolation_x, col = "green4", lty = 3)
  graphics::par(mar = c(4, 4.2, 1, 1))
  sr <- diff(t_ord$data_x) / diff(t_ord$insolation_x) * sedrate_factor
  mid <- (t_ord$insolation_x[-1] + t_ord$insolation_x[-nrow(t_ord)]) / 2
  graphics::plot(range(t_ord$insolation_x), range(c(0, sr), na.rm = TRUE), type = "n", xlim = xr, xlab = "time",
      ylab = "sedimentation rate")
  graphics::segments(t_ord$insolation_x[-nrow(t_ord)], sr, t_ord$insolation_x[-1], sr, lwd = 3, col = "darkorange3")
  graphics::abline(v = t_ord$insolation_x, col = "green4", lty = 3)
  if (isTRUE(show_wavelet)) {
    if (interactive()) {
      try(wr_open_device(width = 10, height = 6), silent = TRUE)
    } else {
      grDevices::pdf(NULL)
    }
    try(plot_wavelet(w, add_data = TRUE, add_avg = TRUE, dev_new = FALSE, lowerPeriod = lo, upperPeriod = up,
        periodlab = "Period (time units)", x_lab = "time"), silent = TRUE)
  }
  invisible(list(tuned = tuned, extracted = ext, correlation = r, sedrate = sr))
}
