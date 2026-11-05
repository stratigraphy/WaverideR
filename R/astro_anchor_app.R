#' @title Interactive Shiny app to tie a proxy record to an astronomical solution (tune live while tying)
#'
#' @description A browser-based version of \code{\link{astro_anchor_interactive}}. One window shows three linked rows:
#' the proxy record, the astronomical solution and, underneath, the record tuned to the current tie points together
#' with a cycle extracted from it (and optionally the astronomical solution and the sedimentation rate). The tuned row
#' is updated while tie points are made, moved or removed, so the effect of every change on the fit can be judged on
#' the go. Tie points are numbered and colour coded in all rows and joined by lines between the two records.
#'
#' @details
#' \strong{What is tied.} Tie points on the astronomical solution always sit on a maximum or minimum (see
#' \code{astrosolution_min_or_max}). On the proxy side they can be placed on the \emph{raw} record or on a cycle that is
#' extracted from it with a Taner filter (\code{\link[astrochron]{taner}}; choose \emph{Cycle extracted from the proxy}
#' and give \code{flow}, \code{fhigh} and the exponent of the roll-off, \code{roll = 10^x}). The tuning is always applied
#' to the raw record. Proxy tie points do not have to lie on a maximum or minimum: with \emph{Snap to peaks} switched off
#' they are placed exactly where you click, and every proxy tie point can be moved left and right with the buttons (fine
#' and coarse steps), by typing its position, by typing its offset from the nearest maximum/minimum, or by dragging.
#' The cycle shown in the tuned record is extracted with the same kind of Taner filter.
#'
#' \strong{Tie} a point: click one record (red) and then the other (in any order). \strong{Select} a tie point by
#' clicking it. \strong{Move} it by clicking a free position of one of the records afterwards, by dragging over the record
#' from the tie point to the new position, by the arrow buttons, or by editing the numbers. \strong{Untie} it by double
#' clicking it or with the \emph{Untie selected} button. Tie points in the wrong order (the tuned time would run
#' backwards) are marked with a black cross; the tuned row waits until they are fixed.
#'
#' \strong{Target and windows.} The card \emph{Astronomical solution} selects the target (a provided solution, a
#' metronome, or a Laskar solution), the time interval it spans and whether maxima or minima are tied. Two sliders above the
#' plot show the visible window (from, to) of the proxy record and of the astronomical solution; type new numbers to set
#' them. With the mouse, the \emph{wheel} zooms the record under the pointer around the pointer and \emph{Shift + drag}
#' (or dragging with the right or middle mouse button) moves it from left to right. \emph{Sync} brings the astronomical window to the partners of the tie points in the proxy window, \emph{Show both records in full} resets both.
#'
#' The sidebar also has buttons to undo changes and to finish (\emph{Done}). Closing the browser tab also finishes and
#' returns the tie points.
#'
#' The app runs on the local computer only (a Shiny app started with \code{\link[shiny]{runApp}}); nothing is sent to
#' other servers. The package \pkg{shiny} has to be installed.
#'
#' The result has the same four columns as the output of \code{\link{astro_anchor}} and can be used by
#' \code{\link{plot_astro_anchor}} and \code{\link{anchor2time}}. \code{data_y} and \code{insolation_y} are the values of
#' the records at the tie positions.
#'
#' @param astro_solution Astronomical solution to tie to (matrix or data frame with two columns: time/age (e.g. kyr)
#'   in column 1, value in column 2), or a named list of such solutions to choose from in the app (the first one is
#'   shown at the start). In the app the target can also be switched to a metronome (a cosine with a period, phase
#'   and time interval you set), or to a Laskar solution (eccentricity, obliquity, climatic precession or ETP of
#'   La2004, La2010a-d, La2011 eccentricity, or 65 N insolation) that is downloaded on request with
#'   \code{\link[astrochron]{getLaskar}} when you press the download button. Required, no default
#'   (\code{Default=NULL} gives an error).
#' @param proxy_signal Proxy record to be tied (matrix or data frame with two columns: depth/time in column 1, proxy
#'   value in column 2; rows with \code{NA} are removed). Required, no default (\code{Default=NULL} gives an error).
#' @param proxy_min_or_max Peaks offered when snapping on the proxy record (character): \code{"max"} (maxima) or
#'   \code{"min"} (minima). \code{Default="max"}.
#' @param astrosolution_min_or_max Peaks on which the tie points of the astronomical solution are placed
#'   (character): \code{"max"} (maxima) or \code{"min"} (minima). \code{Default="max"}.
#' @param pts Number of points to the left and right of a value that the peak detection (\code{\link{max_detect}},
#'   \code{\link{min_detect}}) compares it with (integer, >= 1); larger values give fewer, more certain peaks.
#'   \code{Default=3}.
#' @param ties Existing tie points to start from (data frame with the columns \code{data_x} (proxy depth/time) and
#'   \code{insolation_x} (astronomical age), otherwise the first two columns are used; e.g. the output of
#'   \code{\link{astro_anchor}}). The astronomical positions are moved to the nearest maximum/minimum, the proxy
#'   positions are used as given. \code{NULL} starts without tie points. \code{Default=NULL}.
#' @param time_dir Direction of time in the proxy record (logical): \code{TRUE} if time increases with increasing
#'   depth/time values, \code{FALSE} if it decreases. \code{Default=TRUE}.
#' @param peaks_in_view Approximate number of peaks shown in each window at the start (numeric, > 0).
#'   \code{Default=12}.
#' @param tie_on Proxy record on which the tie points are placed at the start (character): \code{"raw"} (the raw
#'   record) or \code{"cycle"} (a cycle extracted from it with the Taner filter set by \code{proxy_flow},
#'   \code{proxy_fhigh} and \code{proxy_roll}). Can be changed in the app. \code{Default=c("raw", "cycle")}, i.e.
#'   \code{"raw"}.
#' @param proxy_flow,proxy_fhigh Low and high cut-off frequencies of the Taner filter that extracts the proxy cycle,
#'   in cycles per depth/time unit of the proxy record (numeric, 0 < \code{proxy_flow} < \code{proxy_fhigh}; give both).
#'   \code{NULL} uses the strongest period of the wavelet spectrum of the proxy -/+ 25 percent. They can be changed in
#'   the app. \code{Default=NULL}.
#' @param proxy_roll Roll-off of the proxy Taner filter as the exponent \code{x} of \code{roll = 10^x} (numeric; see
#'   \code{\link[astrochron]{taner}}). \code{Default=3}.
#' @param inspect_flow,inspect_fhigh Low and high cut-off frequencies of the Taner filter that extracts the cycle from
#'   the tuned record, in cycles per time unit of the astronomical solution (numeric, 0 < \code{inspect_flow} <
#'   \code{inspect_fhigh}; give both). \code{NULL} uses the median spacing of the peaks of the astronomical solution
#'   -/+ 25 percent. \code{Default=NULL}.
#' @param inspect_roll Roll-off of the Taner filter of the tuned record as the exponent \code{x} of
#'   \code{roll = 10^x} (numeric). \code{Default=3}.
#' @param sedrate_factor Factor applied to the sedimentation rate (depth units per time unit) that is shown (numeric,
#'   > 0), for example 100 to go from m/kyr to cm/kyr. \code{Default=1}.
#' @param launch_browser Open the app in the default web browser (logical, TRUE or FALSE); \code{FALSE} prints the
#'   local address (or uses the RStudio viewer when R is run there). Passed to \code{shiny::runApp} as
#'   \code{launch.browser}. \code{Default=TRUE}.
#' @return A data frame with the tie points sorted by depth/time: \code{data_x}, \code{data_y} (proxy),
#'   \code{insolation_x}, \code{insolation_y} (astronomical solution). Without an interactive session and with
#'   \code{ties} the starting tie points are returned without starting the app (those on the astronomical solution
#'   moved to the nearest peak).
#' @examples
#' \donttest{
#' # Without an interactive session (as here) the starting tie points are returned as they are.
#' astrosignal <- na.omit(astrosignal_example)
#' tie_start <- data.frame(data_x = c(155, 160, 165), insolation_x = c(5040, 5400, 5560))
#' proxy <- grey[grey[, 1] < 175, ]
#' astro_anchor_app(astro_solution = astrosignal, proxy_signal = proxy, ties = tie_start)
#' # in an interactive session, leave out 'ties' to start the app and tie the records with the mouse
#' }
#' @seealso \code{\link{astro_anchor_interactive}}, \code{\link{astro_anchor}}, \code{\link{plot_astro_anchor}},
#'   \code{\link{anchor2time}}
#' @export
astro_anchor_app <- function(astro_solution = NULL,
    proxy_signal = NULL,
    proxy_min_or_max = "max",
    astrosolution_min_or_max = "max",
    pts = 3,
    ties = NULL,
    time_dir = TRUE,
    peaks_in_view = 12,
    tie_on = c("raw", "cycle"),
    proxy_flow = NULL,
    proxy_fhigh = NULL,
    proxy_roll = 3,
    inspect_flow = NULL,
    inspect_fhigh = NULL,
    inspect_roll = 3,
    sedrate_factor = 1,
    launch_browser = TRUE) {
  tie_on <- match.arg(tie_on)
  S <- .aaa_setup(astro_solution, proxy_signal, proxy_min_or_max, astrosolution_min_or_max, pts, ties, peaks_in_view)
  if (!interactive()) {
    if (is.null(ties)) {
      stop("astro_anchor_app() needs an interactive R session (or 'ties' to return existing tie points).", call. = FALSE)
    }
    return(.aaa_output(S$st$ties, S$proxy, S$astro))
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("package 'shiny' is needed for astro_anchor_app(); install it with install.packages(\"shiny\") ",
        "or use astro_anchor_interactive().", call. = FALSE)
  }
  app <- .aaa_app(S, time_dir, sedrate_factor, tie_on, c(proxy_flow, proxy_fhigh, proxy_roll), c(inspect_flow, inspect_fhigh, inspect_roll))
  shiny::runApp(app, launch.browser = launch_browser, quiet = TRUE)
}

# Set up

.aaa_peaks <- function(dat, min_or_max, pts) {
  pk <- if (min_or_max == "max") { max_detect(dat, pts = pts) } else { min_detect(dat, pts = pts) }
  pk <- as.data.frame(pk[, 1:2])
  names(pk) <- c("x", "y")
  pk
}

.aaa_setup <- function(astro_solution, proxy_signal, proxy_min_or_max, astrosolution_min_or_max, pts, ties, peaks_in_view) {
  alist <- if (is.data.frame(astro_solution) || is.matrix(astro_solution)) { list(`Provided solution` = astro_solution) } else { astro_solution }
  if (is.null(names(alist))) {
    names(alist) <- paste("Provided solution", seq_along(alist))
  }
  alist_in <- alist
  alist <- vector("list", length(alist_in))
  names(alist) <- names(alist_in)
  for (i in seq_along(alist_in)) {
    a <- as.data.frame(na.omit(as.data.frame(alist_in[[i]])))[, 1:2]
    a <- a[order(a[, 1]), ]
    names(a) <- c("x", "y")
    alist[[i]] <- a
  }
  astro <- alist[[1]]
  proxy <- as.data.frame(na.omit(as.data.frame(proxy_signal)))[, 1:2]
  proxy <- proxy[order(proxy[, 1]), ]
  names(proxy) <- c("x", "y")
  pk_a <- .aaa_peaks(astro, astrosolution_min_or_max, pts)
  pk_p <- .aaa_peaks(proxy, proxy_min_or_max, pts)
  if (nrow(pk_a) < 2 || nrow(pk_p) < 2) {
    stop("fewer than two peaks found; change 'pts' or the min/max setting")
  }
  st <- list(ties = data.frame(p = numeric(0), a = numeric(0)), pending = NULL, sel = NULL,
      view_p = .aai_first_view(pk_p$x, proxy$x, peaks_in_view),
      view_a = .aai_first_view(pk_a$x, astro$x, peaks_in_view))
  if (!is.null(ties)) {
    ties <- as.data.frame(ties)
    tx <- if ("data_x" %in% names(ties)) { ties$data_x } else { ties[, 1] }
    ta <- if ("insolation_x" %in% names(ties)) { ties$insolation_x } else { ties[, 2] }
    ok <- !is.na(tx) & !is.na(ta)
    ta_ok <- as.numeric(ta[ok])
    ta_snap <- numeric(length(ta_ok))
    for (i in seq_along(ta_ok)) {
      ta_snap[i] <- pk_a$x[which.min(abs(pk_a$x - ta_ok[i]))]
    }
    st$ties <- data.frame(p = as.numeric(tx[ok]), a = ta_snap)
    if (nrow(st$ties) > 0) {
      wide <- function(v, lim, w0) {
        r <- range(v)
        w <- min(max(diff(r) * 1.3, w0), diff(lim))
        .aai_fix_view(mean(r) + c(-0.5, 0.5) * w, lim, w)
      }
      st$view_p <- wide(st$ties$p, range(proxy$x), diff(st$view_p))
      st$view_a <- wide(st$ties$a, range(astro$x), diff(st$view_a))
    }
  }
  list(astro = astro, astro_list = alist, proxy = proxy, pk_a = pk_a, pk_p = pk_p, st = st, pts = pts,
      min_or_max_p = proxy_min_or_max, min_or_max_a = astrosolution_min_or_max)
}

# tie points -> output data frame (values of the records at the tie positions)
.aaa_output <- function(ties, proxy, astro) {
  if (nrow(ties) == 0) {
    return(data.frame(data_x = numeric(0), data_y = numeric(0), insolation_x = numeric(0), insolation_y = numeric(0)))
  }
  out <- data.frame(data_x = ties$p, data_y = stats::approx(proxy$x, proxy$y, ties$p, rule = 2)$y,
      insolation_x = ties$a, insolation_y = stats::approx(astro$x, astro$y, ties$a, rule = 2)$y)
  out <- out[order(out$data_x), ]
  rownames(out) <- NULL
  out
}

# the strongest period of the wavelet spectrum of the proxy record (proxy units)
.aaa_est_period <- function(proxy) {
  dt <- stats::median(diff(proxy$x))
  fallback <- 8 * dt
  tryCatch({
    n <- min(4000, nrow(proxy))
    grid <- seq(min(proxy$x), max(proxy$x), length.out = n)
    reg <- data.frame(x = grid, y = stats::approx(proxy$x, proxy$y, grid)$y)
    d <- grid[2] - grid[1]
    w <- analyse_wavelet(reg, dj = 1 / 20, lowerPeriod = 4 * d, upperPeriod = diff(range(grid)) / 4, verbose = FALSE,
        omega_nr = 6)
    pw <- rowMeans(matrix(as.numeric(w$Power), nrow = w$nr, ncol = w$nc))
    as.numeric(w$Period)[which.max(pw)]
  }, error = function(e) fallback)
}

# Taner filter (astrochron) on a record, interpolated to an even grid; flow/fhigh in cycles per unit, roll = 10^roll_exp
.aaa_taner <- function(dat, flow, fhigh, roll_exp, n_max = 6000) {
  if (!all(is.finite(c(flow, fhigh, roll_exp))) || flow <= 0 || fhigh <= flow) {
    stop("need 0 < flow < fhigh")
  }
  n <- min(n_max, nrow(dat))
  grid <- seq(min(dat$x), max(dat$x), length.out = n)
  reg <- data.frame(x = grid, y = stats::approx(dat$x, dat$y, grid)$y)
  nyq <- 0.5 / (grid[2] - grid[1])
  if (flow >= nyq) {
    stop("flow is above the Nyquist frequency (", signif(nyq, 4), ")")
  }
  fhigh <- min(fhigh, nyq)
  ext <- astrochron::taner(reg, flow = flow, fhigh = fhigh, roll = 10^roll_exp, addmean = FALSE, genplot = FALSE,
      verbose = FALSE, check = FALSE)
  out <- data.frame(x = ext[, 1], y = ext[, 2])
  out[is.finite(out$y), ]
}

# default Taner band around a period: period -/+ 25 percent
.aaa_band <- function(period, width = 0.25) c(flow = 1 / (period * (1 + width)), fhigh = 1 / (period * (1 - width)))

# a cosine metronome: maxima at phase + k * period
.aaa_metronome <- function(period, phase, from, to) {
  if (!all(is.finite(c(period, phase, from, to))) || period <= 0 || to <= from) {
    stop("need period > 0 and from < to")
  }
  dt <- max(period / 24, (to - from) / 50000)
  dt <- period / ceiling(period / dt)                                   # a whole number of steps per period
  x <- phase + dt * seq(ceiling((from - phase) / dt), floor((to - phase) / dt))   # maxima fall exactly on the grid
  data.frame(x = x, y = cos(2 * pi * (x - phase) / period))
}

# Laskar solutions via astrochron::getLaskar (a download); 'getter' can be replaced, e.g. by a cached version
.aaa_laskar_choices <- function() {
  c("Laskar 2004: ETP" = "laskar:la04_etp", "Laskar 2004: eccentricity" = "laskar:la04_ecc",
      "Laskar 2004: obliquity" = "laskar:la04_obl", "Laskar 2004: climatic precession" = "laskar:la04_prec",
      "Laskar 2010a: eccentricity" = "laskar:la10a", "Laskar 2010b: eccentricity" = "laskar:la10b",
      "Laskar 2010c: eccentricity" = "laskar:la10c", "Laskar 2010d: eccentricity" = "laskar:la10d",
      "Laskar 2011: eccentricity" = "laskar:la11", "Laskar 2004: insolation 65 N (summer solstice)" = "laskar:insolation")
}
.aaa_laskar <- function(key, from, to, getter = function(sol) astrochron::getLaskar(sol, verbose = FALSE)) {
  key <- sub("^laskar:", "", key)
  if (!is.finite(from) || !is.finite(to) || to <= from) {
    stop("need from < to")
  }
  if (grepl("^la04_", key)) {
    w <- switch(sub("la04_", "", key), etp = c(1, 1, 1), ecc = c(1, 0, 0), obl = c(0, 1, 0), prec = c(0, 0, 1))
    d <- getter("la04")
    e <- astrochron::etp(tmin = max(0, from), tmax = min(to, max(d[, 1])), dt = 1, eWt = w[1], oWt = w[2], pWt = w[3],
        esinw = TRUE, solution = d, standardize = TRUE, genplot = FALSE, verbose = FALSE)
    out <- data.frame(x = e[, 1], y = e[, 2])
  } else {
    d <- getter(key)
    out <- data.frame(x = d[, 1], y = d[, 2])
    out <- out[out$x >= from & out$x <= to, ]
  }
  out[is.finite(out$y), ]
}

# zoom a window around x by factor f (< 1 zooms in), kept inside lim
.aaa_zoom_at <- function(view, x, f, lim) {
  w <- min(diff(view) * f, diff(lim))
  x <- min(max(x, view[1]), view[2])
  frac <- (x - view[1]) / diff(view)
  .aai_fix_view(c(x - frac * w, x + (1 - frac) * w), lim, w)
}

# Pure tie logic (no shiny)
# state: list(ties = data.frame(p, a) positions in data units, pending = NULL or list(panel, x), sel = NULL or row)

.aaa_snap <- function(x, cands, snap) if (snap && nrow(cands) > 0) { cands$x[which.min(abs(cands$x - x))] } else { x }

# index of the tie point of one record within 'tol' of x (nearest), or integer(0)
.aaa_near <- function(s, panel, x, tol) {
  pos <- s$ties[[if (panel == "proxy") { "p" } else { "a" }]]
  if (length(pos) == 0) {
    return(integer(0))
  }
  k <- which.min(abs(pos - x))
  if (abs(pos[k] - x) <= tol) {
    k
  } else {
    integer(0)
  }
}

# click on position x of record 'panel' ("proxy" or "astro")
.aaa_click <- function(s, panel, x, cands, tol, snap = TRUE) {
  col <- if (panel == "proxy") { "p" } else { "a" }
  hit <- .aaa_near(s, panel, x, tol)
  if (length(hit) > 0) {                      # a tie point: select it (or deselect when it is selected already)
    s$sel <- if (identical(s$sel, hit)) { NULL } else { hit }
    s$pending <- NULL
    return(s)
  }
  xn <- .aaa_snap(x, cands, snap)
  if (!is.null(s$sel)) {                      # free position while a tie point is selected: move that side of the tie
    s$ties[s$sel, col] <- xn
    s$sel <- NULL
    return(s)
  }
  if (is.null(s$pending)) {
    s$pending <- list(panel = panel, x = xn)
  } else if (s$pending$panel == panel) {      # same record again: change the selection
    s$pending <- list(panel = panel, x = xn)
  } else {                                    # other record: tie
    p_x <- if (panel == "proxy") { xn } else { s$pending$x }
    a_x <- if (panel == "astro") { xn } else { s$pending$x }
    s$ties <- rbind(s$ties, data.frame(p = p_x, a = a_x))
    s$pending <- NULL
  }
  s
}

.aaa_untie <- function(s, panel, x, tol) {
  hit <- .aaa_near(s, panel, x, tol)
  if (length(hit) > 0) {
    s$ties <- s$ties[-hit, , drop = FALSE]
  }
  s$sel <- NULL
  s$pending <- NULL
  s
}

# drag from x0 to x1 over one record: the tie point at x0 moves to x1
.aaa_drag <- function(s, panel, x0, x1, cands, tol, snap = TRUE) {
  hit <- .aaa_near(s, panel, x0, tol)
  if (length(hit) == 0) {
    return(s)
  }
  s$ties[hit, if (panel == "proxy") { "p" } else { "a" }] <- .aaa_snap(x1, cands, snap)
  s$sel <- NULL
  s$pending <- NULL
  s
}

# move the selected tie point: to the next peak in direction 'dir' (snap) or by 'dx' data units
.aaa_nudge <- function(s, panel, dx, cands, snap = TRUE) {
  if (is.null(s$sel)) {
    return(s)
  }
  col <- if (panel == "proxy") { "p" } else { "a" }
  cur <- s$ties[s$sel, col]
  if (snap && nrow(cands) > 0) {
    nxt <- if (dx > 0) { cands$x[cands$x > cur + 1e-9] } else { cands$x[cands$x < cur - 1e-9] }
    if (length(nxt) == 0) {
      return(s)
    }
    s$ties[s$sel, col] <- if (dx > 0) { min(nxt) } else { max(nxt) }
  } else {
    s$ties[s$sel, col] <- cur + dx
  }
  s
}

# tie points in the wrong order (the tuned time would run backwards)
.aaa_flag <- function(ties, time_dir) {
  n <- nrow(ties)
  flag <- rep(FALSE, n)
  if (n < 2) {
    return(flag)
  }
  o <- order(ties$p)
  d <- diff(ties$a[o])
  bad <- if (time_dir) { d <= 0 } else { d >= 0 }
  fo <- rep(FALSE, n)
  fo[-1] <- fo[-1] | bad
  fo[-n] <- fo[-n] | bad
  flag[o] <- fo
  flag
}


.aaa_cols <- function(n) {
  if (n == 0) {
    return(character(0))
  }
  grDevices::hsv(h = (seq_len(n) * 0.618034) %% 1, s = 0.75, v = 0.75)
}

# Tuning

.aaa_tune <- function(proxy, astro, pk_a, ties, time_dir, flow, fhigh, roll_exp, sedrate_factor) {
  if (nrow(ties) < 2) {
    stop("tie at least two points")
  }
  t_ord <- ties[order(ties$data_x), ]
  d <- diff(t_ord$insolation_x)
  if (any(if (time_dir) { d <= 0 } else { d >= 0 })) {
    stop("tie points in the wrong order (black crosses)")
  }
  tuned <- as.data.frame(anchor2time(anchor_points = cbind(t_ord$data_x, t_ord$insolation_x), data = proxy))
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
  if (is.null(flow) || is.null(fhigh) || !is.finite(flow) || !is.finite(fhigh)) {
    b <- .aaa_band(stats::median(diff(pk_a$x)))
    flow <- b[["flow"]]
    fhigh <- b[["fhigh"]]
  }
  ext <- .aaa_taner(data.frame(x = reg$time, y = reg$value), flow, fhigh, roll_exp, n_max = 20000)
  a_ext <- stats::approx(astro$x, astro$y, ext$x)$y
  r <- suppressWarnings(stats::cor(ext$y, a_ext, use = "complete.obs"))
  list(reg = reg, ext = as.matrix(ext), r = r, band = c(flow, fhigh), roll = roll_exp, t_ord = t_ord,
      sedrate = diff(t_ord$data_x) / diff(t_ord$insolation_x) * sedrate_factor)
}

# Geometry (normalised device coordinates of the single plot)

.aaa_boxes <- function() {
  list(proxy = c(0.075, 0.985, 0.745, 0.975), astro = c(0.075, 0.985, 0.405, 0.635), tuned = c(0.075, 0.985, 0.045, 0.285))
}

.aaa_yr <- function(dat, view) {
  inside <- dat$x >= view[1] & dat$x <= view[2]
  yr <- range(dat$y[inside], na.rm = TRUE)
  if (!is.finite(diff(yr)) || diff(yr) == 0) {
    yr <- range(dat$y)
  }
  yr + c(-0.08, 0.14) * diff(yr)
}

# data x <-> normalised device x for one record
.aaa_to_ndc <- function(box, view, reversed, x) {
  f <- (x - view[1]) / (view[2] - view[1])
  if (reversed) {
    f <- 1 - f
  }
  box[1] + f * (box[2] - box[1])
}
.aaa_from_ndc <- function(box, view, reversed, xn) {
  f <- (min(max(xn, box[1]), box[2]) - box[1]) / (box[2] - box[1])
  if (reversed) {
    f <- 1 - f
  }
  view[1] + f * (view[2] - view[1])
}
# which record a vertical position (ndc) belongs to
.aaa_panel_at <- function(y) {
  b <- .aaa_boxes()
  for (nm in c("proxy", "astro")) {
    if (y >= b[[nm]][3] - 0.035 && y <= b[[nm]][4] + 0.02) {
      return(nm)
    }
  }
  NULL
}


.aaa_zs <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)


.aaa_draw_record <- function(dat, cands, view, s, panel, flag, reversed, ylab, box, new, snap) {
  col <- if (panel == "proxy") { "p" } else { "a" }
  other <- if (panel == "proxy") { "a" } else { "p" }
  graphics::par(fig = c(0, 1, 0, 1), plt = box, new = new, mgp = c(2.4, 0.7, 0), cex.axis = 0.85)
  yr <- .aaa_yr(dat, view)
  xl <- if (reversed) { rev(view) } else { view }
  graphics::plot(dat$x, dat$y, type = "n", xlim = xl, ylim = yr, xaxs = "i", yaxs = "i", xlab = "", ylab = ylab)
  graphics::rect(graphics::par("usr")[1], yr[1], graphics::par("usr")[2], yr[2], col = "grey98", border = NA)
  graphics::grid(col = "grey88", lty = 1)
  graphics::lines(dat$x, dat$y, col = "grey25", lwd = 1.2)
  graphics::box(col = "grey50")
  if (snap) {
    graphics::points(cands$x, cands$y, pch = 21, bg = "white", col = "grey45", cex = 1.1)
  }
  n <- nrow(s$ties)
  if (n > 0) {
    ord <- rank(s$ties$p, ties.method = "first")             # numbered along the proxy record
    cols <- .aaa_cols(n)[ord]
    xv <- s$ties[[col]]
    yv <- stats::approx(dat$x, dat$y, xv, rule = 2)$y
    vis <- xv >= view[1] & xv <= view[2]
    if (any(vis)) {
      graphics::segments(xv[vis], yr[1], xv[vis], yv[vis], col = cols[vis], lty = 3, lwd = 1.5)
      graphics::points(xv[vis], yv[vis], pch = 21, bg = cols[vis], col = "grey15", cex = 1.9, lwd = 1.2)
      graphics::text(xv[vis], yv[vis], ord[vis], col = "white", cex = 0.7, font = 2)
      graphics::text(xv[vis], rep(yr[1], sum(vis)), paste0(" ", signif(s$ties[[other]][vis], 5)), srt = 90, adj = c(0, -0.4),
          cex = 0.65, col = cols[vis])
    }
    if (any(flag & vis)) {
      graphics::points(xv[flag & vis], yv[flag & vis], pch = 4, cex = 2.8, lwd = 2.5, col = "black")
    }
    if (!is.null(s$sel) && vis[s$sel]) {
      graphics::points(xv[s$sel], yv[s$sel], pch = 1, cex = 3.6, lwd = 3, col = "red")
    }
  }
  if (!is.null(s$pending) && s$pending$panel == panel) {
    yp <- stats::approx(dat$x, dat$y, s$pending$x, rule = 2)$y
    graphics::points(s$pending$x, yp, pch = 21, bg = "red", col = "darkred", cex = 2.1)
  }
  invisible(NULL)
}

# overlay: tie lines between the two records (last plot: its coordinates are the click coordinates)
.aaa_draw_links <- function(s, view_p, view_a, rev_p, flag, rev_a = FALSE) {
  b <- .aaa_boxes()
  graphics::par(fig = c(0, 1, 0, 1), plt = c(0, 1, 0, 1), new = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  n <- nrow(s$ties)
  if (n == 0) {
    return(invisible(NULL))
  }
  ord <- rank(s$ties$p, ties.method = "first")
  cols <- .aaa_cols(n)[ord]
  y1 <- b$proxy[3]
  y2 <- b$astro[4]
  for (k in seq_len(n)) {
    xp <- s$ties$p[k]
    xa <- s$ties$a[k]
    p_in <- xp >= view_p[1] && xp <= view_p[2]
    a_in <- xa >= view_a[1] && xa <= view_a[2]
    lw <- if (identical(s$sel, k)) { 4 } else { 2.2 }
    col <- if (flag[k]) { "black" } else { cols[k] }
    if (p_in && a_in) {
      graphics::segments(.aaa_to_ndc(b$proxy, view_p, rev_p, xp), y1, .aaa_to_ndc(b$astro, view_a, rev_a, xa), y2, lwd = lw, col = col)
    } else if (p_in) {
      xn <- .aaa_to_ndc(b$proxy, view_p, rev_p, xp)
      graphics::segments(xn, y1, xn, y1 - 0.05, lwd = lw, col = col)
      graphics::text(xn, y1 - 0.062, paste0("-> ", signif(xa, 5)), cex = 0.7, col = col)
    } else if (a_in) {
      xn <- .aaa_to_ndc(b$astro, view_a, rev_a, xa)
      graphics::segments(xn, y2, xn, y2 + 0.05, lwd = lw, col = col)
      graphics::text(xn, y2 + 0.062, paste0("-> ", signif(xp, 5)), cex = 0.7, col = col)
    }
  }
  invisible(NULL)
}


.aaa_draw_tuned <- function(tu, astro, s, pk_a, show_astro, box) {
  graphics::par(fig = c(0, 1, 0, 1), plt = box, new = TRUE, mgp = c(2.4, 0.7, 0), cex.axis = 0.85)
  if (!is.null(tu$error)) {
    graphics::plot.new()
    graphics::box(col = "grey70")
    graphics::text(0.5, 0.5, paste0("Tuned record: ", tu$error), col = "grey40", cex = 1.1)
    return(invisible(NULL))
  }
  reg <- tu$reg
  ext <- tu$ext
  xr <- range(reg$time)
  yz <- range(c(.aaa_zs(reg$value), .aaa_zs(ext[, 2])), na.rm = TRUE)
  ai <- astro[astro$x >= xr[1] & astro$x <= xr[2], ]
  if (show_astro && nrow(ai) > 1) {
    yz <- range(c(yz, .aaa_zs(ai$y)), na.rm = TRUE)
  }
  graphics::plot(reg$time, .aaa_zs(reg$value), type = "n", xlim = xr, ylim = yz, xaxs = "i", xlab = "", ylab = "z-score",
      main = sprintf("Tuned record, Taner band %s - %s%s   (r with the astronomical solution = %s)",
          signif(tu$band[1], 3), signif(tu$band[2], 3), if (show_astro) { " and astronomical solution" } else { "" },
          format(round(tu$r, 2), nsmall = 2)),
      cex.main = 0.95, font.main = 1, adj = 0)
  graphics::grid(col = "grey88", lty = 1)
  tie_cols <- .aaa_cols(nrow(tu$t_ord))
  if (length(tie_cols) > 0) {
    graphics::abline(v = tu$t_ord$insolation_x, col = grDevices::adjustcolor(tie_cols, 0.8), lty = 3, lwd = 1.5)
  }
  graphics::lines(reg$time, .aaa_zs(reg$value), col = "grey35")
  if (show_astro && nrow(ai) > 1) {
    graphics::lines(ai$x, .aaa_zs(ai$y), col = grDevices::adjustcolor("dodgerblue3", 0.8), lwd = 1.6)
  }
  graphics::lines(ext[, 1], .aaa_zs(ext[, 2]), col = "red", lwd = 2.2)
  graphics::legend("topright", horiz = TRUE, bty = "n", cex = 0.8, lty = 1, lwd = c(1, 2.2, 1.6),
      col = c("grey35", "red", "dodgerblue3")[c(TRUE, TRUE, show_astro)],
      legend = c("tuned proxy", "extracted cycle", "astronomical solution")[c(TRUE, TRUE, show_astro)])
  graphics::box(col = "grey50")
}

.aaa_draw_sedrate <- function(tu) {
  graphics::par(mar = c(3.4, 4.6, 0.6, 1), mgp = c(2.4, 0.7, 0), cex.axis = 0.85)
  if (!is.null(tu$error)) {
    graphics::plot.new()
    return(invisible(NULL))
  }
  t <- tu$t_ord
  sr <- tu$sedrate
  xr <- range(tu$reg$time)
  graphics::plot(range(t$insolation_x), range(c(0, sr), na.rm = TRUE), type = "n", xlim = xr, xaxs = "i", xlab = "time",
      ylab = "sedimentation rate")
  graphics::grid(col = "grey88", lty = 1)
  graphics::segments(t$insolation_x[-nrow(t)], sr, t$insolation_x[-1], sr, lwd = 4, col = "darkorange3", lend = 1)
  tie_cols <- .aaa_cols(nrow(t))
  graphics::abline(v = t$insolation_x, col = grDevices::adjustcolor(tie_cols, 0.8), lty = 3, lwd = 1.5)
  graphics::box(col = "grey50")
}


# The app

.aaa_app <- function(S, time_dir, sedrate_factor, tie_on = "raw", proxy_tan = NULL, tuned_tan = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("package 'shiny' is needed", call. = FALSE)
  }
  astro <- S$astro
  proxy <- S$proxy
  pk_a <- S$pk_a
  pk_p <- S$pk_p
  # Taner filter settings: c(flow, fhigh, roll exponent)
  tuned_tan <- c(tuned_tan, NA, NA, 3)[seq_len(3)] * 1
  if (length(proxy_tan) < 3) {
    b <- .aaa_band(.aaa_est_period(proxy))
    proxy_tan <- c(b[["flow"]], b[["fhigh"]], if (length(proxy_tan) == 1) { proxy_tan } else { 3 })
  }
  if (length(tuned_tan) < 3 || any(is.na(tuned_tan[1:2]))) {
    b <- .aaa_band(stats::median(diff(pk_a$x)))
    tuned_tan <- c(b[["flow"]], b[["fhigh"]], if (is.finite(tuned_tan[3])) { tuned_tan[3] } else { 3 })
  }
  proxy_tan <- signif(proxy_tan, 4)
  tuned_tan <- signif(tuned_tan, 4)

  css <- "
    body { background:#f4f5f7; }
    .aaa-title { font-weight:600; font-size:20px; margin:10px 0 2px 0; }
    .aaa-sub { color:#667; margin-bottom:10px; }
    .aaa-card { background:#fff; border:1px solid #dde; border-radius:8px; padding:10px 12px; margin-bottom:10px; }
    .aaa-card h4 { margin:0 0 8px 0; font-size:13px; text-transform:uppercase; letter-spacing:.05em; color:#667; }
    .aaa-status { background:#eef3ff; border-left:4px solid #3b6fe0; padding:6px 10px; margin-bottom:8px; border-radius:4px; }
    .aaa-row { margin-bottom:4px; } .aaa-btn { margin:1px 0; }
    .aaa-plot { background:#fff; border:1px solid #dde; border-radius:8px; margin-bottom:8px; padding:2px; }
    .aaa-label { font-size:12px; color:#667; margin:0 0 2px 4px; }
    .aaa-card .form-group { margin-bottom:6px; }
    table.aaa-ties td, table.aaa-ties th { padding:1px 6px; font-size:12px; }
  "
  js <- "
  $(document).on('shiny:connected', function() {
    var el = document.getElementById('plot_main');
    if (!el || el._aaa) return;
    el._aaa = true;
    function rel(e) {
      var img = el.querySelector('img');
      if (!img) return null;
      var r = img.getBoundingClientRect();
      return {x: (e.clientX - r.left) / r.width, y: 1 - (e.clientY - r.top) / r.height};
    }
    el.addEventListener('wheel', function(e) {
      var p = rel(e); if (!p) return;
      e.preventDefault();
      Shiny.setInputValue('wheel_main', {x: p.x, y: p.y, dy: e.deltaY, t: Date.now()}, {priority: 'event'});
    }, {passive: false});
    var drag = null;
    el.addEventListener('mousedown', function(e) {
      if (e.shiftKey || e.button === 1 || e.button === 2) {
        var p = rel(e); if (!p) return;
        drag = p; e.preventDefault(); e.stopPropagation();
      }
    }, true);
    window.addEventListener('mousemove', function(e) {
      if (!drag) return;
      var p = rel(e); if (!p) return;
      Shiny.setInputValue('pan_main', {x: drag.x, y: drag.y, dx: p.x - drag.x, t: Date.now()}, {priority: 'event'});
      drag = p;
    });
    window.addEventListener('mouseup', function() { drag = null; });
    el.addEventListener('contextmenu', function(e) { e.preventDefault(); });
  });
  "
  btn <- function(id, label, ...) shiny::actionButton(id, label, class = "btn-sm aaa-btn", ...)
  numbox <- function(id, label) shiny::div(style = "display:inline-block; width:48%;",
      shiny::numericInput(id, label, NA, step = NA, width = "100%"))
  ui <- shiny::fluidPage(
      shiny::tags$head(shiny::tags$style(shiny::HTML(css))),
      shiny::div(class = "aaa-title", "Tie a proxy record to an astronomical solution"),
      shiny::div(class = "aaa-sub", "Click each record to tie (peaks or anywhere, see 'Snap to peaks'). Click a tie point to select it, ",
          "then use the arrows, type a position, click elsewhere or drag over the record to move it. Double click to untie."),
      shiny::fluidRow(
          shiny::column(3,
              shiny::div(class = "aaa-card", shiny::h4("Astronomical solution"),
                  shiny::selectInput("a_src", NULL, choices = c(stats::setNames(paste0("prov:", names(S$astro_list)), names(S$astro_list)),
                      "Metronome (cosine)" = "metro", .aaa_laskar_choices()), selected = paste0("prov:", names(S$astro_list)[1])),
                  shiny::conditionalPanel("input.a_src == 'metro'",
                      shiny::div(style = "display:flex; gap:6px;",
                          shiny::numericInput("a_period", "period", 405, min = 0, step = NA, width = "50%"),
                          shiny::numericInput("a_phase", "time of a maximum", 0, step = NA, width = "50%")),
                      shiny::selectInput("a_preset", "preset period", c("405 (long eccentricity)" = 405, "100 (short eccentricity)" = 100,
                          "41 (obliquity)" = 41, "23 (precession)" = 23, "19 (precession)" = 19))),
                  shiny::div(style = "display:flex; gap:6px;",
                      shiny::numericInput("a_from", "from (time)", signif(min(S$astro$x), 8), step = NA, width = "50%"),
                      shiny::numericInput("a_to", "to (time)", signif(max(S$astro$x), 8), step = NA, width = "50%")),
                  shiny::conditionalPanel("input.a_src.indexOf('laskar:') === 0",
                      shiny::actionButton("a_load", "Download and load", class = "btn-sm btn-primary"),
                      shiny::div(shiny::tags$small("Downloads the solution over the internet (astrochron::getLaskar); time in ka."))),
                  shiny::radioButtons("pk_a_type", "Tie to", c("maxima" = "max", "minima" = "min"), inline = TRUE,
                      selected = S$min_or_max_a)),
              shiny::div(class = "aaa-card", shiny::h4("Tie on"),
                  shiny::radioButtons("tie_on", NULL, c("Raw proxy record" = "raw", "Cycle extracted from the proxy" = "cycle"),
                      selected = tie_on),
                  shiny::conditionalPanel("input.tie_on == 'cycle'",
                      shiny::tags$small("Taner filter (cycles per proxy unit)"),
                      shiny::div(style = "display:flex; gap:6px;",
                          shiny::numericInput("pflow", "flow", proxy_tan[1], min = 0, step = NA, width = "33%"),
                          shiny::numericInput("pfhigh", "fhigh", proxy_tan[2], min = 0, step = NA, width = "33%"),
                          shiny::numericInput("proll", "roll = 10^x", proxy_tan[3], step = 0.5, width = "33%"))),
                  shiny::radioButtons("pk_p_type", "Peaks of the proxy", c("maxima" = "max", "minima" = "min"), inline = TRUE,
                      selected = S$min_or_max_p),
                  shiny::checkboxInput("snap", "Snap proxy ties to peaks", tie_on == "cycle"),
                  shiny::tags$small("Ties on the astronomical solution always sit on a maximum/minimum.")),
              shiny::div(class = "aaa-card", shiny::h4("Tie points"),
                  shiny::uiOutput("status"),
                  shiny::div(btn("undo", "Undo"), btn("untie", "Untie selected"), btn("clear", "Clear all")),
                  shiny::div(class = "aaa-row", shiny::tags$small("Selected tie point: move left / right")),
                  shiny::div(shiny::tags$small("proxy "), btn("p_cl", "<<"), btn("p_fl", "<"), btn("p_fr", ">"), btn("p_cr", ">>")),
                  shiny::div(shiny::tags$small("astro "), btn("a_cl", "<<"), btn("a_fl", "<"), btn("a_fr", ">"), btn("a_cr", ">>")),
                  shiny::div(numbox("sel_p", "proxy position"), numbox("sel_a", "astro position (peak)")),
                  shiny::div(numbox("sel_off", "proxy offset from nearest peak")),
                  shiny::tableOutput("ties_table")),
              shiny::div(class = "aaa-card", shiny::h4("Axes and reference frame"),
                  shiny::checkboxInput("time_dir", "Time increases with increasing proxy value (depth/height)", time_dir),
                  shiny::checkboxInput("flip_p", "Flip the proxy axis (right to left)", !time_dir),
                  shiny::checkboxInput("flip_a", "Flip the astronomical axis (right to left)", FALSE)),
              shiny::div(class = "aaa-card", shiny::h4("Tuned record"),
                  shiny::checkboxInput("live", "Update while tying", TRUE),
                  shiny::div(shiny::actionButton("tune", "Tune now", class = "btn-sm btn-primary")),
                  shiny::tags$small("Taner filter of the tuned record (cycles per time unit)"),
                  shiny::div(style = "display:flex; gap:6px;",
                      shiny::numericInput("flow", "flow", tuned_tan[1], min = 0, step = NA, width = "33%"),
                      shiny::numericInput("fhigh", "fhigh", tuned_tan[2], min = 0, step = NA, width = "33%"),
                      shiny::numericInput("roll", "roll = 10^x", tuned_tan[3], step = 0.5, width = "33%")),
                  shiny::checkboxInput("show_astro", "Show astronomical solution", TRUE),
                  shiny::checkboxInput("show_sed", "Show sedimentation rate", TRUE)),
              shiny::actionButton("done", "Done - return tie points", class = "btn-success", width = "100%")),
          shiny::column(9,
              shiny::tags$script(shiny::HTML(js)),
              shiny::fluidRow(
                  shiny::column(6, shiny::sliderInput("win_p", "Proxy window", min = min(proxy$x), max = max(proxy$x), value = S$st$view_p,
                      step = signif(diff(range(proxy$x)) / 2000, 2), width = "100%", ticks = FALSE, sep = "")),
                  shiny::column(6, shiny::sliderInput("win_a", "Astronomical window", min = min(S$astro$x), max = max(S$astro$x),
                      value = S$st$view_a, step = signif(diff(range(S$astro$x)) / 2000, 2), width = "100%",
                      ticks = FALSE, sep = ""))),
              shiny::fluidRow(
                  shiny::column(6, shiny::div(style = "display:flex; gap:6px; align-items:flex-end;",
                      shiny::numericInput("wp_from", "Proxy window: from (type)", signif(S$st$view_p[1], 8), step = NA, width = "50%"),
                      shiny::numericInput("wp_to", "to", signif(S$st$view_p[2], 8), step = NA, width = "50%"))),
                  shiny::column(6, shiny::div(style = "display:flex; gap:6px; align-items:flex-end;",
                      shiny::numericInput("wa_from", "Astronomical window: from (type)", signif(S$st$view_a[1], 8), step = NA, width = "50%"),
                      shiny::numericInput("wa_to", "to", signif(S$st$view_a[2], 8), step = NA, width = "50%")))),
              shiny::div(btn("sync", "Sync astro window to the proxy ties"), btn("fit", "Show both records in full")),
              shiny::div(class = "aaa-label", "Mouse: wheel = zoom the record under the pointer, Shift + drag (or right-drag) = move it left/right, plain drag from a tie point = move the tie point"),
              shiny::div(class = "aaa-plot", shiny::plotOutput("plot_main", height = "880px", click = "click_main", dblclick = shiny::dblclickOpts("dbl_main", delay = 250),
                  brush = shiny::brushOpts("brush_main", direction = "x", delay = 500, delayType = "debounce",
                      fill = "#3b6fe0", opacity = 0.15))),
              shiny::uiOutput("sed_ui"))))

  server <- function(input, output, session) {
    s0 <- S$st
    # reference frame: does time increase with the proxy value (tuning and the order check), and how are the axes drawn
    td <- function() if (is.null(input$time_dir)) { time_dir } else { isTRUE(input$time_dir) }
    flip_p <- function() if (is.null(input$flip_p)) { !time_dir } else { isTRUE(input$flip_p) }
    flip_a <- function() isTRUE(input$flip_a)
    rv <- shiny::reactiveValues(ties = s0$ties, pending = NULL, sel = NULL, history = list(),
        view_p = s0$view_p, view_a = s0$view_a)
    state <- function() list(ties = rv$ties, pending = rv$pending, sel = rv$sel)
    set_state <- function(s, push = TRUE) {
      if (push) {
        rv$history <- c(rv$history, list(list(ties = rv$ties, pending = rv$pending, sel = rv$sel)))
      }
      rv$ties <- s$ties
      rv$pending <- s$pending
      rv$sel <- s$sel
    }
    # The astronomical target: provided solution, metronome or Laskar solution, its time interval and its peaks
    raw_a <- shiny::reactiveVal(S$astro_list[[1]])
    lcache <- new.env()
    getter <- function(sol) {
      if (is.null(lcache[[sol]])) {
        lcache[[sol]] <- astrochron::getLaskar(sol, verbose = FALSE)
      }
      lcache[[sol]]
    }
    astro_r <- shiny::reactive({
      d <- raw_a()
      lo <- input$a_from
      hi <- input$a_to
      if (is.null(lo) || is.null(hi) || !is.finite(lo) || !is.finite(hi)) {
        return(d)
      }
      s <- d[d$x >= lo & d$x <= hi, ]
      if (nrow(s) < 20) {
        d
      } else {
        s
      }
    })
    type_a <- function() if (is.null(input$pk_a_type)) { S$min_or_max_a } else { input$pk_a_type }
    type_p <- function() if (is.null(input$pk_p_type)) { S$min_or_max_p } else { input$pk_p_type }
    pka <- shiny::reactive(.aaa_peaks(astro_r(), type_a(), S$pts))
    pkp <- shiny::reactive(.aaa_peaks(proxy, type_p(), S$pts))
    ra <- function() range(astro_r()$x)
    set_interval <- function(d) {
      shiny::updateNumericInput(session, "a_from", value = signif(min(d$x), 8))
      shiny::updateNumericInput(session, "a_to", value = signif(max(d$x), 8))
    }
    build_metronome <- function() {
      d <- tryCatch(.aaa_metronome(input$a_period, input$a_phase, input$a_from, input$a_to),
          error = function(e) {
            shiny::showNotification(paste("Metronome:", conditionMessage(e)), type = "error")
            NULL
          })
      if (!is.null(d)) {
        raw_a(d)
      }
    }
    src_changed <- FALSE
    shiny::observeEvent(input$a_src, {
      src <- input$a_src
      src_changed <<- TRUE
      if (startsWith(src, "prov:")) {
        d <- S$astro_list[[sub("^prov:", "", src)]]
        raw_a(d)
        set_interval(d)
      } else if (src == "metro") {
        build_metronome()
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(list(input$a_period, input$a_phase, input$a_from, input$a_to), {
      if (identical(input$a_src, "metro")) {
        build_metronome()
      }
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$a_preset, shiny::updateNumericInput(session, "a_period", value = as.numeric(input$a_preset)),
        ignoreInit = TRUE)
    shiny::observeEvent(input$a_load, {
      key <- input$a_src
      shiny::withProgress(message = "Downloading / computing the solution ...", value = 0.5, {
        d <- tryCatch(.aaa_laskar(key, input$a_from, input$a_to, getter),
            error = function(e) {
              shiny::showNotification(paste("Could not load the solution:", conditionMessage(e)), type = "error", duration = 10)
              NULL
            })
      })
      if (!is.null(d) && nrow(d) > 20) {
        raw_a(d)
      }
    })
    # a new kind of target removes the tie points; changing only the metronome or the interval moves the astronomical side
    # of the tie points to the nearest peaks when that keeps them apart and inside the new interval
    shiny::observeEvent(raw_a(), {
      d <- raw_a()
      if (nrow(rv$ties) > 0) {
        s <- state()
        pk <- .aaa_peaks(d, type_a(), S$pts)
        new_a <- numeric(length(s$ties$a))
        for (i in seq_along(s$ties$a)) {
          new_a[i] <- pk$x[which.min(abs(pk$x - s$ties$a[i]))]
        }
        if (!src_changed && all(s$ties$a >= min(d$x) & s$ties$a <= max(d$x)) && !anyDuplicated(new_a)) {
          s$ties$a <- new_a
          shiny::showNotification("Astronomical tie points moved to the nearest peaks of the new target.", duration = 5)
        } else {
          s$ties <- s$ties[0, ]
          s$sel <- NULL
          s$pending <- NULL
          shiny::showNotification("Tie points removed: they do not fit the new astronomical target.", duration = 6)
        }
        set_state(s)
      }
      src_changed <<- FALSE
      rv$view_a <- .aai_first_view(.aaa_peaks(d, type_a(), S$pts)$x, d$x, 12)
    }, ignoreInit = TRUE)
    shiny::observeEvent(list(input$a_from, input$a_to, input$pk_a_type), {
      r <- ra()
      rv$view_a <- .aai_fix_view(rv$view_a, r, min(diff(rv$view_a), diff(r)))
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$pk_a_type, {      # maxima <-> minima: re-snap the astronomical tie points
      if (nrow(rv$ties) > 0) {
        s <- state()
        pk <- pka()
        new_a <- numeric(length(s$ties$a))
        for (i in seq_along(s$ties$a)) {
          new_a[i] <- pk$x[which.min(abs(pk$x - s$ties$a[i]))]
        }
        s$ties$a <- new_a
        set_state(s)
      }
    }, ignoreInit = TRUE)
    # visible windows: number boxes, mouse wheel (zoom) and Shift + drag (pan)
    close_to <- function(a, b) all(abs(a - b) <= 1e-6 * max(1, abs(b)))
    shiny::observe({
      shiny::updateNumericInput(session, "wp_from", value = signif(rv$view_p[1], 8))
      shiny::updateNumericInput(session, "wp_to", value = signif(rv$view_p[2], 8))
    })
    shiny::observe({
      shiny::updateNumericInput(session, "wa_from", value = signif(rv$view_a[1], 8))
      shiny::updateNumericInput(session, "wa_to", value = signif(rv$view_a[2], 8))
    })
    shiny::observe({                     # sliders follow the windows (the astronomical range follows the target)
      r <- ra()
      v <- pmin(pmax(rv$view_a, r[1]), r[2])
      shiny::updateSliderInput(session, "win_a", min = r[1], max = r[2], value = v, step = signif(diff(r) / 2000, 2))
    })
    shiny::observe(shiny::updateSliderInput(session, "win_p", value = rv$view_p))
    shiny::observeEvent(input$win_p, set_window("p", input$win_p[1], input$win_p[2]))
    shiny::observeEvent(input$win_a, set_window("a", input$win_a[1], input$win_a[2]))
    set_window <- function(which, lo, hi) {
      if (is.null(lo) || is.null(hi) || !is.finite(lo) || !is.finite(hi) || hi <= lo) {
        return()
      }
      lim <- if (which == "p") { range(proxy$x) } else { ra() }
      v <- c(max(lo, lim[1]), min(hi, lim[2]))
      if (v[2] <= v[1]) {
        return()
      }
      cur <- if (which == "p") { shiny::isolate(rv$view_p) } else { shiny::isolate(rv$view_a) }
      if (close_to(v, cur)) {
        return()
      }
      if (which == "p") {
        rv$view_p <- v
      } else {
        rv$view_a <- v
      }
    }
    shiny::observeEvent(list(input$wp_from, input$wp_to), set_window("p", input$wp_from, input$wp_to), ignoreInit = TRUE)
    shiny::observeEvent(list(input$wa_from, input$wa_to), set_window("a", input$wa_from, input$wa_to), ignoreInit = TRUE)
    panel_for_y <- function(y) {          # the record whose box (plus margins) holds a vertical position
      b <- .aaa_boxes()
      for (nm in c("proxy", "astro")) {
        if (y >= b[[nm]][3] - 0.07 && y <= b[[nm]][4] + 0.02) {
          return(nm)
        }
      }
      NULL
    }
    shiny::observeEvent(input$wheel_main, {
      w <- input$wheel_main
      panel <- panel_for_y(w$y)
      if (is.null(panel)) {
        return()
      }
      b <- .aaa_boxes()
      f <- if (w$dy > 0) { 1.25 } else { 0.8 }
      if (panel == "proxy") {
        xc <- .aaa_from_ndc(b$proxy, rv$view_p, flip_p(), w$x)
        rv$view_p <- .aaa_zoom_at(rv$view_p, xc, f, range(proxy$x))
      } else {
        xc <- .aaa_from_ndc(b$astro, rv$view_a, flip_a(), w$x)
        rv$view_a <- .aaa_zoom_at(rv$view_a, xc, f, ra())
      }
    })
    shiny::observeEvent(input$pan_main, {
      m <- input$pan_main
      panel <- panel_for_y(m$y)
      if (is.null(panel)) {
        return()
      }
      b <- .aaa_boxes()
      if (panel == "proxy") {
        w <- diff(rv$view_p)
        shift <- -m$dx / (b$proxy[2] - b$proxy[1]) * w * (if (flip_p()) { -1 } else { 1 })
        rv$view_p <- .aai_fix_view(rv$view_p + shift, range(proxy$x), w)
      } else {
        w <- diff(rv$view_a)
        rv$view_a <- .aai_fix_view(rv$view_a - m$dx / (b$astro[2] - b$astro[1]) * w * (if (flip_a()) { -1 } else { 1 }), ra(), w)
      }
    })

    out <- function() .aaa_output(rv$ties, proxy, astro_r())

    # what is shown and tied on the proxy side: the raw record or a cycle extracted from it
    cyc <- shiny::reactive({
      if (!identical(input$tie_on, "cycle")) {
        return(NULL)
      }
      tryCatch(list(dat = .aaa_taner(proxy, input$pflow, input$pfhigh, input$proll)),
          error = function(e) list(error = conditionMessage(e)))
    })
    rec_p <- shiny::reactive({
      cy <- cyc()
      if (is.null(cy) || !is.null(cy$error)) {
        return(list(dat = proxy, cands = pkp()))
      }
      list(dat = cy$dat, cands = .aaa_peaks(cy$dat, type_p(), S$pts))
    })
    snap <- function(panel = "proxy") if (panel == "astro") { TRUE } else { isTRUE(input$snap) }   # astro ties: always on peaks
    shiny::observeEvent(input$tie_on, shiny::updateCheckboxInput(session, "snap", value = identical(input$tie_on, "cycle")),
        ignoreInit = TRUE)

    geo <- function() {
      b <- .aaa_boxes()
      list(proxy = list(box = b$proxy, view = rv$view_p, rev = flip_p()),
          astro = list(box = b$astro, view = rv$view_a, rev = flip_a()))
    }
    cands_of <- function(panel) if (panel == "proxy") { rec_p()$cands } else { pka() }
    # map a plot position (normalised device coordinates) to a record and a data value
    locate <- function(xn, yn) {
      panel <- .aaa_panel_at(yn)
      if (is.null(panel)) {
        return(NULL)
      }
      g <- geo()[[panel]]
      list(panel = panel, x = .aaa_from_ndc(g$box, g$view, g$rev, xn), tol = 0.025 * diff(g$view))
    }
    last_brush <- 0
    shiny::observeEvent(input$click_main, {
      if (as.numeric(Sys.time()) - last_brush < 1) {
        return()
      }      # the mouse-up of a drag is not a click
      h <- locate(input$click_main$x, input$click_main$y)
      if (!is.null(h)) {
        set_state(.aaa_click(shiny::isolate(state()), h$panel, h$x, cands_of(h$panel), h$tol, snap(h$panel)))
      }
    })
    shiny::observeEvent(input$dbl_main, {
      h <- locate(input$dbl_main$x, input$dbl_main$y)
      if (is.null(h)) {
        return()
      }
      if (length(.aaa_near(state(), h$panel, h$x, h$tol)) > 0) {
        set_state(.aaa_untie(state(), h$panel, h$x, h$tol))
      }
    })
    shiny::observeEvent(input$brush_main, {
      b <- input$brush_main
      if (is.null(b)) {
        return()
      }
      last_brush <<- as.numeric(Sys.time())
      h0 <- locate(b$xmin, (b$ymin + b$ymax) / 2)
      h1 <- locate(b$xmax, (b$ymin + b$ymax) / 2)
      if (!is.null(h0) && !is.null(h1) && nrow(rv$ties) > 0) {
        s_old <- state()
        pos <- s_old$ties[[if (h0$panel == "proxy") { "p" } else { "a" }]]
        # the tie point nearest to either edge is where the drag started
        start_low <- min(abs(pos - h0$x)) <= min(abs(pos - h1$x))
        x_from <- if (start_low) { h0$x } else { h1$x }
        x_to <- if (start_low) { h1$x } else { h0$x }
        new <- .aaa_drag(s_old, h0$panel, x_from, x_to, cands_of(h0$panel), 3 * h0$tol, snap(h0$panel))
        if (!identical(new$ties, s_old$ties)) {
          set_state(new)
        }
      }
      session$resetBrush("brush_main")
    })

    shiny::observeEvent(input$undo, {
      if (length(rv$history) > 0) {
        last <- rv$history[[length(rv$history)]]
        rv$history <- rv$history[-length(rv$history)]
        rv$ties <- last$ties
        rv$pending <- last$pending
        rv$sel <- last$sel
      }
    })
    shiny::observeEvent(input$clear, {
      s <- state()
      s$ties <- s$ties[0, ]
      s$pending <- NULL
      s$sel <- NULL
      set_state(s)
    })
    shiny::observeEvent(input$untie, {
      if (!is.null(rv$sel)) {
        s <- state()
        s$ties <- s$ties[-s$sel, , drop = FALSE]
        s$sel <- NULL
        set_state(s)
      }
    })
    # arrows: fine and coarse steps (or the next peak when snapping)
    nudge <- function(panel, dir, coarse) {
      v <- if (panel == "proxy") { rv$view_p } else { rv$view_a }
      dx <- dir * diff(v) / (if (coarse) { 40 } else { 400 })
      s <- .aaa_nudge(state(), panel, dx, cands_of(panel), snap(panel))
      if (!identical(s$ties, rv$ties)) {
        set_state(s)
      }
    }
    shiny::observeEvent(input$p_cl, nudge("proxy", -1, TRUE))
    shiny::observeEvent(input$p_fl, nudge("proxy", -1, FALSE))
    shiny::observeEvent(input$p_fr, nudge("proxy", 1, FALSE))
    shiny::observeEvent(input$p_cr, nudge("proxy", 1, TRUE))
    shiny::observeEvent(input$a_cl, nudge("astro", -1, TRUE))
    shiny::observeEvent(input$a_fl, nudge("astro", -1, FALSE))
    shiny::observeEvent(input$a_fr, nudge("astro", 1, FALSE))
    shiny::observeEvent(input$a_cr, nudge("astro", 1, TRUE))
    # typed positions of the selected tie point (the astronomical one snaps to the nearest peak)
    nearest_peak <- function(x) {
      cd <- rec_p()$cands
      if (nrow(cd) == 0) {
        NA_real_
      } else {
        cd$x[which.min(abs(cd$x - x))]
      }
    }
    # the offset is measured from the peak that was nearest when the tie point was selected or moved otherwise
    anchor <- NA_real_
    keep_anchor <- FALSE
    shiny::observe({
      k <- rv$sel
      if (is.null(k) || k > nrow(rv$ties)) {
        for (id in c("sel_p", "sel_a", "sel_off")) {
          shiny::updateNumericInput(session, id, value = NA)
        }
      } else {
        if (!keep_anchor || !is.finite(anchor)) {
          anchor <<- nearest_peak(rv$ties$p[k])
        }
        keep_anchor <<- FALSE
        shiny::updateNumericInput(session, "sel_p", value = signif(rv$ties$p[k], 8))
        shiny::updateNumericInput(session, "sel_a", value = signif(rv$ties$a[k], 8))
        shiny::updateNumericInput(session, "sel_off", value = signif(rv$ties$p[k] - anchor, 5))
      }
    })
    typed <- function(col, v, offset = FALSE) {
      k <- shiny::isolate(rv$sel)
      if (is.null(k) || is.null(v) || !is.finite(v)) {
        return()
      }
      if (offset) {
        if (!is.finite(anchor)) {
          anchor <<- nearest_peak(rv$ties$p[k])
        }
        v <- anchor + v
        keep_anchor <<- TRUE
      }
      if (col == "a") {
        v <- pka()$x[which.min(abs(pka()$x - v))]
      }
      if (abs(v - rv$ties[[col]][k]) > 1e-7 * max(1, abs(v))) {
        s <- state()
        s$ties[k, col] <- v
        set_state(s)
      }
    }
    shiny::observeEvent(input$sel_p, typed("p", input$sel_p), ignoreInit = TRUE)
    shiny::observeEvent(input$sel_a, typed("a", input$sel_a), ignoreInit = TRUE)
    shiny::observeEvent(input$sel_off, typed("p", input$sel_off, TRUE), ignoreInit = TRUE)

    rp <- range(proxy$x)
    shiny::observeEvent(input$fit, {
      rv$view_p <- rp
      rv$view_a <- ra()
    })
    shiny::observeEvent(input$sync, {   # move the astronomical window to the partners of the tie points in the proxy window
      if (nrow(rv$ties) == 0) {
        return()
      }
      in_p <- rv$ties$p >= rv$view_p[1] & rv$ties$p <= rv$view_p[2]
      in_a <- rv$ties$a >= rv$view_a[1] & rv$ties$a <= rv$view_a[2]
      if (any(in_p)) {
        w <- diff(rv$view_a)
        rv$view_a <- .aai_fix_view(mean(range(rv$ties$a[in_p])) + c(-0.5, 0.5) * w, ra(), w)
      } else if (any(in_a)) {
        w <- diff(rv$view_p)
        rv$view_p <- .aai_fix_view(mean(range(rv$ties$p[in_a])) + c(-0.5, 0.5) * w, rp, w)
      }
    })

    flag <- shiny::reactive(.aaa_flag(rv$ties, td()))
    tie_d <- shiny::debounce(shiny::reactive(rv$ties), 400)
    tuned <- shiny::reactive({
      input$tune
      ties <- if (isTRUE(input$live)) { tie_d() } else { shiny::isolate(rv$ties) }
      o <- .aaa_output(ties, proxy, astro_r())
      tryCatch(.aaa_tune(proxy, astro_r(), pka(), o, td(), input$flow, input$fhigh, input$roll, sedrate_factor),
          error = function(e) list(error = conditionMessage(e)))
    })
    output$plot_main <- shiny::renderPlot({
      b <- .aaa_boxes()
      s <- state()
      fl <- flag()
      rp_ <- rec_p()
      lab <- if (!is.null(cyc()) && is.null(cyc()$error)) { "proxy cycle" } else { "proxy" }
      .aaa_draw_record(rp_$dat, rp_$cands, rv$view_p, s, "proxy", fl, flip_p(), lab, b$proxy, FALSE, snap())
      .aaa_draw_record(astro_r(), pka(), rv$view_a, s, "astro", fl, flip_a(), "astronomical solution", b$astro, TRUE, TRUE)
      .aaa_draw_tuned(tuned(), astro_r(), s, pka(), isTRUE(input$show_astro), b$tuned)
      .aaa_draw_links(s, rv$view_p, rv$view_a, flip_p(), fl, flip_a())   # last: its coordinates are the click coordinates
    })
    output$sed_ui <- shiny::renderUI({
      if (!isTRUE(input$show_sed)) {
        return(NULL)
      }
      shiny::tagList(shiny::div(class = "aaa-label", "sedimentation rate between the tie points"),
          shiny::div(class = "aaa-plot", shiny::plotOutput("plot_s", height = "150px")))
    })
    output$plot_s <- shiny::renderPlot(.aaa_draw_sedrate(tuned()))

    output$status <- shiny::renderUI({
      n <- nrow(rv$ties)
      bad <- sum(flag())
      txt <- sprintf("%d tie point%s", n, if (n == 1) { "" } else { "s" })
      if (bad > 0) {
        txt <- paste0(txt, sprintf(", %d in the wrong order (crossed)", bad))
      }
      if (!is.null(cyc()$error)) {
        txt <- paste0(txt, " - filter: ", cyc()$error)
      }
      if (!is.null(rv$pending)) {
        txt <- paste0(txt, sprintf(" - now click the %s",
            if (rv$pending$panel == "proxy") { "astronomical solution" } else { "proxy record" }))
      }
      if (!is.null(rv$sel)) {
        txt <- paste0(txt, sprintf(" - tie %d selected: move it or click a free position",
            rank(rv$ties$p, ties.method = "first")[rv$sel]))
      }
      shiny::div(class = "aaa-status", txt)
    })
    output$ties_table <- shiny::renderTable({
      o <- out()
      if (nrow(o) == 0) {
        return(data.frame(`#` = integer(0), proxy = numeric(0), time = numeric(0), check.names = FALSE))
      }
      data.frame(`#` = seq_len(nrow(o)), proxy = signif(o$data_x, 6), time = signif(o$insolation_x, 6), check.names = FALSE)
    }, digits = 6, spacing = "xs", class = "aaa-ties")

    shiny::observeEvent(input$done, shiny::stopApp(out()))
    session$onSessionEnded(function() shiny::stopApp(shiny::isolate(out())))
  }
  shiny::shinyApp(ui, server)
}
