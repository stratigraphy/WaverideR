#' @title Track the sedimentation rate through an eCOCO map
#'
#' @description
#' Interactively select points in the map of an eCOCO result (sedimentation rate
#' against depth) to trace how the sedimentation rate changes along the record. This is the eCOCO
#' counterpart of \code{\link{track_period}}: the map of one of the eCOCO statistics (by default the combined
#' correlation and significance score \code{pcoco}) is plotted, the local maxima (ridges) are marked, and the
#' sedimentation rate is tracked by clicking ridge points.
#'
#' Selected points are shown as white dots. Clicking a selected point again deselects it (red dot). Press the
#' right mouse button or Esc (depending on the graphics device) to finish. An interactive device
#' that supports \code{\link[graphics]{identify}} is needed (on Windows this works out of the box, on macOS
#' XQuartz or the Quartz device, on Linux an X11 device).
#'
#' The tracking can also be done without a mouse: pass the points as \code{points}. They are snapped to the nearest
#' ridge point, which makes a tracking reproducible and scriptable.
#'
#' @param ecoco An eCOCO result from \code{\link{ecoco_adaptive}}, \code{\link{ecoco_blocked}} or
#' \code{\link{ecoco_interleaved}} (a list with \code{sr}, \code{depth} and rate-by-window matrices).
#' @param map The map used for tracking and plotting: \code{"pcoco"}, \code{"rho"}, \code{"score"},
#' \code{"p_global"} (shown as minus log10 p) or \code{"norbit"}. \code{Default="pcoco"}.
#' @param track_peaks If \code{TRUE} only the local maxima along the sedimentation rate axis can be selected,
#' otherwise every grid cell \code{Default=TRUE}.
#' @param n_peaks Maximum number of peaks per window that can be selected (the highest ones) \code{Default=3}.
#' @param rho_min,p_max,norbit_min Optional criteria a peak must meet: minimum correlation, maximum p-value
#' (global) and minimum number of contributing orbital periods. \code{Default=NULL} (not used).
#' @param points Optional data frame or matrix with columns depth and sedimentation rate. If given, no mouse input
#' is needed: the points are snapped to the nearest selectable point.
#' @param interpolate Linearly interpolate the tracked sedimentation rate to all windows of the eCOCO result
#' (constant outside the tracked interval) \code{Default=FALSE}.
#' @param n.levels Number of colour levels \code{Default=100}.
#' @param plot_horizontal Plot depth on the x-axis (\code{TRUE}) or on the y-axis (\code{FALSE}) \code{Default=TRUE}.
#' @param plot_dir Logical. \code{TRUE} if depth increases along the axis, \code{FALSE} to reverse the depth
#' axis (for example for borehole data with depth increasing downwards) \code{Default=TRUE}.
#' @param x_lab Label of the depth axis \code{Default="depth (metres)"}.
#' @param y_lab Label of the sedimentation rate axis \code{Default="sedimentation rate (cm/kyr)"}.
#' @param palette_name Name of the \code{\link[grDevices]{hcl.colors}} palette \code{Default="viridis"}.
#' @param genplot Plot the result (only needed when \code{points} is used; the interactive mode always plots)
#' \code{Default=FALSE}.
#'
#' @return A data frame with the tracked windows: \code{depth}, \code{sedrate} (cm/kyr) and, when available in
#' \code{ecoco}, \code{rho}, \code{p_global}, \code{norbit} and the value of the map. The result can be used like the
#' output of \code{\link{track_period}} (for example to build an age model).
#'
#' @author Michiel Arts. Inspired by \code{\link{track_period}} and the tracking routines of Acycle
#' (Mingsong Li); see \code{\link{ecoco_track}} for the automatic version.
#'
#' @seealso \code{\link{ecoco_track}}, \code{\link{ecoco_track_ridge}}, \code{\link{track_period}}
#'
#' @examples
#' # synthetic eCOCO-like map: sedimentation rate increases with depth
#' sr <- seq(1, 10, by = 0.25)
#' depth <- seq(0, 20, by = 0.5)
#' true_sr <- 4 + 2 * depth / 20
#' rho <- sapply(true_sr, function(s) 0.8 * exp(-(sr - s)^2 / 0.5))
#' ecoco <- list(sr = sr, depth = depth, rho = rho, pcoco = 2 * rho,
#'               p_global = matrix(0.01, length(sr), length(depth)),
#'               norbit = matrix(9, length(sr), length(depth)))
#' # scripted tracking: points are snapped to the nearest ridge point
#' trk <- track_ecoco(ecoco, points = data.frame(depth = c(2, 10, 18), sr = c(4.2, 5, 5.8)))
#' trk
#' \donttest{
#' # interactive tracking (needs a graphics device with mouse input)
#' if (interactive()) trk <- track_ecoco(ecoco)
#' }
#' @export
#' @importFrom graphics image points identify axis box lines par layout mtext
#' @importFrom grDevices hcl.colors xy.coords
#' @importFrom stats approx
track_ecoco <- function(ecoco = NULL,
                        map = c("pcoco", "rho", "score", "p_global", "norbit"),
                        track_peaks = TRUE,
                        n_peaks = 3,
                        rho_min = NULL,
                        p_max = NULL,
                        norbit_min = NULL,
                        points = NULL,
                        interpolate = FALSE,
                        n.levels = 100,
                        plot_horizontal = TRUE,
                        plot_dir = TRUE,
                        x_lab = "depth (metres)",
                        y_lab = "sedimentation rate (cm/kyr)",
                        palette_name = "viridis",
                        genplot = FALSE) {
  map <- match.arg(map)
  if (is.null(ecoco) || is.null(ecoco$sr) || is.null(ecoco$depth)) {
    stop("ecoco must be an eCOCO result with the elements 'sr' and 'depth'")
  }
  sr <- as.numeric(ecoco$sr)
  depth <- as.numeric(ecoco$depth)
  Z <- ecoco[[map]]
  if (is.null(Z)) stop("the eCOCO result has no element '", map, "'")
  Z <- as.matrix(Z)
  if (!all(dim(Z) == c(length(sr), length(depth)))) {
    stop("'", map, "' must be a matrix with one row per sedimentation rate and one column per window")
  }
  Zplot <- if (map == "p_global") -log10(pmax(Z, .Machine$double.eps)) else Z
  Zplot[!is.finite(Zplot)] <- NA
  nsr <- length(sr)
  nw <- length(depth)

  # --- selectable points -----------------------------------------------------------------------------------------
  get_mat <- function(nm) {
    m <- ecoco[[nm]]
    if (is.null(m)) NULL else as.matrix(m)
  }
  ok <- matrix(TRUE, nsr, nw)
  if (!is.null(rho_min) && !is.null(get_mat("rho"))) ok <- ok & (get_mat("rho") >= rho_min)
  if (!is.null(p_max) && !is.null(get_mat("p_global"))) ok <- ok & (get_mat("p_global") <= p_max)
  if (!is.null(norbit_min) && !is.null(get_mat("norbit"))) ok <- ok & (get_mat("norbit") >= norbit_min)
  ok[is.na(ok)] <- FALSE
  ok <- ok & is.finite(Zplot)
  if (track_peaks) {
    up <- rbind(rep(-Inf, nw), Zplot[-nsr, , drop = FALSE])      # neighbour at lower rate
    dn <- rbind(Zplot[-1, , drop = FALSE], rep(-Inf, nw))        # neighbour at higher rate
    up[is.na(up)] <- -Inf
    dn[is.na(dn)] <- -Inf
    peak <- ok & (Zplot > up) & (Zplot >= dn)
    if (is.finite(n_peaks)) {                                      # keep the highest n_peaks per window
      for (j in seq_len(nw)) {
        idx <- which(peak[, j])
        if (length(idx) > n_peaks) {
          drop <- idx[order(Zplot[idx, j], decreasing = TRUE)[-seq_len(n_peaks)]]
          peak[drop, j] <- FALSE
        }
      }
    }
    sel_mat <- peak
  } else {
    sel_mat <- ok
  }
  cand <- which(sel_mat, arr.ind = TRUE)
  if (nrow(cand) == 0) stop("no selectable points: relax rho_min, p_max or norbit_min")
  cand_depth <- depth[cand[, 2]]
  cand_sr <- sr[cand[, 1]]

  # --- plot helpers -------------------------------------------------------------------------------------------------
  xy_of <- function(d, s) if (plot_horizontal) list(x = d, y = s) else list(x = s, y = d)
  draw_map <- function() {
    cols <- hcl.colors(n.levels, palette_name)
    zr <- range(Zplot, na.rm = TRUE)
    if (zr[1] == zr[2]) zr <- zr + c(-0.5, 0.5)
    brks <- seq(zr[1], zr[2], length.out = n.levels + 1)
    dlim <- range(depth)
    if (!plot_dir) dlim <- rev(dlim)
    slim <- range(sr)
    if (plot_horizontal) {
      image(x = depth, y = sr, z = t(Zplot), col = cols, breaks = brks, useRaster = TRUE,
            xlab = x_lab, ylab = y_lab, xlim = dlim, ylim = slim)
    } else {
      image(x = sr, y = depth, z = Zplot, col = cols, breaks = brks, useRaster = TRUE,
            xlab = y_lab, ylab = x_lab, ylim = dlim, xlim = slim)
    }
    box()
    pxy <- xy_of(cand_depth, cand_sr)
    graphics::points(pxy$x, pxy$y, pch = 1, col = "black", lwd = 0.5)
  }

  # --- selection ------------------------------------------------------------------------------------------------------
  n <- nrow(cand)
  chosen <- rep(FALSE, n)
  if (!is.null(points)) {
    pts <- as.matrix(points)
    if (ncol(pts) < 2) stop("points must have two columns: depth and sedimentation rate")
    rd <- diff(range(depth)); if (rd == 0) rd <- 1
    rs <- diff(range(sr)); if (rs == 0) rs <- 1
    for (k in seq_len(nrow(pts))) {
      dd <- ((cand_depth - pts[k, 1]) / rd)^2 + ((cand_sr - pts[k, 2]) / rs)^2
      chosen[which.min(dd)] <- TRUE
    }
    if (genplot) {
      draw_map()
      pxy <- xy_of(cand_depth[chosen], cand_sr[chosen])
      graphics::points(pxy$x, pxy$y, pch = 19, col = "white")
    }
  } else {
    wr_open_device(width = 12, height = 7, interactive_input = TRUE)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(4, 4, 2, 1))
    draw_map()
    pxy <- xy_of(cand_depth, cand_sr)
    defaultW <- getOption("warn")
    options(warn = -1)
    on.exit(options(warn = defaultW), add = TRUE)
    xy <- xy.coords(pxy$x, pxy$y)
    deselected <- rep(FALSE, n)
    while (sum(chosen) < n) {
      ans <- identify(xy$x, xy$y, n = 1, plot = FALSE, tolerance = 0.1)
      if (!length(ans)) break
      if (!chosen[ans]) {
        chosen[ans] <- TRUE
        deselected[ans] <- FALSE
      } else {
        chosen[ans] <- FALSE
        deselected[ans] <- TRUE
      }
      graphics::points(xy$x[chosen], xy$y[chosen], pch = 19, col = "white")
      graphics::points(xy$x[deselected], xy$y[deselected], pch = 19, col = "red")
    }
  }

  # --- result ---------------------------------------------------------------------------------------------------------
  idx <- cand[chosen, , drop = FALSE]
  if (nrow(idx) == 0) {
    return(data.frame(depth = numeric(0), sedrate = numeric(0)))
  }
  out <- data.frame(depth = depth[idx[, 2]], sedrate = sr[idx[, 1]])
  for (nm in c("rho", "p_global", "norbit")) {
    m <- get_mat(nm)
    if (!is.null(m)) out[[nm]] <- m[idx]
  }
  out[[map]] <- Z[idx]
  out <- out[order(out$depth), , drop = FALSE]
  # one value per window (mean when several peaks of the same window were selected)
  if (anyDuplicated(out$depth)) {
    agg <- stats::aggregate(out[, -1, drop = FALSE], by = list(depth = out$depth), FUN = mean)
    out <- agg
  }
  rownames(out) <- NULL
  if (interpolate && nrow(out) >= 1) {
    full <- data.frame(depth = depth)
    full$sedrate <- if (nrow(out) == 1) rep(out$sedrate, nw) else
      approx(out$depth, out$sedrate, xout = depth, rule = 2)$y
    out <- full
  }
  if (!is.null(points) && genplot) {
    pxy <- xy_of(out$depth, out$sedrate)
    lines(pxy$x, pxy$y, col = "white", lwd = 2)
  }
  out
}
