#' @title Plot Blocked / Interleaved eCOCO maps
#'
#' @description Plots the five eCOCO maps used by Acycle for Blocked and
#' Interleaved eCOCO: consensus correlation, consensus SR-global p and
#' same-rate local p (both as -log10(p) with their own 1-2-5 display range,
#' limited by the Monte Carlo resolution 1/(n_sim+1), with the 0.05 and 0.01
#' contours), the number of contributing orbital periods, and the
#' window-normalised ridge score with the tracked ridge (filled squares:
#' global p <= 0.05; filled circles: local p <= 0.01; open circles: neither).
#' Also accepts cocoR \code{ecoco_adaptive} results. Translated from Acycle
#' v3.1 ecocoplot.m (Mingsong Li); base graphics only.
#'
#' @param ecoco eCOCO result to plot: the list returned by \code{\link{ecoco_blocked}},
#' \code{\link{ecoco_interleaved}} or \code{\link{ecoco_adaptive}} (with \code{sr}, \code{depth}, \code{rho},
#' \code{p_global}, \code{p_local}, \code{norbit}, \code{score} and \code{track}). Required (\code{Default=NULL} gives
#' an error).
#' @param panels Which maps to draw, side by side in the order given (integer vector with values 1 to 5):
#' 1 consensus correlation coefficient (rho), 2 SR-global p, 3 same-rate local p, 4 number of contributing orbital
#' periods, 5 window-normalised ridge score with the tracked ridge. \code{Default=1:5}.
#' @param reverse_depth Reverse the depth axis so that depth increases downwards (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param palette Colours of the colour scale of all maps (character vector of colours, e.g. the output of
#' \code{\link[grDevices]{hcl.colors}}); the first colour is used for the lowest values and the length of the vector
#' sets the number of colour levels. \code{Default=} 64 viridis colours (dark purple to yellow).
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' @references Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking
#' variable sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' dt <- 0.05
#' depth <- seq(0, 30, by = dt)
#' tt <- depth * 100 / 5
#' x <- sin(2 * pi * tt / 405.7) + 0.8 * sin(2 * pi * tt / 98.9) +
#'   0.5 * sin(2 * pi * tt / 40.99) + 0.4 * sin(2 * pi * tt / 23.68) +
#'   0.3 * rnorm(length(tt))
#' res <- ecoco_blocked(cbind(depth, x), window = 10, sr = seq(3, 8, 0.5),
#'                      step_depth = 1, pad = 400, n_sim = 20)
#' plot_ecoco_cv(res)
#' }
#'
#' @return Invisibly, the list of plotted matrices (rho, -log10 p maps,
#' orbit count, normalised score) and the ridge.
#' @export
#' @importFrom graphics par image contour points lines layout axis box
#' @importFrom graphics mtext plot.new plot.window rect
#' @importFrom grDevices hcl.colors
plot_ecoco_cv <- function(ecoco = NULL, panels = 1:5, reverse_depth = FALSE,
    palette = .wr_hcl_cols(64, "viridis")){
  # Internal helper functions (local to this function)

  ecocox_plog <- function(p, pfloor){
    if (is.null(p)){
      return(NULL)
    }
    out <- p * NaN
    ok <- is.finite(p) & p > 0 & p <= 1
    v <- p[ok]
    if (is.finite(pfloor)){
      v <- pmax(v, pfloor)
    }
    out[ok] <- -log10(v)
    if (is.finite(pfloor)){
      out[is.finite(p) & p == 0] <- -log10(pfloor)
    }
    out
  }

  ecocox_zrange <- function(z, kind, maps, pfloor = NaN){
    if (kind == 1){
      v <- c(maps$rho[is.finite(maps$rho)])
      if (length(v) < 1){
        return(c(0, 1))
      }
      r <- range(v)
      if (r[1] == r[2]){
        r <- r + c(-1, 1) * max(1, abs(r[1])) * 1e-6
      }
      return(r)
    }
    if (kind %in% c(2, 3)){
      v <- z[is.finite(z)]
      if (length(v) < 1){
        return(c(0, 1))
      }
      pmin_v <- 10^(-max(v))
      ex <- floor(log10(pmin_v))
      man <- pmin_v / 10^ex
      lead <- 1
      if (man >= 5){
        lead <- 5
      } else if (man >= 2){
        lead <- 2
      }
      dm <- lead * 10^ex
      if (is.finite(pfloor)){
        dm <- max(dm, pfloor)
      }
      if (dm >= 1){
        dm <- 0.1
      }
      return(c(0, -log10(dm)))
    }
    if (kind == 4){
      v <- z[is.finite(z)]
      mx <- 9
      if (length(v) > 0){
        mx <- max(1, ceiling(max(v)))
      }
      return(c(0, mx))
    }
    c(0, 1)
  }

  ecocox_ridge_values <- function(M, sr, rate){
    out <- rep(NaN, length(rate))
    if (is.null(M)){
      return(out)
    }
    for (j in seq_len(min(length(rate), ncol(M)))){
      if (!is.finite(rate[j])){
        next
      }
      i <- which.min(abs(sr - rate[j]))
      if (is.finite(M[i, j])){
        out[j] <- M[i, j]
      }
    }
    out
  }

  # Main code
  if (is.null(ecoco)){
    stop("ecoco result required")
  }
  sr <- ecoco$sr
  depth <- ecoco$depth
  nr <- length(sr)
  nw <- length(depth)
  score <- ecoco$score
  ns <- matrix(NaN, nr, nw)
  for (j in seq_len(nw)){
    col <- score[, j]
    ok <- is.finite(col)
    if (!any(ok)){
      next
    }
    lo <- min(col[ok])
    hi <- max(col[ok])
    if (hi > lo){
      ns[ok, j] <- (col[ok] - lo) / (hi - lo)
    } else {
      ns[ok, j] <- 0
    }
  }
  pfloor <- ecoco$p_floor
  if (is.null(pfloor) || !is.finite(pfloor)){
    pfloor <- NaN
    if (!is.null(ecoco$nsim) && ecoco$nsim > 0){
      pfloor <- 1 / (ecoco$nsim + 1)
    }
  }
  maps <- list(rho = ecoco$rho,
      p_global = ecocox_plog(ecoco$p_global, pfloor),
      p_local = ecocox_plog(ecoco$p_local, pfloor),
      norbit = ecoco$norbit, score = ns)
  titles <- c("Consensus correlation", "Global p", "Local p",
      "Contributing orbital periods", "Ridge score")
  if (is.null(ecoco$forward)){
    titles[1] <- "Correlation coefficient"
  }
  ridge <- ecoco$track
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  np <- length(panels)
  graphics::layout(matrix(c(seq_len(np), np + seq_len(np)), 2, np,
      byrow = TRUE), heights = c(5, 1))
  graphics::par(mar = c(4, 4, 2.5, 1))
  ylim <- range(depth)
  if (reverse_depth){
    ylim <- rev(ylim)
  }
  zl <- vector("list", np)
  for (k in seq_len(np)){
    pi_k <- panels[k]
    z <- maps[[pi_k]]
    zr <- ecocox_zrange(z, pi_k, maps, pfloor)
    zl[[k]] <- zr
    zz <- z
    if (all(!is.finite(zz))){
      zz <- matrix(0, nr, nw)
    }
    zz[is.finite(zz) & zz > zr[2]] <- zr[2]
    zz[is.finite(zz) & zz < zr[1]] <- zr[1]
    ylab <- ""
    if (k == 1){
      ylab <- "Depth (m)"
    }
    graphics::image(sr, depth, zz, col = palette, zlim = zr, ylim = ylim,
        xlab = "Sedimentation rate (cm/kyr)", ylab = ylab,
        main = titles[pi_k], useRaster = FALSE)
    if (pi_k %in% c(2, 3) && nr > 1 && nw > 1){
      lev <- -log10(0.01)
      if (pi_k == 2){
        lev <- -log10(0.05)
      }
      zc <- z
      zc[!is.finite(zc)] <- 0
      if (any(zc > lev) && any(zc < lev)){
        graphics::contour(sr, depth, zc, levels = lev, add = TRUE,
            drawlabels = FALSE, lwd = 0.6)
      }
    }
    if (pi_k == 5 && !is.null(ridge)){
      rr <- ridge$sr
      dd <- ridge$depth
      ok <- is.finite(rr) & is.finite(dd)
      if (any(ok)){
        pg <- ecocox_ridge_values(ecoco$p_global, sr, rr)
        pl <- ecocox_ridge_values(ecoco$p_local, sr, rr)
        gs <- ok & is.finite(pg) & pg <= 0.05
        ls <- ok & !gs & is.finite(pl) & pl <= 0.01
        nsg <- ok & !gs & !ls
        graphics::lines(rr, dd, col = "red", lwd = 1)
        graphics::points(rr[nsg], dd[nsg], pch = 1, col = "red", cex = 0.6)
        graphics::points(rr[ls], dd[ls], pch = 16, col = "red", cex = 0.7)
        graphics::points(rr[gs], dd[gs], pch = 15, col = "red", cex = 0.7)
      }
    }
  }
  graphics::par(mar = c(3, 4, 0.5, 1))
  for (k in seq_len(np)){
    pi_k <- panels[k]
    zr <- zl[[k]]
    graphics::plot.new()
    graphics::plot.window(xlim = zr, ylim = c(0, 1), xaxs = "i")
    brk <- seq(zr[1], zr[2], length.out = length(palette) + 1)
    graphics::rect(brk[-length(brk)], 0, brk[-1], 1, col = palette,
        border = NA)
    graphics::box()
    if (pi_k %in% c(2, 3)){
      at <- pretty(zr)
      at <- at[at >= zr[1] & at <= zr[2]]
      graphics::axis(1, at = at, labels = signif(10^(-at), 3))
      graphics::mtext("p", side = 1, line = 2, cex = 0.7)
    } else {
      graphics::axis(1)
      lab <- c("rho", "", "", "Number of periods", "Window-normalized score")
      graphics::mtext(lab[pi_k], side = 1, line = 2, cex = 0.7)
    }
  }
  invisible(list(maps = maps, ridge = ridge))
}
