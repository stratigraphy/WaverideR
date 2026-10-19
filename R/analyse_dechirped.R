#' De-chirped scalogram of any transform: the record warped to time with a track, transformed, mapped back to depth
#'
#' The universal form of the de-chirping of \code{\link{analyse_chirp_eha}}: the record is converted to the time of a
#' tracked sedimentation rate (any \code{wr_track}: \code{\link{track_chirp_hmm}}, \code{\link{track_astro_hmm}},
#' \code{\link{track_astro_auto}}), so that all cycles become stationary where the track is right, the chosen
#' transform (Morlet wavelet, superlet-SST, multisynchrosqueezing, chirplet-SST or multitaper EHA) is computed in that
#' time domain (the "tuned" scalogram, with the periods in kyr), and the result is mapped back onto the depth axis
#' with the period axis in depth units at every position (period in depth = period in kyr x local rate), so that
#' the cycles are shown as sharp lines that follow the changes of the rate. Where the track is wrong the lines do not
#' straighten in time and do not sharpen in depth, which makes the map a check of the track.
#'
#' With \code{n_sim > 0} the same warp and transform are applied to AR1 surrogates of the record and the power is
#' compared with the 95th percentile of the surrogates cell by cell (the warp changes the noise as well, so wavy
#' structure in the de-chirped map is only evidence where it exceeds the warped noise).
#'
#' @param data Record: first column depth, second column proxy.
#' @param track A \code{wr_track} (or \code{\link{track_astro_auto}} result) with \code{depth}, \code{sedrate} and
#' \code{time}.
#' @param transform \code{"wavelet"}, \code{"superlet_sst"}, \code{"wavelet_sst"} (multisynchrosqueezing),
#' \code{"chirplet_sst"}, \code{"eha"} or \code{"combined"}: the geometric mean of the rank-normalised Morlet,
#' multitaper EHA and superlet-SST maps (the Morlet is best at the long periods, the EHA at the short ones and the
#' superlet-SST localises the periods best; the product keeps what all three agree on).
#' @param lowerPeriod,upperPeriod Period range in kyr \code{Default=10 and 600}.
#' @param dj Period spacing \code{Default=1/50}.
#' @param omega_nr Morlet omega_0 for \code{"wavelet"} (and the orders of the superlets / chirplets when a vector)
#' \code{Default=10}.
#' @param window Window length in kyr for \code{"eha"} \code{Default=20 x lowerPeriod}.
#' @param n_sim Number of AR1 surrogates for the significance \code{Default=0}.
#' @param seed Seed of the surrogates.
#' @param ... Further arguments of the transform.
#' @return A list of class \code{wr_dechirped} with \code{tuned} (the transform in time: an \code{analyze.wavelet}
#' object with \code{Period} in kyr, usable in \code{\link{plot_wavelet}}), \code{depth} (an \code{analyze.wavelet}
#' object on the depth axis with \code{Period} in depth units: the de-chirped scalogram), \code{Power_norm} (the
#' depth map with every column divided by its smooth continuum), and with \code{n_sim > 0} \code{Power_sig} (ratio
#' to the 95th percentile of the warped surrogates) and \code{Power_pval}.
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 2400, by = 1)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.7 * sin(2 * pi * t / 600)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1)
#' res <- track_astro_auto(data.frame(z, x), sedrate = c(0.8, 5), astro_periods = c(41, 23.7, 19))
#' dc <- analyse_dechirped(data.frame(z, x), res$track, "superlet_sst", lowerPeriod = 10,
#'     upperPeriod = 120)
#' plot_dechirped(dc)
#' }
#' @export
analyse_dechirped <- function(data, track, transform = c("wavelet", "superlet_sst", "wavelet_sst", "chirplet_sst", "eha", "combined"),
                              lowerPeriod = 10, upperPeriod = 600, dj = 1 / 50, omega_nr = 10, window = NULL,
                              n_sim = 0, seed = 1, ...) {
  transform <- match.arg(transform)
  if (inherits(track, "wr_track_auto")) track <- track$track
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dz <- stats::median(diff(dat[, 1]))
  z <- seq(dat[1, 1], dat[nrow(dat), 1], by = dz)
  x <- stats::approx(dat[, 1], dat[, 2], z)$y
  # time of every depth from the track (extrapolated with the end rates beyond its range)
  tz <- stats::approx(track$depth, track$time, z, rule = 2)$y
  rz <- stats::approx(track$depth, track$sedrate, z, rule = 2)$y
  below <- z < min(track$depth); above <- z > max(track$depth)
  if (any(below)) tz[below] <- track$time[1] - (min(track$depth) - z[below]) / track$sedrate[1]
  if (any(above)) tz[above] <- utils::tail(track$time, 1) + (z[above] - max(track$depth)) / utils::tail(track$sedrate, 1)
  if (any(diff(tz) <= 0)) tz <- cummax(tz + seq_along(tz) * 1e-9)
  n <- length(z)
  dtau <- diff(range(tz)) / (n - 1)
  tau <- seq(min(tz), max(tz), by = dtau)
  run <- function(xx) {
    xt <- stats::approx(tz, xx, tau)$y
    d <- data.frame(tau, xt)
    switch(transform,
      wavelet = analyse_wavelet(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], verbose = FALSE, ...),
      superlet_sst = analyse_superlet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj,
                                          omega_nr = if (length(omega_nr) > 1) omega_nr else c(8, 12, 16), ...),
      wavelet_sst = analyse_wavelet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], n_iter = 4, ...),
      chirplet_sst = analyse_chirplet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj,
                                          omega_nr = if (length(omega_nr) > 1) omega_nr else c(8, 12, 16), ...),
      eha = {
        w <- if (is.null(window)) 20 * lowerPeriod else window
        e <- analyse_chirp_eha(d, window = w, chirp = 0, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, ...)
        e$Power <- e$Power_norm; e
      },
      combined = {
        # consensus of the three: Morlet (best at the long periods), multitaper EHA (best at the short ones) and
        # superlet-SST (best period localisation), each rank-normalised, combined by their geometric mean on the
        # grid of the wavelet (the superlet principle applied across methods)
        wv <- analyse_wavelet(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], verbose = FALSE)
        su <- analyse_superlet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj)
        w <- if (is.null(window)) 20 * lowerPeriod else window
        e <- analyse_chirp_eha(d, window = w, chirp = 0, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj)
        ongrid <- function(P, ax1, ax2) {                 # interpolate a map onto the wavelet grid
          P1 <- t(apply(P, 1, function(v) stats::approx(ax1, v, wv$axis.1, rule = 2)$y))
          apply(P1, 2, function(v) stats::approx(ax2, v, wv$axis.2, rule = 2)$y)
        }
        rk <- function(P) { P[!is.finite(P)] <- 0; matrix(rank(P) / length(P), nrow(P), ncol(P)) }
        Pw <- rk(wv$Power / rowMeans(wv$Power))
        Pe <- rk(ongrid(e$Power_norm, e$axis.1, e$axis.2))
        Ps <- rk(ongrid(su$Power, su$axis.1, su$axis.2))
        wv$Power <- (Pw * Pe * Ps)^(1 / 3)
        wv$Power_components <- list(wavelet = Pw, eha = Pe, superlet_sst = Ps)
        wv
      })
  }
  tuned <- run(x)
  # map back to depth: columns at tau(z), period axis shifted by log2 of the local rate
  lp_t <- tuned$axis.2; ax_t <- tuned$axis.1
  lp_z <- seq(min(lp_t) + log2(min(rz)), max(lp_t) + log2(max(rz)), by = dj)
  to_depth <- function(P) {
    Pz <- t(apply(P, 1, function(v) stats::approx(ax_t, v, tz, rule = 2)$y))    # periods(kyr) x depth
    out <- matrix(NA_real_, length(lp_z), n)
    for (j in seq_len(n)) out[, j] <- stats::approx(lp_t + log2(rz[j]), Pz[, j], lp_z, rule = 1)$y
    out
  }
  Pd <- to_depth(tuned$Power)
  kb <- max(3, 2 * floor(0.5 / dj) + 1)
  Pn <- apply(Pd, 2, function(v) {
    ok <- is.finite(v) & v > 0
    if (sum(ok) < 10) return(v * NA)
    bg <- rep(NA_real_, length(v)); bg[ok] <- exp(stats::lowess(lp_z[ok], log(v[ok]), f = 0.4, iter = 2)$y)
    v / bg
  })
  depth_obj <- tuned
  depth_obj$Power <- Pd; depth_obj$Period <- 2^lp_z; depth_obj$axis.2 <- lp_z; depth_obj$axis.1 <- z
  depth_obj$Wave <- sqrt(pmax(Pd, 0)) + 0i; depth_obj$Ampl <- sqrt(pmax(Pd, 0)); depth_obj$Phase <- NULL
  depth_obj$Scale <- depth_obj$Period / depth_obj$fourier_factor; depth_obj$nr <- length(lp_z); depth_obj$nc <- n
  depth_obj$coi.1 <- NULL; depth_obj$coi.2 <- NULL; depth_obj$x <- z; depth_obj$y <- x
  depth_obj$Power.avg <- rowMeans(Pd, na.rm = TRUE); depth_obj$dt <- dz
  out <- list(tuned = tuned, depth = depth_obj, Power_norm = Pn, transform = transform, track = track, rate = rz, time = tz)
  if (n_sim > 0) {
    set.seed(seed)
    xc <- x - mean(x); rho <- stats::cor(xc[-1], xc[-n])
    sims <- array(NA_real_, c(dim(Pn), n_sim))
    sparse <- transform %in% c("superlet_sst", "wavelet_sst", "chirplet_sst")
    for (i in seq_len(n_sim)) {
      e <- as.numeric(stats::arima.sim(list(ar = rho), n = n, sd = sqrt(1 - rho^2))); e <- e / stats::sd(e) * stats::sd(xc)
      ts <- run(e)
      Ps <- to_depth(ts$Power)
      Pns <- if (sparse) Ps else apply(Ps, 2, function(v) { ok <- is.finite(v) & v > 0; if (sum(ok) < 10) return(v * NA)
        bg <- rep(NA_real_, length(v)); bg[ok] <- exp(stats::lowess(lp_z[ok], log(v[ok]), f = 0.4, iter = 2)$y); v / bg })
      sims[, , i] <- Pns
    }
    # the AR1 surrogates are stationary in depth, so for every period (row) the surrogate values of all positions and
    # all surrogates form one null distribution: far more stable than per cell with a few surrogates. Squeezed maps
    # are mostly empty, so there a ridge is compared with the ridges noise produces (the positive cells).
    obs <- if (sparse) Pd else Pn
    q95 <- rep(NA_real_, nrow(obs)); pv <- matrix(NA_real_, nrow(obs), ncol(obs))
    for (r in seq_len(nrow(obs))) {
      v <- as.numeric(sims[r, , ]); v <- v[is.finite(v)]
      if (sparse) v <- v[v > 0]
      if (length(v) < 20) next
      q95[r] <- stats::quantile(v, 0.95, names = FALSE)
      Fr <- stats::ecdf(v)
      pv[r, ] <- 1 - Fr(obs[r, ]) + 1 / length(v)
    }
    pv[!is.finite(obs)] <- NA
    if (sparse) pv[is.finite(obs) & obs <= 0] <- 1
    out$Power_sig <- obs / pmax(q95, 1e-12 * max(obs, na.rm = TRUE))
    out$Power_pval <- pmin(pv, 1)
  }
  class(out) <- "wr_dechirped"
  out
}

#' Plot a de-chirped scalogram
#'
#' @param x Result of \code{\link{analyse_dechirped}}.
#' @param style \code{"outline"} (de-chirped map with the significant cells outlined, when surrogates were run),
#' \code{"significance"} (ratio to the 95th percentile of the warped surrogates), \code{"plain"} (the de-chirped map)
#' or \code{"tuned"} (the transform in time, periods in kyr).
#' @param show Astronomical periods (kyr) to draw at the tracked rate.
#' @param level Significance level of the outline \code{Default=0.95}.
#' @param palette \code{hcl.colors} palette \code{Default="viridis"}.
#' @param main Title.
#' @param xlab,ylab Axis labels.
#' @param scale Colour scale as in all image plots of the package (\code{"auto"}, \code{"log"}, \code{"quantile"},
#' \code{"sqrt"}, \code{"linear"}; see \code{\link{plot_wavelet}}) \code{Default="auto"}. The significance style
#' always uses a fixed diverging scale centred on the 95th percentile of the surrogates.
#' @param outline_lwd Line width of the significance outline \code{Default=0.6}.
#' @return Invisibly, the matrix shown.
#' @export
plot_dechirped <- function(x, style = c("outline", "significance", "plain", "tuned"), show = NULL, level = 0.95,
                           palette = "viridis", main = NULL, xlab = "Depth", ylab = "Period", scale = "auto",
                           outline_lwd = 0.6) {
  style <- match.arg(style)
  sparse <- x$transform %in% c("superlet_sst", "wavelet_sst", "chirplet_sst")
  if (style %in% c("outline", "significance") && is.null(x$Power_sig)) {
    message("plot_dechirped: no surrogates were run (n_sim = 0 in analyse_dechirped), drawing style = 'plain'")
    style <- "plain"
  }
  if (style == "tuned") {
    M <- x$tuned$Power; if (!sparse) M <- M / rowMeans(M, na.rm = TRUE)     # flatten the red-noise slope
    ax1 <- x$tuned$axis.1; ax2 <- x$tuned$axis.2
    if (is.null(main)) main <- paste(x$transform, "of the record in the time of the track (tuned): cycles should be horizontal lines")
    xl <- "Time of the track (kyr)"; yl <- "Period (kyr)"
  } else {
    M <- if (style == "significance") log10(pmax(x$Power_sig, 1e-6)) else if (sparse) x$depth$Power else x$Power_norm
    ax1 <- x$depth$axis.1; ax2 <- x$depth$axis.2
    if (is.null(main)) main <- switch(style,
      significance = paste(x$transform, "de-chirped by the track, ratio to the 95th percentile of warped surrogates (log10; > 0 significant)"),
      outline = paste(x$transform, "de-chirped by the track, cells significant against warped surrogates outlined"),
      plain = paste(x$transform, "de-chirped by the track (lines against the continuum of each position)"))
    xl <- xlab; yl <- ylab
  }
  # empty space outside the mapped period range is drawn in the lowest colour, not left white
  M <- .wr_fill_na(M)
  if (style == "significance") {
    # diverging scale with the threshold (the 95th percentile of the surrogates, log10 ratio 0) in the middle:
    # blue below (down to a tenth of the threshold), red above
    .wr_image_signif(ax1, ax2, M, yaxt = "n", xlab = xl, ylab = yl, main = main, xaxs = "i")
  } else {
    .wr_image(ax1, ax2, M, scale = scale, palette = palette, sparse = if (sparse) TRUE else NULL, fill_na = FALSE,
              yaxt = "n", xlab = xl, ylab = yl, main = main, xaxs = "i")
  }
  ticks <- 2^(floor(min(ax2)):ceiling(max(ax2)))
  graphics::axis(2, at = log2(ticks), labels = ticks)
  if (style == "outline") {
    sig <- (x$Power_pval <= 1 - level) * 1; sig[is.na(sig)] <- 0; sig[is.na(x$Power_pval)] <- 0
    graphics::contour(ax1, ax2, t(sig), levels = 0.5, drawlabels = FALSE, add = TRUE, col = "white", lwd = outline_lwd)
  }
  if (!is.null(show)) {
    if (style == "tuned") graphics::abline(h = log2(show), col = "white", lty = 3, lwd = 0.9)
    else for (p in show) graphics::lines(x$depth$axis.1, log2(p * x$rate), col = if (style == "significance") "black" else "white", lwd = 0.9, lty = 3)
  }
  invisible(M)
}
