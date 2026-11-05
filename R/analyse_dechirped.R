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
#' @param data Input record in the depth domain: a matrix or data frame with depth in column 1 and the proxy value
#' in column 2 (further columns are ignored). Incomplete rows are removed, the rows are sorted and the record is
#' interpolated to the median depth spacing. Required, no default.
#' @param track Sedimentation-rate track used to warp the record to time: a \code{wr_track} (from
#' \code{\link{track_chirp_hmm}} or \code{\link{track_astro_hmm}}) or a \code{\link{track_astro_auto}} result
#' (whose \code{track} element is used), i.e. a list with \code{depth}, \code{sedrate} (depth units per kyr) and
#' \code{time} (kyr). Beyond the depth range of the track the time is extrapolated with the end rates. Required, no
#' default.
#' @param transform Transform computed in the time domain (character), one of:
#' \itemize{
#' \item \code{"wavelet"}: Morlet wavelet (\code{\link{analyse_wavelet}});
#' \item \code{"superlet_sst"}: superlet synchrosqueezing (\code{\link{analyse_superlet_sst}});
#' \item \code{"wavelet_sst"}: multisynchrosqueezing (\code{\link{analyse_wavelet_sst}} with \code{n_iter = 4});
#' \item \code{"chirplet_sst"}: chirplet synchrosqueezing (\code{\link{analyse_chirplet_sst}});
#' \item \code{"eha"}: multitaper EHA (\code{\link{analyse_chirp_eha}} with \code{chirp = 0});
#' \item \code{"combined"}: the geometric mean of the rank-normalised Morlet, multitaper EHA and superlet-SST maps
#' (the Morlet is best at the long periods, the EHA at the short ones and the superlet-SST localises the periods
#' best; the product keeps what all three agree on).
#' }
#' \code{Default="wavelet"} (the first element of the list of choices).
#' @param lowerPeriod,upperPeriod Shortest and longest period of the transform in the time domain, in kyr (numeric,
#' 0 < \code{lowerPeriod} < \code{upperPeriod}). Defaults: \code{lowerPeriod}: \code{Default=10};
#' \code{upperPeriod}: \code{Default=600}.
#' @param dj Spacing of the period axis in octaves (numeric, > 0), used both in time and for the period axis of
#' the depth map. \code{Default=1/50}.
#' @param omega_nr Number of cycles omega_0 of the Morlet wavelet for \code{"wavelet"}, \code{"wavelet_sst"} and
#' \code{"combined"} (numeric; the first value is used). For \code{"superlet_sst"} and \code{"chirplet_sst"} a
#' vector of several values gives the orders; a single value uses the orders \code{c(8, 12, 16)}.
#' \code{Default=10}.
#' @param window Length of the moving window of the \code{"eha"} and \code{"combined"} transforms, in kyr (numeric,
#' > 0); \code{NULL} uses \code{20 * lowerPeriod}. \code{Default=NULL}.
#' @param n_sim Number of AR(1) surrogates of the record that are warped and transformed in the same way to test
#' the significance of the de-chirped power (integer, >= 0; 0 skips the test). \code{Default=0}.
#' @param seed Seed of the random number generator for the surrogates (integer), so that the significance test can
#' be reproduced; only used when \code{n_sim > 0}. \code{Default=1}.
#' @param ... Further arguments passed to the transform function chosen with \code{transform}
#' (\code{\link{analyse_wavelet}}, \code{\link{analyse_superlet_sst}}, \code{\link{analyse_wavelet_sst}},
#' \code{\link{analyse_chirplet_sst}} or \code{\link{analyse_chirp_eha}}); not used for \code{"combined"}.
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
  if (inherits(track, "wr_track_auto")) {
    track <- track$track
  }
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dz <- stats::median(diff(dat[, 1]))
  z <- seq(dat[1, 1], dat[nrow(dat), 1], by = dz)
  x <- stats::approx(dat[, 1], dat[, 2], z)$y
  # time of every depth from the track (extrapolated with the end rates beyond its range)
  tz <- stats::approx(track$depth, track$time, z, rule = 2)$y
  rz <- stats::approx(track$depth, track$sedrate, z, rule = 2)$y
  below <- z < min(track$depth)
  above <- z > max(track$depth)
  if (any(below)) {
    tz[below] <- track$time[1] - (min(track$depth) - z[below]) / track$sedrate[1]
  }
  if (any(above)) {
    tz[above] <- utils::tail(track$time, 1) + (z[above] - max(track$depth)) / utils::tail(track$sedrate, 1)
  }
  if (any(diff(tz) <= 0)) {
    tz <- cummax(tz + seq_along(tz) * 1e-9)
  }
  n <- length(z)
  dtau <- diff(range(tz)) / (n - 1)
  tau <- seq(min(tz), max(tz), by = dtau)
  run <- function(xx) {
    xt <- stats::approx(tz, xx, tau)$y
    d <- data.frame(tau, xt)
    switch(transform,
        wavelet = analyse_wavelet(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], verbose = FALSE, ...),
        superlet_sst = analyse_superlet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj,
            omega_nr = if (length(omega_nr) > 1) { omega_nr } else { c(8, 12, 16) }, ...),
        wavelet_sst = analyse_wavelet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], n_iter = 4, ...),
        chirplet_sst = analyse_chirplet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj,
            omega_nr = if (length(omega_nr) > 1) { omega_nr } else { c(8, 12, 16) }, ...),
        eha = {
          w <- if (is.null(window)) { 20 * lowerPeriod } else { window }
          e <- analyse_chirp_eha(d, window = w, chirp = 0, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, ...)
          e$Power <- e$Power_norm
          e
        },
        combined = {
          # consensus of the three: Morlet (best at the long periods), multitaper EHA (best at the short ones) and
          # superlet-SST (best period localisation), each rank-normalised, combined by their geometric mean on the
          # grid of the wavelet (the superlet principle applied across methods)
          wv <- analyse_wavelet(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj, omega_nr = omega_nr[1], verbose = FALSE)
          su <- analyse_superlet_sst(d, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj)
          w <- if (is.null(window)) { 20 * lowerPeriod } else { window }
          e <- analyse_chirp_eha(d, window = w, chirp = 0, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, dj = dj)
          ongrid <- function(P, ax1, ax2) {                 # interpolate a map onto the wavelet grid
            P1 <- matrix(NA_real_, nrow(P), length(wv$axis.1))
            for (i in seq_len(nrow(P))) {
              P1[i, ] <- stats::approx(ax1, P[i, ], wv$axis.1, rule = 2)$y
            }
            P2 <- matrix(NA_real_, length(wv$axis.2), ncol(P1))
            for (j in seq_len(ncol(P1))) {
              P2[, j] <- stats::approx(ax2, P1[, j], wv$axis.2, rule = 2)$y
            }
            P2
          }
          rk <- function(P) {
            P[!is.finite(P)] <- 0
            matrix(rank(P) / length(P), nrow(P), ncol(P))
          }
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
  lp_t <- tuned$axis.2
  ax_t <- tuned$axis.1
  lp_z <- seq(min(lp_t) + log2(min(rz)), max(lp_t) + log2(max(rz)), by = dj)
  to_depth <- function(P) {
    Pz <- matrix(NA_real_, nrow(P), length(tz))                                  # periods(kyr) x depth
    for (i in seq_len(nrow(P))) {
      Pz[i, ] <- stats::approx(ax_t, P[i, ], tz, rule = 2)$y
    }
    out <- matrix(NA_real_, length(lp_z), n)
    for (j in seq_len(n)) {
      out[, j] <- stats::approx(lp_t + log2(rz[j]), Pz[, j], lp_z, rule = 1)$y
    }
    out
  }
  Pd <- to_depth(tuned$Power)
  kb <- max(3, 2 * floor(0.5 / dj) + 1)
  Pn <- matrix(NA_real_, nrow(Pd), ncol(Pd))
  for (j in seq_len(ncol(Pd))) {
    v <- Pd[, j]
    ok <- is.finite(v) & v > 0
    if (sum(ok) < 10) {
      Pn[, j] <- v * NA
    } else {
      bg <- rep(NA_real_, length(v))
      bg[ok] <- exp(stats::lowess(lp_z[ok], log(v[ok]), f = 0.4, iter = 2)$y)
      Pn[, j] <- v / bg
    }
  }
  depth_obj <- tuned
  depth_obj$Power <- Pd
  depth_obj$Period <- 2^lp_z
  depth_obj$axis.2 <- lp_z
  depth_obj$axis.1 <- z
  depth_obj$Wave <- sqrt(pmax(Pd, 0)) + 0i
  depth_obj$Ampl <- sqrt(pmax(Pd, 0))
  depth_obj$Phase <- NULL
  depth_obj$Scale <- depth_obj$Period / depth_obj$fourier_factor
  depth_obj$nr <- length(lp_z)
  depth_obj$nc <- n
  depth_obj$coi.1 <- NULL
  depth_obj$coi.2 <- NULL
  depth_obj$x <- z
  depth_obj$y <- x
  depth_obj$Power.avg <- rowMeans(Pd, na.rm = TRUE)
  depth_obj$dt <- dz
  out <- list(tuned = tuned, depth = depth_obj, Power_norm = Pn, transform = transform, track = track, rate = rz, time = tz)
  if (n_sim > 0) {
    set.seed(seed)
    xc <- x - mean(x)
    rho <- stats::cor(xc[-1], xc[-n])
    sims <- array(NA_real_, c(dim(Pn), n_sim))
    sparse <- transform %in% c("superlet_sst", "wavelet_sst", "chirplet_sst")
    for (i in seq_len(n_sim)) {
      e <- as.numeric(stats::arima.sim(list(ar = rho), n = n, sd = sqrt(1 - rho^2)))
      e <- e / stats::sd(e) * stats::sd(xc)
      ts <- run(e)
      Ps <- to_depth(ts$Power)
      if (sparse) {
        Pns <- Ps
      } else {
        Pns <- matrix(NA_real_, nrow(Ps), ncol(Ps))
        for (j in seq_len(ncol(Ps))) {
          v <- Ps[, j]
          ok <- is.finite(v) & v > 0
          if (sum(ok) < 10) {
            Pns[, j] <- v * NA
          } else {
            bg <- rep(NA_real_, length(v))
            bg[ok] <- exp(stats::lowess(lp_z[ok], log(v[ok]), f = 0.4, iter = 2)$y)
            Pns[, j] <- v / bg
          }
        }
      }
      sims[, , i] <- Pns
    }
    # the AR1 surrogates are stationary in depth, so for every period (row) the surrogate values of all positions and
    # all surrogates form one null distribution: far more stable than per cell with a few surrogates. Squeezed maps
    # are mostly empty, so there a ridge is compared with the ridges noise produces (the positive cells).
    obs <- if (sparse) { Pd } else { Pn }
    q95 <- rep(NA_real_, nrow(obs))
    pv <- matrix(NA_real_, nrow(obs), ncol(obs))
    for (r in seq_len(nrow(obs))) {
      v <- as.numeric(sims[r, , ])
      v <- v[is.finite(v)]
      if (sparse) {
        v <- v[v > 0]
      }
      if (length(v) < 20) {
        next
      }
      q95[r] <- stats::quantile(v, 0.95, names = FALSE)
      Fr <- stats::ecdf(v)
      pv[r, ] <- 1 - Fr(obs[r, ]) + 1 / length(v)
    }
    pv[!is.finite(obs)] <- NA
    if (sparse) {
      pv[is.finite(obs) & obs <= 0] <- 1
    }
    out$Power_sig <- obs / pmax(q95, 1e-12 * max(obs, na.rm = TRUE))
    out$Power_pval <- pmin(pv, 1)
  }
  class(out) <- "wr_dechirped"
  out
}

#' Plot a de-chirped scalogram
#'
#' @param x Result of \code{\link{analyse_dechirped}} (list of class \code{wr_dechirped}). Required, no default.
#' @param style What to draw (character), one of:
#' \itemize{
#' \item \code{"outline"}: the de-chirped depth map with the cells that are significant against the warped
#' surrogates outlined in white (needs \code{n_sim > 0} in \code{\link{analyse_dechirped}});
#' \item \code{"significance"}: the ratio to the 95th percentile of the warped surrogates (log10, diverging colour
#' scale; needs \code{n_sim > 0});
#' \item \code{"plain"}: the de-chirped depth map;
#' \item \code{"tuned"}: the transform in the time of the track (periods in kyr).
#' }
#' Without surrogates \code{"outline"} and \code{"significance"} fall back to \code{"plain"} with a message.
#' \code{Default="outline"} (the first element of \code{c("outline", "significance", "plain", "tuned")}).
#' @param show Astronomical periods in kyr to draw as dotted lines (numeric vector, e.g. \code{c(405, 100, 41)}):
#' in the depth styles at the period they have at the tracked rate (period x local rate), in the \code{"tuned"}
#' style as horizontal lines. \code{NULL} draws none. \code{Default=NULL}.
#' @param level Significance level of the outline (numeric, between 0 and 1): cells with a p-value of at most
#' \code{1 - level} are outlined. Only used when \code{style="outline"}. \code{Default=0.95}.
#' @param palette Name of a \code{\link[grDevices]{hcl.colors}} palette for the colour scale, a character string, for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}; run \code{grDevices::hcl.pals()} to list all
#' 115 names. Low values get the first colour of the palette. Names of the viridis family (\code{"viridis"},
#' \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}) give the colours of the
#' viridis package, as in the wavelet plots. \code{Default="viridis"}.
#' @param main Title of the plot (character string); \code{NULL} gives an automatic title that names the transform
#' and the style. \code{Default=NULL}.
#' @param xlab,ylab Labels of the depth and period axes (character strings) of the depth styles; the
#' \code{"tuned"} style always uses "Time of the track (kyr)" and "Period (kyr)". Defaults: \code{xlab}:
#' \code{Default="Depth"}; \code{ylab}: \code{Default="Period"}.
#' @param scale Colour scale as in all image plots of the package (character): \code{"auto"}, \code{"log"},
#' \code{"quantile"}, \code{"sqrt"} or \code{"linear"} (see \code{\link{plot_wavelet}}). The significance style
#' always uses a fixed diverging scale centred on the 95th percentile of the surrogates. \code{Default="auto"}.
#' @param outline_lwd Line width of the significance outline (numeric, > 0). \code{Default=0.6}.
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
    M <- x$tuned$Power
    if (!sparse) {
      M <- M / rowMeans(M, na.rm = TRUE)
    }     # flatten the red-noise slope
    ax1 <- x$tuned$axis.1
    ax2 <- x$tuned$axis.2
    if (is.null(main)) {
      main <- paste(x$transform, "of the record in the time of the track (tuned): cycles should be horizontal lines")
    }
    xl <- "Time of the track (kyr)"
    yl <- "Period (kyr)"
  } else {
    M <- if (style == "significance") { log10(pmax(x$Power_sig, 1e-6)) } else if (sparse) { x$depth$Power } else { x$Power_norm }
    ax1 <- x$depth$axis.1
    ax2 <- x$depth$axis.2
    if (is.null(main)) {
      main <- switch(style,
          significance = paste(x$transform, "de-chirped by the track, ratio to the 95th percentile of warped surrogates (log10; > 0 significant)"),
          outline = paste(x$transform, "de-chirped by the track, cells significant against warped surrogates outlined"),
          plain = paste(x$transform, "de-chirped by the track (lines against the continuum of each position)"))
    }
    xl <- xlab
    yl <- ylab
  }
  # empty space outside the mapped period range is drawn in the lowest colour, not left white
  M <- .wr_fill_na(M)
  if (style == "significance") {
    # diverging scale with the threshold (the 95th percentile of the surrogates, log10 ratio 0) in the middle:
    # blue below (down to a tenth of the threshold), red above
    .wr_image_signif(ax1, ax2, M, yaxt = "n", xlab = xl, ylab = yl, main = main, xaxs = "i")
  } else {
    .wr_image(ax1, ax2, M, scale = scale, palette = palette, sparse = if (sparse) { TRUE } else { NULL }, fill_na = FALSE,
        yaxt = "n", xlab = xl, ylab = yl, main = main, xaxs = "i")
  }
  ticks <- 2^(floor(min(ax2)):ceiling(max(ax2)))
  graphics::axis(2, at = log2(ticks), labels = ticks)
  if (style == "outline") {
    sig <- (x$Power_pval <= 1 - level) * 1
    sig[is.na(sig)] <- 0
    sig[is.na(x$Power_pval)] <- 0
    graphics::contour(ax1, ax2, t(sig), levels = 0.5, drawlabels = FALSE, add = TRUE, col = "white", lwd = outline_lwd)
  }
  if (!is.null(show)) {
    if (style == "tuned") {
      graphics::abline(h = log2(show), col = "white", lty = 3, lwd = 0.9)
    } else {
      for (p in show) {
        graphics::lines(x$depth$axis.1, log2(p * x$rate), col = if (style == "significance") { "black" } else { "white" }, lwd = 0.9, lty = 3)
      }
    }
  }
  invisible(M)
}
