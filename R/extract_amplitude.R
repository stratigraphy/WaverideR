#' @title Extract the instantaneous amplitude (envelope) of a signal using the CWT
#'
#' @description
#' The wavelet equivalent of the Hilbert transform. A complex (analytic) Morlet
#' wavelet is used, so the wavelet coefficients already are the analytic
#' representation of a real signal; summing them back up gives the analytic
#' signal \eqn{z(t)} whose real part is the input and whose modulus is the
#' instantaneous amplitude.
#'
#' Called with only a signal, it returns the envelope of the whole record and
#' behaves like a Hilbert transform (median ratio 1.000 on a ZB18a eccentricity window), but computed scale by
#' scale, so a period range can be
#' imposed and the ends of the record do not have to be wrapped around.
#'
#' Given \code{cycle}, it instead returns the amplitude of that one cycle: the
#' analytic signal is rebuilt only from the scales inside \code{band} around the
#' given period, so other cycles, noise and trends do not enter the amplitude.
#' This is the amplitude-modulation curve used in cyclostratigraphy.
#'
#' Calibration is absolute in both modes, with no rescaling against peak
#' statistics: a signal of constant amplitude 3 returns 3, independent of dt,
#' dj, omega_nr and the scale range. In whole-record mode the reconstruction
#' constant is fixed by requiring \eqn{\mathrm{Re}(z)} to reproduce the input;
#' in cycle mode it is the analytic response of the wavelet family to a unit
#' cosine at the cycle's frequency.
#'
#' @param signal Matrix/data frame: column 1 = depth/time, column 2 = proxy.
#' @param cycle What to measure the amplitude of. \code{NULL} (default) = the
#'   whole signal, i.e. the Hilbert-equivalent envelope. A single number = the
#'   cycle of that period (e.g. \code{cycle = 405}). A 2-column matrix
#'   (depth/time, period) = a tracked period curve, e.g. from
#'   \code{track_period_wavelet()} + \code{completed_series()}, for records
#'   whose sedimentation rate varies. \code{"ridge"} = follow the strongest
#'   cycle at each depth (only for records dominated by one cycle).
#' @param band Period band kept around \code{cycle}, as factors of it.
#'   \code{Default = c(0.5, 2)}. Narrow it (e.g. \code{c(0.8, 1.2)}) for noisy
#'   records or cycles close together; a wide band is needed when the cycle is
#'   strongly modulated. Ignored when \code{cycle = NULL}.
#' @param periods Period range analysed, \code{c(lower, upper)}. Default:
#'   \code{2*dt} to half the record. Use it to keep noise or a trend out of the
#'   whole-record envelope, or to stop \code{cycle = "ridge"} wandering.
#' @param omega_nr Cycles in the Morlet wavelet; >= 5 keeps it analytic.
#'   \code{Default = 6}.
#' @param dj Scale spacing. \code{Default = 1/100}.
#' @param pad End treatment: \code{"AR"} (default) extends the record with a
#'   forward/backward linear-prediction forecast so the wavelet does not see a
#'   step to zero, \code{"reflect"} mirrors it, \code{"zero"} is the classic
#'   zero padding. Falls back automatically if the AR fit fails.
#' @param mask_coi Return NA inside the cone of influence. \code{Default = FALSE}.
#' @param genplot Plot the signal with its +/- envelope, and the envelope.
#' @param ver_results Print the median ratio to the Hilbert envelope (outside the cone of influence).
#' @param keep_editable Keep plot parameters editable. \code{Default = FALSE}.
#' @param ... Ignored, with a note; lets old calls using \code{pts},
#'   \code{remean}, \code{tracked_cycle_curve}, \code{lowerPeriod} etc. keep
#'   working (the old names are mapped onto the new ones where possible).
#'
#' @return Matrix, column 1 = depth/time, column 2 = amplitude. In cycle mode
#'   the period actually used is attached as \code{attr(x, "period")}.
#'
#' @examples
#' \donttest{
#' # envelope of the whole record - nothing else to set (insolation-like signal, time in kyr)
#' amp <- extract_amplitude(astrosignal_example, genplot = TRUE)
#'
#' # amplitude of one cycle only: the precession band (~21 kyr) of the example signal;
#' # its amplitude is the eccentricity modulation
#' amp21 <- extract_amplitude(astrosignal_example, cycle = 21, band = c(0.75, 1.35),
#'                            genplot = TRUE)
#'
#' # depth domain with a varying sedimentation rate: hand it the tracked 405-kyr period curve
#' # (magnetic susceptibility of Pas et al. 2018, tracked 405-kyr cycle in metres)
#' mag_wt <- analyse_wavelet(mag, dj = 1/100, lowerPeriod = 0.1, upperPeriod = 254,
#'                           verbose = FALSE, omega_nr = 10)
#' mag_track <- completed_series(wavelet = mag_wt, tracked_curve = mag_track_solution,
#'                               period_up = 1.2, period_down = 0.8, extrapolate = TRUE,
#'                               genplot = FALSE)
#' amp405 <- extract_amplitude(mag, cycle = mag_track, band = c(0.8, 1.25), genplot = TRUE)
#' }
#' @export
#' @importFrom stats fft approx median na.omit ar predict sd
extract_amplitude <- function(signal = NULL,
                              cycle = NULL,
                              band = c(0.5, 2),
                              periods = NULL,
                              omega_nr = 6,
                              dj = 1/100,
                              pad = c("AR", "reflect", "zero"),
                              mask_coi = FALSE,
                              genplot = FALSE,
                              ver_results = FALSE,
                              keep_editable = FALSE,
                              ...) {

  ## ---- 0. accept the old argument names ------------------------------------
  dots <- list(...)
  if (length(dots)) {
    if (is.null(cycle) && !is.null(dots$tracked_cycle_curve)) {
      cycle <- dots$tracked_cycle_curve
      if (!is.null(dots$extract_cycle) && !is.null(dots$tracked_cycle_period))
        cycle[, 2] <- cycle[, 2] * (dots$extract_cycle / dots$tracked_cycle_period)
    }
    if (is.null(cycle) && !is.null(dots$extract_cycle)) cycle <- dots$extract_cycle
    if (is.null(periods) && (!is.null(dots$lowerPeriod) || !is.null(dots$upperPeriod)))
      periods <- c(if (is.null(dots$lowerPeriod)) NA else dots$lowerPeriod,
                   if (is.null(dots$upperPeriod)) NA else dots$upperPeriod)
    if (!is.null(dots$period_down) || !is.null(dots$period_up))
      band <- c(if (is.null(dots$period_down)) band[1] else dots$period_down,
                if (is.null(dots$period_up))   band[2] else dots$period_up)
    old <- intersect(names(dots), c("pts", "remean", "add_mean", "wavelet", "method"))
    if (length(old))
      message("extract_amplitude(): argument(s) ", paste(old, collapse = ", "),
              " are no longer used and were ignored.")
  }
  pad <- match.arg(pad)
  if (omega_nr < 5)
    warning("omega_nr < 5: the Morlet wavelet is no longer analytic; ",
            "amplitudes will be biased. Use omega_nr >= 5.")
  band <- sort(as.numeric(band))

  ## ---- 1. regular axis ------------------------------------------------------
  dat <- stats::na.omit(as.matrix(signal))
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dt  <- stats::median(diff(dat[, 1]))
  xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
  dat <- as.data.frame(stats::approx(dat[, 1], dat[, 2], xout, method = "linear"))
  x_axis <- dat[, 1]; x_raw <- dat[, 2]
  mean_signal <- mean(x_raw)
  x <- x_raw - mean_signal
  nx <- length(x)
  span <- x_axis[nx] - x_axis[1]

  lowerPeriod <- 2 * dt
  upperPeriod <- span / 2
  if (!is.null(periods)) {
    if (length(periods) != 2) stop("'periods' must be c(lower, upper).")
    if (!is.na(periods[1])) lowerPeriod <- periods[1]
    if (!is.na(periods[2])) upperPeriod <- periods[2]
  }
  # a named cycle must fit inside the analysed range
  if (is.numeric(cycle) && length(cycle) == 1) {
    lowerPeriod <- min(lowerPeriod, cycle * band[1] * 0.9)
    upperPeriod <- max(upperPeriod, cycle * band[2] * 1.1)
    lowerPeriod <- max(lowerPeriod, 2 * dt)
  }
  if (upperPeriod <= lowerPeriod) stop("'periods' upper must exceed lower.")

  ## ---- 2. CWT, Torrence & Compo Morlet, input not standardised --------------
  omega0 <- omega_nr
  ff <- (2 * pi) / omega0                       # period = ff * scale
  J <- as.integer(log2((upperPeriod/ff)/(lowerPeriod/ff)) / dj)
  scales <- (lowerPeriod/ff) * 2^((0:J) * dj)
  ns <- length(scales)
  per <- ff * scales

  pot2 <- trunc(log2(nx) + 0.5)
  pad.length <- 2^(pot2 + 1) - nx
  pad.l <- floor(pad.length/2); pad.r <- pad.length - pad.l

  ar_extend <- function(v, k, backwards = FALSE) {
    if (k <= 0) return(numeric(0))
    vv <- if (backwards) rev(v) else v
    om <- max(1, min(100, floor(length(vv)/4)))
    fit <- try(stats::ar(vv, order.max = om, aic = TRUE), silent = TRUE)
    if (inherits(fit, "try-error")) stop("ar failed")
    p <- as.numeric(stats::predict(fit, n.ahead = k)$pred)
    if (anyNA(p)) stop("ar failed")
    if (backwards) rev(p) else p
  }
  reflect_pad <- function() {
    lf <- if (pad.l > 0) rev(x[2:min(pad.l + 1, nx)]) else numeric(0)
    if (length(lf) < pad.l) lf <- c(rep(lf[1], pad.l - length(lf)), lf)
    rt <- if (pad.r > 0) rev(x[max(1, nx - pad.r):(nx - 1)]) else numeric(0)
    if (length(rt) < pad.r) rt <- c(rt, rep(rt[length(rt)], pad.r - length(rt)))
    list(s = c(lf, x, rt), off = pad.l)
  }
  padded <- NULL
  if (pad == "AR") {
    padded <- tryCatch(list(s = c(ar_extend(x, pad.l, TRUE), x, ar_extend(x, pad.r)),
                            off = pad.l),
                       error = function(e) {
                         warning("AR padding failed; using reflective padding."); NULL })
    if (!is.null(padded) && length(padded$s) != nx + pad.length) padded <- NULL
  }
  if (is.null(padded) && pad %in% c("AR", "reflect"))
    padded <- tryCatch(reflect_pad(), error = function(e) NULL)
  if (is.null(padded)) padded <- list(s = c(x, rep(0, pad.length)), off = 0)

  N <- length(padded$s); off <- padded$off
  ok <- 1:floor(N/2); ok <- ok * (2*pi)/(N*dt)
  ok <- c(0, ok, -ok[floor((N - 1)/2):1])
  fx <- stats::fft(padded$s)

  Wave <- matrix(0 + 0i, nrow = ns, ncol = nx)
  for (i in seq_len(ns)) {
    nf <- pi^(1/4) * sqrt(2 * scales[i]/dt)
    dau <- nf * exp(-((scales[i]*ok - omega0)^2/2) * (ok > 0)) * (ok > 0)
    Wave[i, ] <- (stats::fft(fx * dau, inverse = TRUE)/N)[(off + 1):(off + nx)]
  }
  Wn <- Wave / sqrt(scales)              # the 1/sqrt(s) of the inverse transform

  ## ---- 3. which scales contribute at each depth -----------------------------
  centre <- NULL
  if (is.null(cycle)) {
    mode <- "whole"
  } else if (is.character(cycle) && cycle[1] == "ridge") {
    mode <- "cycle"
    k <- apply(Mod(Wn), 2, which.max)
    lp <- log(per); centre <- per[k]
    okr <- k > 1 & k < ns
    if (any(okr)) {                      # parabolic refinement in log-period
      j <- k[okr]; cl <- which(okr)
      y1 <- Mod(Wn)[cbind(j - 1, cl)]; y2 <- Mod(Wn)[cbind(j, cl)]
      y3 <- Mod(Wn)[cbind(j + 1, cl)]
      den <- y1 - 2*y2 + y3
      dd <- pmax(pmin(ifelse(den != 0, 0.5*(y1 - y3)/den, 0), 1), -1)
      centre[okr] <- exp(lp[j] + dd*(lp[j + 1] - lp[j - 1])/2)
    }
  } else if (is.numeric(cycle) && length(cycle) == 1) {
    mode <- "cycle"; centre <- rep(cycle, nx)
  } else {
    tc <- stats::na.omit(as.matrix(cycle))
    if (ncol(tc) < 2) stop("'cycle' must be NULL, a period, a 2-column curve, or \"ridge\".")
    mode <- "cycle"
    centre <- stats::approx(tc[, 1], tc[, 2], xout = x_axis,
                            method = "linear", rule = 2)$y
  }

  ## ---- 4. amplitude ---------------------------------------------------------
  if (mode == "whole") {
    # analytic signal from all scales; the reconstruction constant is fixed by
    # requiring Re(K*z) to reproduce the input signal (exact, and independent of
    # dj / omega_nr / scale range, so no tabulated C_delta is needed)
    z <- colSums(Wn)
    K <- sum(x * Re(z)) / sum(Re(z)^2)
    ampl <- Mod(K * z)
    rec  <- Re(K * z)
    fit  <- 1 - sum((x - rec)^2)/sum(x^2)        # how much variance is captured
    if (fit < 0.9)
      message(sprintf(paste0("The analysed band reproduces %.0f%% of the signal variance; ",
                             "widen 'periods' if the envelope looks too low."), 100*fit))
  } else {
    # band around the cycle; calibrated with the wavelet family's exact response
    # to a unit cosine at that period:
    #   sum_s W/sqrt(s) = (A/2)*pi^(-1/4)*sqrt(2*pi/dt)*S(wc)*exp(i*wc*t)
    cal0 <- 0.5 * pi^(-1/4) * sqrt(2*pi/dt)
    ampl <- numeric(nx)
    const_band <- length(unique(centre)) == 1
    if (const_band) {
      idx <- which(per >= centre[1]*band[1] & per <= centre[1]*band[2])
      if (!length(idx)) idx <- which.min(abs(per - centre[1]))
      S <- sum(exp(-(scales[idx]*(2*pi/centre[1]) - omega0)^2/2))
      ampl <- Mod(colSums(Wn[idx, , drop = FALSE])) / (cal0 * S)
    } else {
      for (i in seq_len(nx)) {
        wc <- 2*pi/centre[i]
        idx <- which(per >= centre[i]*band[1] & per <= centre[i]*band[2])
        if (!length(idx)) idx <- which.min(abs(per - centre[i]))
        S <- sum(exp(-(scales[idx]*wc - omega0)^2/2))
        ampl[i] <- Mod(sum(Wn[idx, i])) / (cal0 * S)
      }
    }
  }

  ## ---- 5. cone of influence -------------------------------------------------
  coi <- ff * sqrt(2) * dt *
    c(1e-5, 1:((nx + 1)/2 - 1), rev(1:(nx/2 - 1)), 1e-5)[seq_len(nx)]
  ref <- if (mode == "cycle") centre else rep(upperPeriod, nx)
  in.coi <- ref > coi
  if (mask_coi) ampl[in.coi] <- NA

  out <- cbind(x_axis, ampl)
  colnames(out) <- c("x_axis", "amplitude")
  if (mode == "cycle") attr(out, "period") <- centre

  ## ---- 6. optional checks and plot ------------------------------------------
  if (ver_results) {
    hb <- .ea_hilbert(cbind(x_axis, x_raw))
    g <- !in.coi & !is.na(ampl)
    cat(sprintf("median(CWT amplitude / Hilbert amplitude) outside the COI = %.3f\n",
                stats::median(ampl[g]/hb[g, 2])))
    if (mode == "cycle")
      cat("  (in cycle mode a ratio below 1 is expected: the Hilbert envelope also\n",
          "   contains the other cycles, the noise and the trend.)\n")
  }

  if (genplot) {
    if (!keep_editable) { oldpar <- graphics::par(no.readonly = TRUE)
                          on.exit(graphics::par(oldpar)) }
    graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    plot(x_axis, x_raw, type = "l", col = "grey45", xlab = "",
         ylab = "signal",
         main = if (mode == "whole") "signal and envelope"
                else sprintf("signal and amplitude of the %s cycle",
                             if (length(unique(centre)) == 1)
                               signif(centre[1], 4) else "tracked"))
    graphics::lines(x_axis, mean_signal + ampl, col = "red", lwd = 2)
    graphics::lines(x_axis, mean_signal - ampl, col = "red", lwd = 2)
    plot(x_axis, ampl, type = "l", lwd = 2, xlab = "depth/time", ylab = "amplitude")
  }

  out
}

# Hilbert envelope by FFT (used only for ver_results)
.ea_hilbert <- function(data) {                    # used only for ver_results
  x <- data[, 2] - mean(data[, 2]); n <- length(x); h <- rep(0, n)
  if (n %% 2 == 0) { h[1] <- 1; h[n/2 + 1] <- 1; h[2:(n/2)] <- 2
  } else            { h[1] <- 1; h[2:((n + 1)/2)] <- 2 }
  cbind(data[, 1], Mod(stats::fft(stats::fft(x) * h, inverse = TRUE)/n))
}
