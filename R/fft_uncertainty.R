#' @title Frequency uncertainty of the Fourier transform (FFT) for a record length
#'
#' @description Calculator for the frequency uncertainty of a Fourier based spectrum. The practical
#' resolution of an FFT is the bin spacing, the inverse of the record length (\eqn{\Delta f = 1/T}).
#' The smallest frequency uncertainty that is possible according to the Heisenberg-Gabor
#' uncertainty principle (\eqn{\Delta t \Delta f \ge 1/(4\pi)}) with \eqn{\Delta t = T} is
#' \eqn{\Delta f = 1/(4\pi T)}, a factor \eqn{4\pi} smaller. The relative uncertainty at a frequency f is
#' \eqn{\Delta f / f = 1/(f T)}. For a record of N samples with sampling frequency fs, \eqn{T = N / f_s}.
#'
#' It is the Fourier counterpart of \code{\link{wavelet_uncertainty}}, which gives the uncertainty of the
#' wavelet transform. Units are free: with the record length in metres (or kyr) the frequencies are in
#' cycles per metre (or per kyr).
#'
#' @param T_length Length of the record (time or depth). Alternatively give \code{N} and \code{fs}.
#' @param N Number of samples (used with \code{fs} when \code{T_length} is not given).
#' @param fs Sampling frequency (samples per unit time or depth), \code{Default=1}.
#' @param f Optional frequency (or frequencies) of interest in cycles per unit. Gives the relative
#' uncertainty and the corresponding period range.
#' @param period Optional period(s) of interest, an alternative to \code{f} (\code{f = 1 / period}).
#'
#' @return A data frame with the record length, the practical uncertainty \code{df_practical},
#' the Heisenberg-Gabor limit \code{df_gabor}, their ratio and, when \code{f} or \code{period} is given, the
#' frequency, the relative uncertainties (\code{rel_practical}, \code{rel_gabor}) and the period ranges
#' implied by \code{f - df} and \code{f + df} for both definitions.
#'
#' @references
#' Gabor, D. (1946). Theory of communication. Part 1: The analysis of information. Journal of the Institution of
#' Electrical Engineers - Part III 93, 429-441.
#'
#' @examples
#' # 2 s record, 100 Hz: practical 0.5 Hz, Heisenberg-Gabor limit 0.040 Hz, 0.5 percent
#' fft_frequency_uncertainty(T_length = 2, f = 100)
#' # 1000 samples every 0.5 kyr and the 405 kyr eccentricity cycle
#' fft_frequency_uncertainty(N = 1000, fs = 1 / 0.5, period = 405)
#' @export
fft_frequency_uncertainty <- function(T_length = NULL, N = NULL, fs = 1, f = NULL, period = NULL) {
  if (is.null(T_length)) {
    if (is.null(N)) stop("give T_length or N (and fs)")
    T_length <- N / fs
  }
  if (any(T_length <= 0)) stop("the record length must be positive")
  if (is.null(f) && !is.null(period)) f <- 1 / period
  n <- max(length(T_length), length(f))
  T_length <- rep_len(T_length, n)
  out <- data.frame(T_length = T_length,
                    df_practical = 1 / T_length,
                    df_gabor = 1 / (4 * pi * T_length))
  out$ratio <- out$df_practical / out$df_gabor
  if (!is.null(f)) {
    out$f <- rep_len(f, n)
    out$rel_practical <- out$df_practical / out$f
    out$rel_gabor <- out$df_gabor / out$f
    inv <- function(x) ifelse(x > 0, 1 / x, Inf)
    out$period <- 1 / out$f
    out$period_min_practical <- inv(out$f + out$df_practical)
    out$period_max_practical <- inv(out$f - out$df_practical)
    out$period_min_gabor <- inv(out$f + out$df_gabor)
    out$period_max_gabor <- inv(out$f - out$df_gabor)
  }
  out
}

#' @title Calculate the uncertainty associated with a tracked cycle based on the frequency uncertainty of the FFT
#'
#' @description Fourier counterpart of \code{\link{wavelet_uncertainty}}. A tracked cycle (depth, period) gets an
#' uncertainty in frequency that equals the bin spacing of an FFT of the analysed record (\eqn{1/T}, the
#' practical resolution) or the Heisenberg-Gabor limit (\eqn{1/(4 \pi T)}), see
#' \code{\link{fft_frequency_uncertainty}}. In contrast to the wavelet the uncertainty does not depend on the
#' frequency, only on the length of the record (or of the window of an evolutive spectrum). The uncertain
#' frequencies are converted to curves of period (sedimentation rate) and with
#' \code{\link{curve2time}} to uncertain time models.
#'
#' @param tracked_cycle Curve of the tracked cycle: a matrix or data frame whose first column is depth or time and the
#' second column the period.
#' @param period_of_tracked_cycle Period of the tracked cycle in time (e.g. 405 kyr).
#' @param signal_length Length of the record (or window) that sets the frequency uncertainty, in the unit of the
#' depth column. \code{Default=NULL}: the extent of \code{tracked_cycle}.
#' @param type \code{"practical"} (1/T, default) or \code{"gabor"} (1/(4 pi T), the theoretical limit).
#' @param multi Multiple of the uncertainty that is used \code{Default=1}.
#' @param verbose Print text \code{Default=FALSE}.
#' @param genplot_time Plot the time curves with the upper and lower uncertainty \code{Default=FALSE}.
#' @param genplot_uncertainty Plot the period curves with the upper and lower uncertainty \code{Default=FALSE}.
#' @param keep_editable Keep the plot settings after plotting \code{Default=FALSE}.
#'
#' @return A data frame with the columns "depth", "period", "frequency", "uncertainty in frequency", "uncertainty in
#' frequency x_times", "time mean", "time plus x_times" and "time min x_times" (the same layout as
#' \code{\link{wavelet_uncertainty}}, with the Fourier uncertainty), and the attribute \code{"type"}.
#' @examples
#' tr <- data.frame(depth = seq(0, 100, by = 1), period = 10 + 0 * seq(0, 100, by = 1))
#' unc <- fft_uncertainty(tr, period_of_tracked_cycle = 100, type = "practical")
#' tail(unc, 2)
#' @export
fft_uncertainty <- function(tracked_cycle = NULL,
                            period_of_tracked_cycle = NULL,
                            signal_length = NULL,
                            type = c("practical", "gabor"),
                            multi = 1,
                            verbose = FALSE,
                            genplot_time = FALSE,
                            genplot_uncertainty = FALSE,
                            keep_editable = FALSE) {
  type <- match.arg(type)
  data <- as.data.frame(tracked_cycle[, c(1, 2)])
  if (is.null(signal_length)) signal_length <- diff(range(data[, 1]))
  u <- fft_frequency_uncertainty(T_length = signal_length)
  df <- if (type == "practical") u$df_practical else u$df_gabor
  data$f0 <- 1 / data[, 2]
  data$df <- rep(df, nrow(data))
  data$df_multi <- data$df * multi
  f_lo <- pmax(data$f0 - data$df_multi, .Machine$double.eps)   # lower frequency = longer period = lower sedimentation rate
  f_hi <- data$f0 + data$df_multi
  sed_min <- cbind(data[, 1], 1 / f_lo)
  sed_plus <- cbind(data[, 1], 1 / f_hi)
  sed <- cbind(data[, 1], 1 / data$f0)
  time_min <- curve2time(tracked_cycle_curve = sed_min, tracked_cycle_period = period_of_tracked_cycle, genplot = FALSE)
  time_plus <- curve2time(tracked_cycle_curve = sed_plus, tracked_cycle_period = period_of_tracked_cycle, genplot = FALSE)
  time <- curve2time(tracked_cycle_curve = sed, tracked_cycle_period = period_of_tracked_cycle, genplot = FALSE)

  if (verbose) {
    cat("total duration is: ", time[nrow(time), 2], "\n")
    cat("uncertainty of the duration is: ", time_plus[nrow(time_plus), 2] - time[nrow(time), 2], "\n")
    cat("frequency uncertainty (", type, "): ", df, " cycles per unit, defined as ", multi, " times\n")
  }
  if (genplot_time) {
    if (!keep_editable) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    wr_open_device(width = 7, height = 7)
    plot(time_min, type = "l", col = "red", ylim = c(0, max(time_plus[, 2], time_min[, 2])),
         xlab = "depth", ylab = "time")
    lines(time_plus, col = "blue")
    lines(time, col = "black", lwd = 3)
  }
  if (genplot_uncertainty) {
    if (!keep_editable) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    wr_open_device(width = 7, height = 7)
    yl <- range(1 / f_lo, 1 / f_hi)
    plot(data[, 1], 1 / data$f0, type = "l", ylim = yl, lwd = 2, xlab = "depth", ylab = "period")
    lines(data[, 1], 1 / f_lo, col = "red", lwd = 2)
    lines(data[, 1], 1 / f_hi, col = "blue", lwd = 2)
  }
  data$time_mean <- time[, 2]
  data$time_min <- time_min[, 2]
  data$time_plus <- time_plus[, 2]
  colnames(data) <- c("depth", "period", "frequency", "uncertainty in frequency",
                      "uncertainty in frequency x_times", "time mean",
                      "time min x_times", "time plus x_times")
  attr(data, "type") <- type
  data
}
