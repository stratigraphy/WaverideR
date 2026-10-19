#' @title Evolutive maximum entropy spectral analysis (evolutive MESA)
#'
#' @description Slides a window along a record and computes a maximum entropy (Burg) spectrum
#' (\code{\link{mesa_spectrum}}) in every window position, giving the change of the spectrum with
#' depth or time. Several window lengths can be given at once, because the best window depends on
#' the cycle: a short window follows changes in the cycle period closely but resolves long cycles
#' badly, a long window resolves them but smooths out changes. Near the top and bottom of the
#' record the window extends beyond the data; this part is padded (see \code{pad}), so that the
#' evolutive spectrum spans the entire record. The spectral peaks of every window position are
#' extracted and can be followed with \code{\link{track_mesa}} and plotted with
#' \code{\link{plot_mesa_evolutive}}.
#'
#' @param data Input is a time series with the first column being depth or time and the second
#' column being a proxy.
#' @param window Window length (in the units of the first column). A vector gives several windows,
#' for example \code{c(30, 40, 50)}. The default is a fifth of the record.
#' @param step Distance between the centres of successive windows. The default is a tenth of the
#' shortest window.
#' @param order AR order used in every window. The default NULL chooses the order in each window with
#' \code{criterion}; a fixed order gives windows of the same spectral smoothness, which is easier
#' to compare along the record.
#' @param order_max,criterion,demean,detrend,regularise see \code{\link{mesa_spectrum}}. The
#' mean (or trend) is removed in every window.
#' @param pad Filling of the part of a window that lies outside the record: \code{"reflect"} (mirror
#' the record at its edge; default), \code{"mean"}, \code{"zero"} or \code{"none"} (the window is
#' cut at the edge of the record; windows with less than \code{min_cover} of their length inside
#' the record are left empty).
#' @param min_cover For \code{pad = "none"}: smallest fraction of a window that must lie in the record.
#' @param nfreq Number of frequencies of the spectra.
#' @param fmin,fmax Frequency range. Defaults 0 and the Nyquist frequency.
#' @param n_peaks Maximum number of peaks stored per window position (the strongest).
#' @param min_prominence Peaks with a prominence (peak power divided by the higher of its two
#' neighbouring minima) below this value are not stored.
#'
#' @details Powers are not normalised in the object; \code{\link{plot_mesa_evolutive}} and
#' \code{\link{track_mesa}} normalise every window position by its maximum, so that weak and strong
#' intervals of the record are comparable.
#'
#' @return A list of class \code{mesa_evolutive} with
#' \describe{
#' \item{windows}{one element per window length with \code{window}, \code{power} (matrix, frequencies in rows
#' and window positions in columns), \code{order} (AR order of each position) and \code{peaks}
#' (data frame with position, frequency, period, power, relative power and prominence).}
#' \item{frequency}{frequencies of the rows.}
#' \item{position}{centres of the windows.}
#' \item{dt, range}{sampling interval and the range of the record.}
#' }
#'
#' @examples
#' # a cycle that doubles its period half way
#' set.seed(2)
#' t <- seq(0, 800, by = 2)
#' y <- sin(2 * pi * cumsum(1 / ifelse(t < 400, 40, 80)) * 2) + rnorm(length(t), sd = 0.4)
#' ev <- mesa_evolutive(data.frame(t, y), window = c(150, 250), order = 20)
#' plot_mesa_evolutive(ev, vertical = FALSE)
#'
#' @importFrom stats median
#' @export
mesa_evolutive <- function(data = NULL,
                           window = NULL,
                           step = NULL,
                           order = NULL,
                           order_max = NULL,
                           criterion = c("fpe", "aic", "mdl", "cat"),
                           demean = TRUE,
                           detrend = FALSE,
                           regularise = TRUE,
                           pad = c("reflect", "mean", "zero", "none"),
                           min_cover = 0.5,
                           nfreq = 512,
                           fmin = 0,
                           fmax = NULL,
                           n_peaks = 8,
                           min_prominence = 1.5) {
  criterion <- match.arg(criterion)
  pad <- match.arg(pad)
  prep <- .mesa_prepare(data, regularise = regularise, demean = FALSE, detrend = FALSE)
  t <- prep$t
  x <- prep$x
  dt <- prep$dt
  rng <- range(t)
  len <- diff(rng)
  if (is.null(window)){
    window <- len / 5
  }
  window <- sort(unique(window))
  if (any(window < 8 * dt)){
    stop("every window must hold at least 8 samples")
  }
  if (is.null(step)){
    step <- max(dt, min(window) / 10)
  }
  nyq <- 1 / (2 * dt)
  if (is.null(fmax)){
    fmax <- nyq
  }
  fmax <- min(fmax, nyq)
  if (fmin < 0 || fmin >= fmax){
    stop("fmin must be >= 0 and smaller than fmax")
  }
  freq <- seq(fmin, fmax, length.out = nfreq)
  pos <- seq(rng[1], rng[2], by = step)
  if (utils::tail(pos, 1) < rng[2] - 1e-9 * len){
    pos <- c(pos, rng[2])
  }
  # remove the overall mean (or trend) once; the window is demeaned again below
  if (detrend){
    x <- x - stats::lm.fit(cbind(1, seq_along(x)), x)$fitted.values
  } else if (demean){
    x <- x - mean(x)
  }
  res <- vector("list", length(window))
  for (w in seq_along(window)){
    wl <- window[w]
    pw <- matrix(NA_real_, nfreq, length(pos))
    ord <- rep(NA_integer_, length(pos))
    pk_list <- vector("list", length(pos))
    for (j in seq_along(pos)){
      seg <- .mesa_segment(x, t[1], dt, pos[j] - wl / 2, pos[j] + wl / 2, pad)
      if (pad == "none" && length(seg) < max(8, min_cover * wl / dt)){
        next
      }
      if (length(seg) < 8){
        next
      }
      if (demean || detrend){
        seg <- seg - mean(seg)
      }
      if (all(seg == 0)){
        next
      }
      fit <- .mesa_fit(seg, order = order, order_max = order_max, criterion = criterion)
      p <- .mesa_psd(fit$a, fit$noise_variance, dt, freq)
      pw[, j] <- p
      ord[j] <- fit$order
      pk <- .mesa_peaks(freq, p, n_peaks)
      pk <- pk[pk$prominence >= min_prominence, , drop = FALSE]
      if (nrow(pk)){
        pk$position <- pos[j]
        pk$rel_power <- pk$power / max(p)
        pk_list[[j]] <- pk[, c("position", "frequency", "period", "power", "rel_power", "prominence")]
      }
    }
    pk_all <- do.call(rbind, pk_list)
    if (is.null(pk_all)){
      pk_all <- data.frame(position = numeric(0), frequency = numeric(0), period = numeric(0),
                           power = numeric(0), rel_power = numeric(0), prominence = numeric(0))
    }
    res[[w]] <- list(window = wl, power = pw, order = ord, peaks = pk_all)
  }
  names(res) <- paste0("window_", format(window, trim = TRUE))
  out <- list(windows = res, frequency = freq, position = pos, dt = dt, range = rng,
              pad = pad, criterion = if (is.null(order)) criterion else "fixed")
  class(out) <- "mesa_evolutive"
  out
}

#' @export
print.mesa_evolutive <- function(x, ...){
  cat("Evolutive MESA: record", signif(x$range[1], 6), "to", signif(x$range[2], 6),
      ", dt =", signif(x$dt, 5), ", padding:", x$pad, "\n")
  for (w in x$windows){
    cat(sprintf("  window %-8s %d positions, AR order %s, %d peaks\n", format(w$window),
                length(x$position),
                if (all(is.na(w$order))) "-" else paste(range(w$order, na.rm = TRUE), collapse = "-"),
                nrow(w$peaks)))
  }
  invisible(x)
}

#' @title Average spectrum of an evolutive MESA
#'
#' @description Averages the spectra of all window positions of an evolutive MESA
#' (\code{\link{mesa_evolutive}}): the average amplitude spectrum, shown beside the spectral
#' image in many cyclostratigraphic workflows. Amplitude is the square root of the power.
#'
#' @param x output of \code{\link{mesa_evolutive}}.
#' @param window number (position in \code{x$windows}) or length of the window.
#' @param normalise If TRUE (default) every window position is divided by its maximum before
#' averaging, so that all positions weigh equally.
#' @return data frame with \code{frequency}, \code{period} and \code{amplitude}.
#' @export
mesa_average <- function(x, window = 1, normalise = TRUE){
  w <- .mesa_pick_window(x, window)
  pw <- x$windows[[w]]$power
  if (normalise){
    pw <- sweep(pw, 2, apply(pw, 2, max, na.rm = TRUE), "/")
  }
  amp <- rowMeans(sqrt(pw), na.rm = TRUE)
  data.frame(frequency = x$frequency, period = ifelse(x$frequency > 0, 1 / x$frequency, Inf),
             amplitude = amp)
}

# index of a window given its number or its length
.mesa_pick_window <- function(x, window){
  lens <- vapply(x$windows, function(w) w$window, numeric(1))
  hit <- which(abs(lens - window) < 1e-9 * max(1, abs(window)))
  if (length(hit)){
    return(hit[1])
  }
  if (window >= 1 && window <= length(lens) && window == round(window)){
    return(as.integer(window))
  }
  stop("window not found; available window lengths: ", paste(format(lens), collapse = ", "))
}
