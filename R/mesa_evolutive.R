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
#' @param data Record to analyse: a two-column matrix or data frame with depth or time in column 1 and the proxy
#' value in column 2; an unevenly sampled record is interpolated when \code{regularise = TRUE}. Required
#' (\code{Default=NULL}).
#' @param window Length of the moving window, in the unit of the first column of \code{data} (numeric, each value
#' at least 8 sampling intervals). A vector gives several window lengths at once, e.g. \code{c(30, 40, 50)}.
#' \code{NULL} uses a fifth of the record length. \code{Default=NULL}.
#' @param step Distance between the centres of successive window positions, in the unit of the first column of
#' \code{data} (numeric, > 0). \code{NULL} uses a tenth of the shortest window (at least one sampling interval).
#' \code{Default=NULL}.
#' @param order AR order used in every window (integer, >= 1). \code{NULL} chooses the order in each window with
#' \code{criterion}; a fixed order gives windows of the same spectral smoothness, which are easier to compare along
#' the record. \code{Default=NULL}.
#' @param order_max Largest AR order considered when the order is selected (integer, >= 1); \code{NULL} uses
#' \code{min(floor(N / 3), 300)} for N samples in the window (see \code{\link{mesa_spectrum}}). \code{Default=NULL}.
#' @param criterion Order-selection criterion (character), one of \code{"fpe"} (final prediction error),
#' \code{"aic"} (Akaike information criterion), \code{"mdl"} (minimum description length) or \code{"cat"}
#' (criterion autoregressive transfer); ignored when \code{order} is given (see \code{\link{mesa_spectrum}}).
#' \code{Default=c("fpe", "aic", "mdl", "cat")}, i.e. \code{"fpe"}.
#' @param demean Remove the mean of the record, and again of every window, before the analysis (logical, TRUE or
#' FALSE). \code{Default=TRUE}.
#' @param detrend Remove a linear trend of the whole record before the analysis instead of only the mean; every
#' window is then demeaned (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param regularise Linearly interpolate an unevenly sampled record to its median sampling interval (logical, TRUE
#' or FALSE); with \code{FALSE} an unevenly sampled record gives an error. \code{Default=TRUE}.
#' @param pad Filling of the part of a window that lies outside the record (character): \code{"reflect"} (mirror
#' the record at its edge), \code{"mean"} (the mean of the part of the window inside the record), \code{"zero"}
#' (zero, i.e. the mean of the record after demeaning) or \code{"none"} (the window is cut at the edge of the record;
#' windows with less than \code{min_cover} of their length inside the record are left empty).
#' \code{Default=c("reflect", "mean", "zero", "none")}, i.e. \code{"reflect"}.
#' @param min_cover Only for \code{pad = "none"}: smallest fraction of a window that must lie inside the record for
#' its spectrum to be computed (numeric, between 0 and 1). \code{Default=0.5}.
#' @param nfreq Number of frequencies, evenly spaced between \code{fmin} and \code{fmax}, of every spectrum (integer,
#' > 1). \code{Default=512}.
#' @param fmin Lowest frequency of the spectra, in cycles per unit of the first column of \code{data} (numeric,
#' >= 0 and < \code{fmax}). \code{Default=0}.
#' @param fmax Highest frequency of the spectra, in cycles per unit of the first column of \code{data} (numeric;
#' values above the Nyquist frequency are reduced to it). \code{NULL} uses the Nyquist frequency.
#' \code{Default=NULL}.
#' @param n_peaks Maximum number of spectral peaks stored per window position, the strongest ones (integer, > 0).
#' \code{Default=8}.
#' @param min_prominence Smallest prominence of a stored peak (numeric, >= 1): the prominence is the peak power divided
#' by the higher of its two neighbouring minima; weaker peaks are not stored. \code{Default=1.5}.
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
      pad = pad, criterion = if (is.null(order)) { criterion } else { "fixed" })
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
        if (all(is.na(w$order))) { "-" } else { paste(range(w$order, na.rm = TRUE), collapse = "-") },
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
#' @param x Evolutive MESA, the output of \code{\link{mesa_evolutive}}. Required (no default).
#' @param window Which window length of \code{x} is averaged (numeric): either the window length itself (as given to
#' \code{mesa_evolutive}) or its number, i.e. its position in \code{x$windows}. \code{Default=1}.
#' @param normalise Divide the spectrum of every window position by its maximum before averaging, so that all
#' positions weigh equally (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @return data frame with \code{frequency}, \code{period} and \code{amplitude}.
#' @export
mesa_average <- function(x, window = 1, normalise = TRUE){
  w <- .mesa_pick_window(x, window)
  pw <- x$windows[[w]]$power
  if (normalise){
    cmax <- numeric(ncol(pw))
    for (j in seq_len(ncol(pw))) {
      cmax[j] <- max(pw[, j], na.rm = TRUE)
    }
    pw <- sweep(pw, 2, cmax, "/")
  }
  amp <- rowMeans(sqrt(pw), na.rm = TRUE)
  data.frame(frequency = x$frequency, period = ifelse(x$frequency > 0, 1 / x$frequency, Inf),
      amplitude = amp)
}

# index of a window given its number or its length
.mesa_pick_window <- function(x, window){
  lens <- numeric(length(x$windows))
  for (i in seq_along(x$windows)) {
    lens[i] <- x$windows[[i]]$window
  }
  names(lens) <- names(x$windows)
  hit <- which(abs(lens - window) < 1e-9 * max(1, abs(window)))
  if (length(hit)){
    return(hit[1])
  }
  if (window >= 1 && window <= length(lens) && window == round(window)){
    return(as.integer(window))
  }
  stop("window not found; available window lengths: ", paste(format(lens), collapse = ", "))
}
