#' @title Track a spectral peak through an evolutive MESA
#'
#' @description Follows one spectral ridge (for example the 100 ka cycle) through the window
#' positions of an evolutive MESA (\code{\link{mesa_evolutive}}). A dynamic-programming search finds
#' the path through the spectral image that has the highest summed log relative power, with the
#' frequency allowed to change by at most \code{max_jump} (as a fraction of the frequency) from one
#' position to the next and a penalty for large changes, so that the track does not hop between
#' neighbouring peaks. Optionally the track is read off the stored peaks instead of the image.
#'
#' @param x output of \code{\link{mesa_evolutive}}.
#' @param window number or length of the window to track in.
#' @param period Period (wavelength) band \code{c(lower, upper)} in which the ridge is searched,
#' in the units of the first column of the record. Give \code{period} or \code{frequency}.
#' @param frequency Frequency band \code{c(lower, upper)}, alternative to \code{period}.
#' @param max_jump Largest relative change of the frequency between neighbouring window positions
#' (default 0.1, i.e. 10 percent).
#' @param penalty Penalty (in units of log10 relative power) for a change of the frequency by
#' \code{max_jump}; smaller changes cost proportionally less. Default 0.3.
#' @param start Optional period at the first position (or last, if \code{start_at = "end"})
#' that fixes where the track begins; the search then starts at the peak nearest to it.
#' @param start_at \code{"top"} (first position, default) or \code{"end"}.
#' @param snap If TRUE (default) every position of the track is moved to the nearest local maximum
#' of the spectrum within one frequency step, so that the track follows the spectral peak and not
#' the shoulder of a broad one.
#'
#' @return A data frame of class \code{mesa_track} with \code{position}, \code{frequency},
#' \code{period} and \code{rel_power} (power relative to the maximum at that position), and the
#' attribute \code{window} giving the window length.
#'
#' @examples
#' set.seed(2)
#' t <- seq(0, 800, by = 2)
#' y <- sin(2 * pi * cumsum(1 / ifelse(t < 400, 40, 80)) * 2) + rnorm(length(t), sd = 0.4)
#' ev <- mesa_evolutive(data.frame(t, y), window = 200, order = 20)
#' tr <- track_mesa(ev, period = c(25, 120))
#' plot_mesa_evolutive(ev, track = tr, vertical = FALSE)
#' @export
track_mesa <- function(x,
                       window = 1,
                       period = NULL,
                       frequency = NULL,
                       max_jump = 0.1,
                       penalty = 0.3,
                       start = NULL,
                       start_at = c("top", "end"),
                       snap = TRUE) {
  start_at <- match.arg(start_at)
  w <- .mesa_pick_window(x, window)
  pw <- x$windows[[w]]$power
  freq <- x$frequency
  if (is.null(period) && is.null(frequency)){
    stop("give a period or a frequency band to search in")
  }
  if (is.null(frequency)){
    frequency <- sort(1 / period)
  }
  frequency <- sort(frequency)
  band <- which(freq >= frequency[1] & freq <= frequency[2] & freq > 0)
  if (length(band) < 3){
    stop("the band holds fewer than 3 frequencies; widen it")
  }
  valid <- which(colSums(is.finite(pw)) > 0)
  if (!length(valid)){
    stop("no spectra in this window")
  }
  rel <- sweep(pw, 2, apply(pw, 2, max, na.rm = TRUE), "/")
  sc <- log10(pmax(rel[band, valid, drop = FALSE], 1e-12))
  nb <- length(band)
  nt <- length(valid)
  fb <- freq[band]
  # allowed transitions: |f_j - f_i| <= max_jump * f_i
  cost <- function(i){
    j <- which(abs(fb - fb[i]) <= max_jump * fb[i])
    list(j = j, c = penalty * abs(log(fb[j] / fb[i])) / log1p(max_jump))
  }
  tr <- lapply(seq_len(nb), cost)
  dir <- if (start_at == "top") seq_len(nt) else rev(seq_len(nt))
  score <- matrix(-Inf, nb, nt)
  from <- matrix(NA_integer_, nb, nt)
  s0 <- sc[, dir[1]]
  if (!is.null(start)){
    i0 <- which.min(abs(fb - 1 / start))
    near <- abs(fb - fb[i0]) <= max_jump * fb[i0]
    s0[!near] <- -Inf
  }
  score[, 1] <- s0
  for (tt in 2:nt){
    prev <- score[, tt - 1]
    for (i in seq_len(nb)){
      cand <- prev[tr[[i]]$j] - tr[[i]]$c
      k <- which.max(cand)
      score[i, tt] <- cand[k] + sc[i, dir[tt]]
      from[i, tt] <- tr[[i]]$j[k]
    }
  }
  path <- integer(nt)
  path[nt] <- which.max(score[, nt])
  for (tt in nt:2){
    path[tt - 1] <- from[path[tt], tt]
  }
  idx <- band[path]
  if (snap){
    pws <- pw[, valid[dir], drop = FALSE]
    for (tt in seq_len(nt)){
      lo <- max(1, idx[tt] - 2)
      hi <- min(length(freq), idx[tt] + 2)
      idx[tt] <- lo - 1 + which.max(pws[lo:hi, tt])
    }
  }
  pos <- x$position[valid[dir]]
  relp <- vapply(seq_len(nt), function(tt) rel[idx[tt], valid[dir][tt]], numeric(1))
  out <- data.frame(position = pos, frequency = freq[idx], period = 1 / freq[idx],
                    rel_power = relp)
  out <- out[order(out$position), ]
  rownames(out) <- NULL
  attr(out, "window") <- x$windows[[w]]$window
  class(out) <- c("mesa_track", "data.frame")
  out
}
