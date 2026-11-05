#' @title Track a spectral peak through an evolutive MESA
#'
#' @description Follows one spectral ridge (for example the 100 ka cycle) through the window
#' positions of an evolutive MESA (\code{\link{mesa_evolutive}}). A dynamic-programming search finds
#' the path through the spectral image that has the highest summed log relative power, with the
#' frequency allowed to change by at most \code{max_jump} (as a fraction of the frequency) from one
#' position to the next and a penalty for large changes, so that the track does not hop between
#' neighbouring peaks. Optionally the track is read off the stored peaks instead of the image.
#'
#' @param x Evolutive MESA to track in: the output of \code{\link{mesa_evolutive}} (list with the elements
#' \code{windows}, \code{frequency} and \code{position}). Required, no default.
#' @param window Which window length of \code{x} to track in (numeric): either the window length itself (in the units
#' of the first column of the record, matched exactly) or the index of the window in \code{x$windows} (integer, 1 is
#' the first window length given to \code{\link{mesa_evolutive}}). A length is matched first. \code{Default=1}.
#' @param period Period (wavelength) band in which the ridge is searched (numeric vector \code{c(lower, upper)}, in the
#' units of the first column of the record; the band must hold at least three frequencies of \code{x}). Give
#' \code{period} or \code{frequency}; \code{frequency} takes precedence when both are given. \code{Default=NULL}.
#' @param frequency Frequency band in which the ridge is searched (numeric vector \code{c(lower, upper)}, in cycles per
#' unit of the first column of the record), alternative to \code{period}. \code{Default=NULL}.
#' @param max_jump Largest relative change of the frequency between neighbouring window positions (numeric, > 0;
#' 0.1 allows 10 percent); larger values let the track follow faster changes but make hopping between neighbouring
#' peaks easier. \code{Default=0.1}.
#' @param penalty Penalty for a change of the frequency by \code{max_jump}, in units of log10 relative power
#' (numeric, >= 0); smaller changes cost proportionally less (linear in the log of the frequency ratio), 0 means no
#' penalty within \code{max_jump}. \code{Default=0.3}.
#' @param start Optional period (numeric, same units as \code{period}) at the first position (or at the last, if
#' \code{start_at = "end"}) that fixes where the track begins: only frequencies within \code{max_jump} of it are allowed
#' there. \code{NULL} lets the search start anywhere in the band. \code{Default=NULL}.
#' @param start_at Direction of the search (character, one of \code{"top"}, \code{"end"}): \code{"top"} starts at the
#' first window position, \code{"end"} at the last; this also decides where \code{start} applies.
#' \code{Default="top"}.
#' @param snap Move every position of the track to the largest value of the spectrum within two frequency steps of the
#' path (logical, TRUE or FALSE), so that the track follows the spectral peak and not the shoulder of a broad one.
#' \code{Default=TRUE}.
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
  cmax <- numeric(ncol(pw))
  for (j in seq_len(ncol(pw))) {
    cmax[j] <- max(pw[, j], na.rm = TRUE)
  }
  rel <- sweep(pw, 2, cmax, "/")
  sc <- log10(pmax(rel[band, valid, drop = FALSE], 1e-12))
  nb <- length(band)
  nt <- length(valid)
  fb <- freq[band]
  # allowed transitions: |f_j - f_i| <= max_jump * f_i
  cost <- function(i){
    j <- which(abs(fb - fb[i]) <= max_jump * fb[i])
    list(j = j, c = penalty * abs(log(fb[j] / fb[i])) / log1p(max_jump))
  }
  tr <- vector("list", nb)
  for (i in seq_len(nb)) {
    tr[[i]] <- cost(i)
  }
  dir <- if (start_at == "top") { seq_len(nt) } else { rev(seq_len(nt)) }
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
  relp <- numeric(nt)
  for (tt in seq_len(nt)) {
    relp[tt] <- rel[idx[tt], valid[dir][tt]]
  }
  out <- data.frame(position = pos, frequency = freq[idx], period = 1 / freq[idx],
      rel_power = relp)
  out <- out[order(out$position), ]
  rownames(out) <- NULL
  attr(out, "window") <- x$windows[[w]]$window
  class(out) <- c("mesa_track", "data.frame")
  out
}
