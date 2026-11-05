#' @title Build dynamic-filter control points from a tracked period curve
#'
#' @description The \code{\link{dynamic_filter_controls}} function converts a
#' tracked period curve [depth, period] into pass band control points for
#' \code{\link{dynamic_filter}}. This replaces the interactive clicking of
#' lower/upper frequency boundaries on an evolutionary spectrum in Acycle
#' (dynamic_filter_lang.m) by a reproducible, non-interactive step. The
#' curve can come directly from the WaverideR functions
#' \code{track_period_wavelet} (first two columns: depth, period) or
#' \code{completed_series} (depth, period), so cycles already tracked in a
#' wavelet scalogram can be filtered with the Acycle dynamic filter.
#' The centre frequency is \code{1 / period}; the pass band is
#' \code{f_centre * (1 -/+ rel_bandwidth / 2)}, i.e. \code{rel_bandwidth} is
#' the half-power bandwidth relative to the centre frequency (the band is
#' symmetric in frequency, as the Acycle Gaussian response is).
#' Written for the R translation of Acycle v3.1 dynamic filtering (Mingsong
#' Li); no MATLAB equivalent.
#'
#' @param tracked Tracked period curve of the cycle (matrix or data frame whose first two columns are depth/time and
#' period in the same depth/time unit, period > 0; further columns are ignored), e.g. the output of
#' \code{\link{track_period_wavelet}} or \code{\link{completed_series}}. Rows with non-finite values are dropped,
#' rows are sorted by depth and duplicated depths are removed (first kept). Required (\code{Default=NULL}).
#' @param rel_bandwidth Half-power bandwidth of the pass band relative to the centre frequency \code{1 / period}
#' (numeric, > 0 and < 2): the band is \code{f_centre * (1 -/+ rel_bandwidth / 2)}, so larger values give a wider
#' band. \code{Default=0.2}.
#' @param n_points Maximum number of control points (integer, >= 1, or \code{NULL}). When the tracked curve has more
#' rows, it is linearly interpolated to \code{n_points} evenly spaced depths spanning the (clipped) curve
#' (\code{n_points=1} gives one point at the middle). \code{NULL} keeps all rows. \code{Default=NULL}.
#' @param data Series that will be filtered (two-column matrix or data frame, depth/time in column 1, or
#' \code{NULL}). When given, control points outside its depth range are removed, because
#' \code{\link{dynamic_filter}} requires control depths inside the data range. \code{Default=NULL}.
#'
#' @author Michiel Arts
#'
#' @references
#' Arts, M., 2023. WaverideR: Extracting signals from wavelet spectra.
#' <https://CRAN.R-project.org/package=WaverideR>
#'
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' @examples
#' \donttest{
#' depth <- seq(0, 300, by = 1)
#' period <- 5 + 10 * depth / 300
#' ctrl <- dynamic_filter_controls(tracked = cbind(depth, period),
#'   rel_bandwidth = 0.3, n_points = 20)
#' head(ctrl)
#' }
#'
#' @return A matrix with columns \code{depth}, \code{f_low}, \code{f_high}
#' and \code{f_centre}, directly usable as the \code{controls} argument of
#' \code{\link{dynamic_filter}}.
#'
#' @export
#' @importFrom stats approx

dynamic_filter_controls <- function(tracked = NULL,
    rel_bandwidth = 0.2,
    n_points = NULL,
    data = NULL){
  # Internal helper functions (local to this function)

  # @title Positive finite scalar check
  # @keywords internal
  # @noRd
  dynf_positive_scalar <- function(value = NULL, name = "value"){
    if (is.null(value) || length(value) != 1 || !is.numeric(value) ||
        !is.finite(value) || value <= 0){
      stop(paste0(name, " must be a finite positive number"))
    }
    as.numeric(value)
  }

  # Main code
  if (is.null(tracked)){
    stop("tracked must be a matrix with columns [depth, period]")
  }
  tracked <- as.matrix(tracked)
  if (ncol(tracked) < 2){
    stop("tracked must have at least two columns: [depth, period]")
  }
  tracked <- matrix(as.numeric(tracked[, 1:2]), ncol = 2)
  keep <- is.finite(tracked[, 1]) & is.finite(tracked[, 2])
  tracked <- tracked[keep, , drop = FALSE]
  tracked <- tracked[order(tracked[, 1]), , drop = FALSE]
  tracked <- tracked[!duplicated(tracked[, 1]), , drop = FALSE]
  if (any(tracked[, 2] <= 0)){
    stop("tracked periods must be positive")
  }
  rel_bandwidth <- dynf_positive_scalar(rel_bandwidth, "rel_bandwidth")
  if (rel_bandwidth >= 2){
    stop("rel_bandwidth must be smaller than 2 (f_low must stay positive)")
  }
  if (!is.null(data)){
    data <- as.matrix(data)
    d_range <- range(as.numeric(data[, 1]), na.rm = TRUE)
    inside <- tracked[, 1] >= d_range[1] & tracked[, 1] <= d_range[2]
    tracked <- tracked[inside, , drop = FALSE]
  }
  if (nrow(tracked) < 1){
    stop("no tracked points left (inside the data range)")
  }
  if (!is.null(n_points)){
    n_points <- floor(dynf_positive_scalar(n_points, "n_points"))
    if (nrow(tracked) > n_points){
      if (n_points < 2){
        new_depth <- mean(range(tracked[, 1]))
      } else {
        new_depth <- seq(tracked[1, 1], tracked[nrow(tracked), 1],
            length.out = n_points)
      }
      new_period <- approx(tracked[, 1], tracked[, 2], xout = new_depth,
          method = "linear", rule = 2)$y
      tracked <- cbind(new_depth, new_period)
    }
  }
  f_centre <- 1 / tracked[, 2]
  f_low <- f_centre * (1 - rel_bandwidth / 2)
  f_high <- f_centre * (1 + rel_bandwidth / 2)
  cbind(depth = tracked[, 1], f_low = f_low, f_high = f_high,
      f_centre = (f_low + f_high) / 2)
}
