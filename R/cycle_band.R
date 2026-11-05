# One section, one value per cycle
#
# With a single number per astronomical cycle there is nothing to sample a
# sedimentation rate against, so the Monte Carlo consistency test of Zeeden
# et al. (2023) is not the right tool. Two things are, and both work with as
# little as one or two cycles:
#
#   k_from_cycles()  each observed cycle is read as a named term, giving
#                    k = f_observed - f_term directly (the cycle must be in
#                    kyr, i.e. a tuned record or a known sedimentation rate)
#   k_from_ratio()   the RATIO of two cycles in the same section gives k
#                    without any time scale at all: it is identical in the
#                    depth and time domain, so no tuning and no sedimentation
#                    rate is needed
#
# Uncertainty on each cycle comes from the spectrum: a peak cannot be located
# better than the frequency resolution of the record, 1/(record length),
# which cycle_band() turns into a period range.

#' @title Period range implied by the frequency resolution of a record
#'
#' @param period Measured period of the cycle (single numeric value, > 0): in kyr for a tuned record or in metres
#' for a depth series. Required, no default.
#' @param record_length Length of the analysed record, in the same unit as \code{period} (numeric, > 0; kyr for a
#' tuned record, metres for a depth series). The half-width of the band is \code{mult / record_length} in
#' frequency (the Rayleigh resolution); a warning is given when the record holds fewer than two repeats of the
#' cycle. \code{NULL} uses \code{tol} instead. \code{Default=NULL}.
#' @param tol Relative half-width of the period range (numeric, between 0 and 1; \code{0.05} means +/-5\%), used
#' only when \code{record_length} is \code{NULL}; \code{NULL} then uses 0.05. \code{Default=NULL}.
#' @param mult Number of Rayleigh bandwidths used as half-width when \code{record_length} is given (numeric, > 0;
#' \code{1} means +/-1 bandwidth). \code{Default=1}.
#' @return c(min, max) period range
#' @export
cycle_band <- function(period, record_length = NULL, tol = NULL, mult = 1){
  if (!is.null(record_length)){
    if (record_length < 2 * period){
      warning(sprintf(paste("A %.4g cycle in a %.4g-long record is not resolved",
          "(fewer than two repeats); its frequency band is",
          "open-ended, so read the result as a limit."),
          period, record_length), call. = FALSE)
    }
    f <- 1 / period
    df <- mult / record_length
    sort(c(1 / (f + df), 1 / max(f - df, .Machine$double.eps)))
  } else {
    if (is.null(tol)){
      tol <- 0.05
    }
    sort(period * c(1 - tol, 1 + tol))
  }
}
