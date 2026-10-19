# =============================================================================
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
# =============================================================================

#' @title Period range implied by the frequency resolution of a record
#'
#' @param period the measured cycle period (kyr) or wavelength (m)
#' @param record_length length of the analysed record, same units as `period`
#'   (kyr for a tuned record, m for a depth series). The half-width used is
#'   `mult / record_length` in frequency.
#' @param tol alternative: a relative half-width (0.05 = +-5\%)
#' @param mult multiples of the Rayleigh resolution to use (1 = +-1 bandwidth)
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
