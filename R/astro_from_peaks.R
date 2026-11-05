#' @title Astronomical parameters from spectral peaks
#'
#' @description Wraps ctest() and ctest_earthmoon(): give the wavelength ranges
#' (m) of the cycles you identified in a depth-domain spectrum and an
#' approximate age, and get back the secular frequencies g2-g5 and s3, the
#' precession constant k, the Earth-Moon distance and the length of day, with
#' uncertainties.
#'
#' @param age Approximate age of the record in Ma (numeric): \code{c(min, max)}, or a single value, which is
#'   widened by \code{age_tol} Ma on either side. Required, no default.
#' @param e405,e125,e95 Wavelength ranges in metres of the 405-kyr, ~125-kyr and ~95-kyr eccentricity cycles
#'   identified in the depth-domain spectrum (each a numeric vector \code{c(min, max)}). All three are required,
#'   no default.
#' @param obl Wavelength range in metres of the obliquity cycle (k+s3), if identified (numeric vector
#'   \code{c(min, max)}); \code{NULL} if not identified. \code{Default=NULL}.
#' @param prec_long,prec_mid,prec_short Wavelength ranges in metres of the precession cycles, from longest to
#'   shortest (k+g5; k+g2; k+g4/k+g3), each a numeric vector \code{c(min, max)} or \code{NULL}. All three must be
#'   given for the precession test to be used. \code{Default=NULL}.
#' @param age_tol Half-width in Ma of the age range used when \code{age} is a single value (numeric, >= 0).
#'   \code{Default=5}.
#' @param nsim Number of accepted Monte Carlo samples to collect (integer, > 0); passed to \code{\link{ctest}}.
#'   \code{Default=2000}.
#' @param ptest Range of precession constants k to sample, in years (numeric vector \code{c(min, max)}); passed to
#'   \code{\link{ctest}}. \code{Default=c(10000, 28000)}.
#' @param mode Which cycles are used in the consistency test (character): \code{"auto"} picks the mode from the
#'   cycles supplied, \code{"PT"} precession and obliquity (tilt), \code{"T"} obliquity only, \code{"P"} precession
#'   only. Passed to \code{\link{ctest}}. \code{Default=c("auto", "PT", "T", "P")}, i.e. \code{"auto"}.
#' @param sampleorb How the secular frequencies g2-g5 are sampled, passed to \code{\link{ctest}} (character):
#'   \code{"U"} uniform within 2 sigma (as in the paper) or \code{"N"} normal. \code{Default="U"}.
#' @param method How the sampled precession constants are converted into an Earth-Moon state (distance, length of
#'   day) by \code{\link{ctest_earthmoon}} (character): \code{"angular_momentum"} (Earth-Moon angular-momentum
#'   conservation) or \code{"astrogeo22"} (inversion of the AstroGeo22 curve).
#'   \code{Default=c("angular_momentum", "astrogeo22")}, i.e. \code{"angular_momentum"}.
#' @param s3_uncertainty How the uncertainty of s3 is sampled, passed to \code{\link{ctest}} (character):
#'   \code{"published"} follows the published CTest code, in which the relative uncertainty is applied as an
#'   absolute one in years (s3 is effectively fixed; use it to reproduce the paper); \code{"scaled"} applies it
#'   relative to the s3 period, as for g2-g5, giving s3 a realistic spread.
#'   \code{Default=c("published", "scaled")}, i.e. \code{"published"}.
#' @param seed Seed of the random number generator (integer), passed to \code{\link{ctest}} so that the
#'   Monte Carlo results can be reproduced; \code{NULL} leaves the current random state unchanged.
#'   \code{Default=NULL}.
#' @param verbose Print progress messages of the Monte Carlo sampling (logical, TRUE or FALSE).
#'   \code{Default=TRUE}.
#' @param data_dir Folder that contains the AstroGeo22 table \code{AstroGeo22.txt} (character path);
#'   \code{NULL} uses the copy installed with the package. Passed to \code{\link{ctest_earthmoon}} and
#'   \code{\link{astrogeo22}}. \code{Default=NULL}.
#' @return object of class `astro_from_peaks`: `$frequencies` (g, s, k in both
#'   arcsec/yr and yr), `$earthmoon` (distance, day length, months, days per
#'   year), `$periods` (the astronomical periods implied), `$reference` (the
#'   AstroGeo22 model at the same age) and the underlying `$ctest`
#' @export
astro_from_peaks <- function(age, e405, e125, e95, obl = NULL,
    prec_long = NULL, prec_mid = NULL, prec_short = NULL,
    age_tol = 5, nsim = 2000, ptest = c(10000, 28000),
    mode = c("auto", "PT", "T", "P"), sampleorb = "U",
    method = c("angular_momentum", "astrogeo22"),
    s3_uncertainty = c("published", "scaled"),
    seed = NULL, verbose = TRUE, data_dir = NULL){
  # Internal helper functions (local to this function)

  # One-call wrapper: spectral peaks + approximate age -> g, s, k, Earth-Moon

  # @title Mean, sd and 95% interval of a sample, NA when nothing is finite
  # @keywords internal
  # @noRd
  .ctest_ci95 <- function(v){
    if (!any(is.finite(v))){
      return(c(mean = NA_real_, sd = NA_real_, `2.5%` = NA_real_, `97.5%` = NA_real_))
    }
    c(mean = mean(v), sd = stats::sd(v), stats::quantile(v, c(0.025, 0.975)))
  }

  # @title Frequency, period and implied-cycle tables from accepted samples
  # @description Shared by astro_from_peaks() and astro_from_periods().
  # @keywords internal
  # @noRd
  .ctest_tables <- function(samples){
    per <- samples[, c("g2", "g3", "g4", "g5", "s3", "k")]
    sgn <- c(1, 1, 1, 1, -1, 1)          # s3 is a negative frequency
    nper <- length(per)
    freqs <- matrix(NA_real_, nper, 4)
    pers <- matrix(NA_real_, nper, 4)
    for (i in seq_len(nper)){
      fi <- .ctest_ci95(sgn[i] * 360 * 3600 / per[[i]])
      pi_ <- .ctest_ci95(per[[i]])
      freqs[i, ] <- fi
      pers[i, ] <- pi_
      if (i == 1){
        colnames(freqs) <- names(fi)
        colnames(pers) <- names(pi_)
      }
    }
    rownames(pers) <- names(per)
    frequencies <- data.frame(term = names(per),
        freq_arcsec_yr = freqs[, "mean"], freq_sd = freqs[, "sd"],
        freq_lo95 = freqs[, 3], freq_hi95 = freqs[, 4],
        period_yr = pers[, "mean"], period_sd_yr = pers[, "sd"],
        row.names = NULL)
    frequencies$term[frequencies$term == "k"] <- "k (precession constant)"

    cols <- c("ecc405", "ecc130", "ecc124", "ecc99", "ecc95",
        "obl_ks3", "prec_kg5", "prec_kg2", "prec_kg3", "prec_kg4")
    use <- logical(length(cols))
    for (j in seq_along(cols)){
      use[j] <- any(is.finite(samples[[cols[j]]]))
    }
    cols <- cols[use]
    periods <- matrix(NA_real_, length(cols), 4)
    for (j in seq_along(cols)){
      cj <- .ctest_ci95(samples[[cols[j]]])
      periods[j, ] <- cj
      if (j == 1){
        colnames(periods) <- names(cj)
      }
    }
    periods <- data.frame(cycle = cols, period_kyr = periods[, "mean"],
        sd_kyr = periods[, "sd"], lo95 = periods[, 3],
        hi95 = periods[, 4], row.names = NULL)
    list(frequencies = frequencies, periods = periods)
  }

  # Main code
  mode <- match.arg(mode)
  method <- match.arg(method)
  s3_uncertainty <- match.arg(s3_uncertainty)
  if (length(age) == 1){
    age <- c(age - age_tol, age + age_tol)
  }
  has_obl <- !is.null(obl)
  has_prec <- !is.null(prec_long) && !is.null(prec_mid) && !is.null(prec_short)
  if (mode == "auto"){
    if (has_obl && has_prec){
      mode <- "PT"
    } else if (has_obl){
      mode <- "T"
    } else if (has_prec){
      mode <- "P"
    } else {
      stop("Give an obliquity range, the three precession ranges, or both")
    }
  }
  na2 <- c(NA, NA)
  obl_in <- na2
  if (has_obl){
    obl_in <- obl
  }
  p1_in <- na2
  p2_in <- na2
  p3_in <- na2
  if (has_prec){
    p1_in <- prec_long
    p2_in <- prec_mid
    p3_in <- prec_short
  }
  input <- ctest_input(age = age, e405 = e405, e125 = e125, e95 = e95,
      obl = obl_in, p1 = p1_in, p2 = p2_in, p3 = p3_in)
  fit <- ctest(input, nsim = nsim, mode = mode, ptest = ptest,
      sampleorb = sampleorb, seed = seed, verbose = verbose,
      s3_uncertainty = s3_uncertainty)
  mid_age <- mean(age)
  em <- ctest_earthmoon(fit, age_ma = mid_age, method = method, data_dir = data_dir)

  # secular frequencies and precession constant, as periods and frequencies
  tabs <- .ctest_tables(fit$samples)

  out <- list(frequencies = tabs$frequencies, periods = tabs$periods,
      earthmoon = em$summary, sedimentation_rate_m_kyr = .ctest_ci95(fit$samples$SR),
      reference = astrogeo22(mid_age, data_dir = data_dir),
      ctest = fit, em = em, age_ma = age, mode = mode, method = method)
  class(out) <- "astro_from_peaks"
  out
}

#' @title Print an astro_from_peaks result
#' @param x Result to print: an object of class \code{astro_from_peaks}, as returned by
#'   \code{\link{astro_from_peaks}} or \code{\link{astro_from_periods}}. Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{print()}.
#' @return x, invisibly
#' @export
print.astro_from_peaks <- function(x, ...){
  cat(sprintf("\nAstronomical parameters from spectral peaks (CTest %s, %d samples, age %g-%g Ma)\n",
      x$mode, x$ctest$nsim, x$age_ma[1], x$age_ma[2]))
  if (grepl("time", x$mode)){
    cat(sprintf("Time-scale stretch vs the tuning target: %.4f +- %.4f\n",
        x$sedimentation_rate_m_kyr["mean"], x$sedimentation_rate_m_kyr["sd"]))
  } else {
    cat(sprintf("Sedimentation rate: %.3f +- %.3f m/kyr\n",
        x$sedimentation_rate_m_kyr["mean"], x$sedimentation_rate_m_kyr["sd"]))
  }
  if (!is.null(x$ctest$constraints)){
    cat(sprintf("Cycles used: %s\n", paste(names(x$ctest$constraints), collapse = ", ")))
  }
  cat("\n")
  cat("Secular frequencies and precession constant\n")
  print(data.frame(term = x$frequencies$term,
      `arcsec/yr` = sprintf("%8.3f +- %.3f", x$frequencies$freq_arcsec_yr,
          x$frequencies$freq_sd),
      `period (yr)` = sprintf("%9.0f +- %.0f", x$frequencies$period_yr,
          x$frequencies$period_sd_yr),
      check.names = FALSE), row.names = FALSE)
  cat("\nAstronomical periods implied (kyr)\n")
  print(data.frame(cycle = x$periods$cycle,
      kyr = sprintf("%6.2f +- %.2f", x$periods$period_kyr, x$periods$sd_kyr),
      check.names = FALSE), row.names = FALSE)
  cat(sprintf("\nEarth-Moon system (%s inversion)\n", x$method))
  print(round(x$earthmoon, 3))
  cat("\nAstroGeo22 model at the same age, for comparison\n")
  r <- x$reference
  cat(sprintf("  precession %.2f kyr | distance %.2f R_E (%.0f km) | solar day %.2f h | obliquity %.2f deg\n\n",
      r$p_kyr, r$a_RE, r$a_km, r$lod_solar_h, r$obliquity_deg))
  invisible(x)
}
