#' @title Astronomical parameters from spectral peaks
#'
#' @description Wraps ctest() and ctest_earthmoon(): give the wavelength ranges
#' (m) of the cycles you identified in a depth-domain spectrum and an
#' approximate age, and get back the secular frequencies g2-g5 and s3, the
#' precession constant k, the Earth-Moon distance and the length of day, with
#' uncertainties.
#'
#' @param age approximate age in Ma: c(min, max), or a single value (then
#'   `age_tol` Ma is used on either side)
#' @param e405,e125,e95 c(min, max) wavelength (m) of the 405-kyr, ~125-kyr and
#'   ~95-kyr eccentricity cycles (all three are required)
#' @param obl c(min, max) wavelength (m) of the obliquity cycle, if identified
#' @param prec_long,prec_mid,prec_short c(min, max) wavelength (m) of the
#'   precession cycles, longest to shortest (k+g5; k+g2; k+g4/k+g3). Give what
#'   you have; all three are needed for the precession test.
#' @param age_tol half-width (Ma) used when `age` is a single value
#' @param nsim,ptest,sampleorb,s3_uncertainty,seed,verbose passed to ctest()
#' @param mode "auto" picks PT / T / P from what you supplied
#' @param method how k is converted to an Earth-Moon state:
#'   "angular_momentum" or "astrogeo22"
#' @param data_dir folder with the AstroGeo22 table
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
  # ---- internal helper functions (local to this function) ----

  # ---------------------------------------------------------------------------
  # One-call wrapper: spectral peaks + approximate age -> g, s, k, Earth-Moon
  # ---------------------------------------------------------------------------

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

  # ---- main code ----
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
#' @param x an astro_from_peaks object
#' @param ... ignored
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
