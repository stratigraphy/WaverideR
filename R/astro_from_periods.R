#' @title Astronomical parameters from a tuned record (cycles already in kyr)
#'
#' @description Same output as astro_from_peaks(), but the input is cycle
#' periods in kyr and the time scale is set by the tuning target rather than by
#' a sampled sedimentation rate. Give only the terms you actually resolved.
#'
#' @param cycles named list of c(min, max) periods in kyr (see astro_mc())
#' @param age approximate age in Ma
#' @param anchor period assumed for the tuning target, kyr (405 by default)
#' @param anchor_term the term you tuned to
#' @param nsim,ptest,s3_uncertainty,seed,verbose passed to astro_mc()
#' @param method "angular_momentum" or "astrogeo22"
#' @param data_dir folder with the AstroGeo22 table
#' @param ... further arguments to astro_mc()
#' @return object of class `astro_from_peaks` (see astro_from_peaks())
#' @export
astro_from_periods <- function(cycles, age, anchor = 405, anchor_term = "e405",
                               nsim = 2000, ptest = c(10000, 28000),
                               method = c("angular_momentum", "astrogeo22"),
                               s3_uncertainty = c("published", "scaled"),
                               seed = NULL, verbose = TRUE, data_dir = NULL, ...){
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
  method <- match.arg(method)
  s3_uncertainty <- match.arg(s3_uncertainty)
  fit <- astro_mc(cycles, age = age, domain = "time", anchor = anchor,
                  anchor_term = anchor_term, nsim = nsim, ptest = ptest,
                  s3_uncertainty = s3_uncertainty, seed = seed, verbose = verbose, ...)
  mid_age <- mean(fit$age_ma)
  em <- ctest_earthmoon(fit, age_ma = mid_age, method = method, data_dir = data_dir)
  tabs <- .ctest_tables(fit$samples)
  out <- list(frequencies = tabs$frequencies, periods = tabs$periods, earthmoon = em$summary,
              sedimentation_rate_m_kyr = .ctest_ci95(fit$samples$scale),
              reference = astrogeo22(mid_age, data_dir = data_dir),
              ctest = fit, em = em, age_ma = fit$age_ma,
              mode = fit$mode, method = method)
  class(out) <- "astro_from_peaks"
  out
}
