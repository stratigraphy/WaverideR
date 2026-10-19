## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

## ---- 8. worked example -----------------------------------------------------

#' The Pridoli example: the four terms resolved at Cellon (Carnic Alps)
#' on the astrochronozone time scale of Arts and Da Silva.
#' @param file Name of the example file.
example_pridoli <- function(file = NULL){
  # ---- internal helper functions (local to this function) ----

  period_kyr  <- function(freq_arcsec){
    ARCSEC_PER_CYCLE / abs(freq_arcsec) / 1000
  }

  astrogeo22_nearest <- function(age_ma = NULL, file = NULL, url = ASTROGEO22_URL){
    if (is.null(file)){
      file <- file.path(tempdir(), "AstroGeo22.txt")
      if (!file.exists(file)){
        try(utils::download.file(url, file, quiet = TRUE), silent = TRUE)
      }
    }
    if (!file.exists(file)){
      message("AstroGeo22.txt not available (no local copy and no network); ",
              "comparisons with the tidal model are skipped.")
      return(NULL)
    }
    A <- utils::read.table(file, comment.char = "%")
    out <- data.frame(age_ma = A[[1]] * 1000, a_RE = A[[2]], lod_solar_h = A[[5]],
                      obliquity_deg = A[[8]], p_arcsec = A[[11]])
    out$k_kyr <- period_kyr(out$p_arcsec)
    out$days_per_year <- days_per_year(out$lod_solar_h)
    if (is.null(age_ma)){
      return(out)
    }
    idx <- numeric(length(age_ma))
    for (i in seq_along(age_ma)){
      idx[i] <- which.min(abs(out$age_ma - age_ma[[i]]))
    }
    names(idx) <- names(age_ma)
    out[idx, ]
  }

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  ## ---- 6. reference models ---------------------------------------------------

  ## AstroGeo22 (Farhat et al. 2022).  Supply a local copy of AstroGeo22.txt, or
  ## let the function fetch it from the Acycle repository.  Columns used:
  ## 1 time (Ga), 2 a (R_E), 5 length of day (h), 8 obliquity (deg),
  ## 11 precession frequency (arcsec/yr).
  ASTROGEO22_URL <- paste0("https://raw.githubusercontent.com/mingsongli/acycle/",
                           "master/code/package/secular_periods/AstroGeo22.txt")

  ## Zeeden et al. (2023), their table 1, combined Devonian estimate (1 sigma)
  ZEEDEN2023 <- list(k_kyr = c(mean = 20.70, sd = 0.22),
                     periods = c("k+g5" = 19.35, "k+g2" = 18.47,
                                 "k+g4" = 16.13, "k+g3" = 16.18, "k+s3" = 29.50),
                     sd = c("k+g5" = 0.20, "k+g2" = 0.18, "k+g4" = 0.23,
                            "k+g3" = 0.16, "k+s3" = 0.46))

  # ---- main code ----
  windows <- list("k+s3" = c(29.0, 30.4),    # obliquity
                  "k+g5" = c(19.5, 20.6),
                  "k+g2" = c(18.7, 19.5),
                  "k+g4" = c(16.2, 17.1))    # k+g4 / k+g3 unresolved
  peaks <- c(29.68, 20.02, 19.11, 16.94, 16.48)   # measured MTM peaks, kyr
  ## the consistency test is the estimator used in the paper
  fit <- astro_ctest(windows, age_ma = c(418, 424), nsim = 3000, seed = 41)
  cat("\nPrecession constant k (kyr), consistency test\n")
  print(round(fit$k, 3))
  ## the simple averaging estimator, for comparison
  fit2 <- astro_invert(windows, nsim = 5000, obl_deg = fit$obl_deg, seed = 41)
  cat("\nPrecession constant k (kyr), simple average over terms\n")
  print(round(fit2$k, 3))
  cat("\nEarth-Moon state\n")
  print(round(fit$state, 3))
  ag <- astrogeo22_nearest(421)
  if (!is.null(ag)){
    cat("\nAstroGeo22 at 421 Ma (Farhat et al. 2022)\n")
    print(round(ag[, c("k_kyr", "a_RE", "lod_solar_h", "days_per_year",
                       "obliquity_deg")], 3))
    cat("\nnine principal terms at that k (kyr), as calculate_orbit9 returns them\n")
    print(orbit9(ag$k_kyr)[, c("term", "period_kyr")], row.names = FALSE, digits = 6)
  }
  cat("\nZeeden et al. (2023) combined Devonian: k =",
      ZEEDEN2023$k_kyr[["mean"]], "+/-", ZEEDEN2023$k_kyr[["sd"]], "kyr\n")
  plot_astro(fit, peaks = peaks, age_ma = 421, file = file)
  invisible(fit)
}
