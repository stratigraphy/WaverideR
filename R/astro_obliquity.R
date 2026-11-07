## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

#' Obliquity at a given age, from AstroGeo22 if available, else today's value
#' @param age_ma Age in Ma at which the obliquity is wanted (numeric, single value; if a vector is given, only its
#'   first value is used). The nearest row of the AstroGeo22 table installed with WaverideR is taken (if that copy
#'   is missing it is downloaded to \code{tempdir()}); if no table is available today's obliquity (23.2545 degrees)
#'   is returned.
#'   Required, no default.
#' @export
astro_obliquity <- function(age_ma){
  # Internal helper functions (local to this function)

  period_kyr  <- function(freq_arcsec){
    ARCSEC_PER_CYCLE / abs(freq_arcsec) / 1000
  }

  astrogeo22_nearest <- function(age_ma = NULL, file = NULL, url = ASTROGEO22_URL){
    if (is.null(file)){
      # the copy installed with WaverideR; downloaded only if that is missing
      file <- system.file("extdata", "AstroGeo22.txt", package = "WaverideR")
      if (!nzchar(file)){
        file <- file.path(tempdir(), "AstroGeo22.txt")
        if (!file.exists(file)){
          try(utils::download.file(url, file, quiet = TRUE), silent = TRUE)
        }
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

  AC <- local({
    AU  <- 1.495978707e11                 # m
    day <- 86400                          # s
    G   <- 6.67430e-11                    # m^3 kg^-1 s^-2
    yr  <- 365.25 * day                   # Julian year, s
    GMS  <- 2.959122082853813556248e-4 * AU^3 / day^2      # G * M_sun
    GMEM <- 8.997011395221144381906e-10 * AU^3 / day^2     # G * (M_E + M_M)
    massratio <- 81.30056789872074318737                   # M_E / M_M
    GMM <- GMEM / (1 + massratio)
    GME <- GMEM - GMM
    list(AU = AU, day = day, G = G, yr = yr, GMS = GMS, GMEM = GMEM,
        GMM = GMM, GME = GME, ME = GME / G, MM = GMM / G,
        RE = 6378.1366e3,                       # Earth radius, m
        aE = AU,                                # Earth semi-major axis, m
        inertia_factor = 0.3306947357075918999972,  # C0 / (M_E R_E^2)
        k2f = 0.93,                             # fluid Love number
        Ed0 = 0.003243,                         # dynamical ellipticity, today
        a0_RE = 60.142611,                      # Earth-Moon distance today, R_E
        lod0_sid_h = 23.934468,                 # sidereal day today, h
        obl0_deg = 23.2545,                     # obliquity today, degrees
        p0_arcsec = 50.467718,                  # precession frequency today
        arcsec_per_rad = 180 * 3600 / pi)
  })

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  # Reference models

  ## AstroGeo22 (Farhat et al. 2022).  Supply a local copy of AstroGeo22.txt, or
  ## let the function fetch it from the Acycle repository.  Columns used:
  ## 1 time (Ga), 2 a (R_E), 5 length of day (h), 8 obliquity (deg),
  ## 11 precession frequency (arcsec/yr).
  ASTROGEO22_URL <- paste0("https://raw.githubusercontent.com/mingsongli/acycle/",
      "master/code/package/secular_periods/AstroGeo22.txt")

  # Main code
  ag <- astrogeo22_nearest(age_ma)
  if (is.null(ag)){
    AC$obl0_deg
  } else {
    ag$obliquity_deg[1]
  }
}
