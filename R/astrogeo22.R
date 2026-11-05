# AstroGeo22 (Farhat et al. 2022 nominal model)

#' Earth-Moon state from the AstroGeo22 model at a given age
#'
#' @param age_ma Age(s) in Ma at which the Earth-Moon state is returned (numeric vector, 0 to 4429; one output row
#'   per value). Ages outside the table give \code{NA} when \code{interpolate=TRUE}. Required, no default.
#' @param interpolate How values are taken from the 1-Myr rows of the AstroGeo22 table (logical): \code{TRUE}
#'   interpolates linearly between rows, \code{FALSE} takes the nearest row (as calculate_orbit9.m of Acycle does).
#'   \code{Default=TRUE}.
#' @param data_dir Folder that contains the file \code{AstroGeo22.txt} (character path); \code{NULL} uses the copy
#'   installed with the package (or the folder set with \code{options(cocoR.data = ...)}). \code{Default=NULL}.
#' @return data.frame: age_ma, a_RE, a_km, lod_h, obliquity_deg, p_arcsec,
#'   p_kyr and the model's min/max envelope for each
#' @export
astrogeo22 <- function(age_ma, interpolate = TRUE, data_dir = NULL){
  # Internal helper functions (local to this function)

  coco_data_dir <- function(){
    if (!is.null(.coco_env$src)){
      src_data <- file.path(dirname(.coco_env$src), "..", "data")
    } else {
      src_data <- NULL
    }
    cand <- c(getOption("cocoR.data"), system.file("extdata", package = "WaverideR"), "data", file.path("cocoR", "data"), src_data)
    for (d in cand){
      if (!is.null(d) && file.exists(file.path(d, "AstroGeo22.txt"))){
        return(d)
      }
    }
    stop("Cannot find AstroGeo22.txt / La04gisi.csv. Reinstall WaverideR or set options(cocoR.data = '/path/to/data').")
  }

  astro_constants <- function(){
    AU <- 1.495978707e11            # m
    day <- 86400                    # s
    G <- 6.67430e-11                # m^3 kg^-1 s^-2
    yr <- 365.25 * day              # Julian year, s
    GMS <- 2.959122082853813556248e-4 * AU^3 / day^2
    GMEM <- 8.997011395221144381906e-10 * AU^3 / day^2   # G(M_E + M_M)
    massratio <- 81.30056789872074318737                 # M_E / M_M
    GMM <- GMEM / (1 + massratio)
    GME <- GMEM - GMM
    RE <- 6378.1366e3               # m
    list(AU = AU, day = day, G = G, yr = yr,
        GMS = GMS, GMEM = GMEM, GMM = GMM, GME = GME,
        ME = GME / G, MM = GMM / G, massratio = massratio,
        RE = RE, aE = 1.495978707e11,          # Earth semi-major axis, m
        inertia_factor = 0.3306947357075918999972,       # C0 / (M_E R_E^2)
        k2f = 0.93,                                      # fluid Love number
        Ed0 = 0.003243,                                  # dynamical ellipticity today
        a0_RE = 60.142611,                               # Moon, present, R_E
        lod0_sid_h = 23.934468,                          # sidereal day, h
        obliquity0_deg = 23.2545,
        p0_arcsec = 50.467718,
        arcsec_per_rad = 180 * 3600 / pi)
  }

  .astrogeo22_table <- function(data_dir = NULL){
    if (is.null(data_dir)){
      data_dir <- coco_data_dir()
    }
    tab <- utils::read.table(file.path(data_dir, "AstroGeo22.txt"), comment.char = "%")
    names(tab)[1:13] <- c("t_ga", "a", "amin", "amax", "lod", "lodmin", "lodmax",
        "obli", "oblimin", "oblimax", "pre", "premin", "premax")
    tab$age_ma <- tab$t_ga * 1000
    tab
  }

  # Astronomical periods at a given age (calculate_orbit9.m)
  # Nine principal periods from the La2004 secular frequencies g_i, s_i
  # (La04gisi.csv) and the precession constant k(t) of AstroGeo22
  # (AstroGeo22.txt, column 11, arcsec/yr), both bundled in Acycle
  # code/package/secular_periods and copied to acycle2R/data.
  #   1 g2-g5 (long ecc.)   2 g3-g2   3 g4-g2   4 g3-g5   5 g4-g5 (short ecc.)
  #   6 k+s3 (obliquity)    7 k+g5    8 k+g2    9 k+g4 (precession)
  # The age is rounded to the nearest integer Ma (AstroGeo22 has 1-Myr steps,
  # 0-4429 Ma). g and s are held at their present La2004 values.
  .coco_env <- new.env()

  .coco_env$src <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e){
    NULL
  })

  .AC <- astro_constants()

  # Main code
  tab <- .astrogeo22_table(data_dir)
  f <- function(col){
    if (interpolate){
      stats::approx(tab$age_ma, tab[[col]], xout = age_ma, rule = 1)$y
    } else {
      idx <- integer(length(age_ma))
      for (i in seq_along(age_ma)){
        idx[i] <- which.min(abs(tab$age_ma - age_ma[[i]]))
      }
      names(idx) <- names(age_ma)
      tab[[col]][idx]
    }
  }
  out <- data.frame(age_ma = age_ma,
      a_RE = f("a"), a_min_RE = f("amin"), a_max_RE = f("amax"),
      lod_h = f("lod"), lod_min_h = f("lodmin"), lod_max_h = f("lodmax"),
      obliquity_deg = f("obli"), obliquity_min_deg = f("oblimin"),
      obliquity_max_deg = f("oblimax"),
      p_arcsec = f("pre"), p_min_arcsec = f("premin"),
      p_max_arcsec = f("premax"))
  out$a_km <- out$a_RE * .AC$RE / 1000
  out$p_kyr <- 360 * 3600 / out$p_arcsec / 1000
  out$lod_solar_h <- out$lod_h
  out$lod_sidereal_h <- lod_hours(earth_omega(out$lod_h, "solar"), "sidereal")
  out$days_per_year <- days_per_year(out$lod_h)
  out$sidereal_months_per_year <- months_per_year(out$a_RE, "sidereal")
  out$synodic_months_per_year <- months_per_year(out$a_RE, "synodic")
  out
}
