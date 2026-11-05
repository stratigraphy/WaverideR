#' Invert the AstroGeo22 precession curve: Earth-Moon state for a measured
#' precession frequency / period.
#'
#' The model curve is not strictly monotonic (tidal resonances), so every
#' crossing is returned.
#'
#' @param p_kyr Period of the axial precession, in kyr (numeric, scalar or vector, > 0; about 25.7 kyr today), i.e.
#' \code{360 * 3600 / (p_kyr * 1000)} arcsec/yr. Ignored when \code{p_arcsec} is given; one of the two is
#' required. \code{Default=NULL}.
#' @param age_range Age range to which the search along the AstroGeo22 curve is restricted, in Ma (numeric vector
#' \code{c(min, max)}, or \code{NULL} for the whole curve). Only used with \code{method="astrogeo22"}.
#' \code{Default=NULL}.
#' @param p_arcsec Precession constant (rate of axial precession), in arcsec/yr (numeric, scalar or vector, > 0;
#' about 50.47 today). Takes precedence over \code{p_kyr}. \code{Default=NULL}.
#' @param method How the Earth-Moon state is found (character): \code{"astrogeo22"} finds every age at which the
#' AstroGeo22 precession curve crosses the requested value (the curve is not monotonic, so several ages may be
#' returned) and returns the AstroGeo22 state there; \code{"angular_momentum"} solves for the Earth-Moon distance
#' (25 to 80 Earth radii) that gives the requested precession with \code{\link{precession_frequency}}, assuming
#' conservation of the Earth-Moon angular momentum. \code{Default=c("astrogeo22", "angular_momentum")}, i.e.
#' \code{"astrogeo22"}.
#' @param obliquity_deg Obliquity of the Earth, in degrees (numeric), used by \code{method="angular_momentum"};
#' \code{NULL} uses the present-day value 23.2545. \code{Default=NULL}.
#' @param data_dir Directory that contains the file \code{AstroGeo22.txt} (character, or \code{NULL}); \code{NULL}
#' uses the data installed with the package (or \code{options(cocoR.data = ...)}). \code{Default=NULL}.
#' @export
earthmoon_from_precession <- function(p_kyr = NULL, p_arcsec = NULL,
    method = c("astrogeo22", "angular_momentum"),
    obliquity_deg = NULL, age_range = NULL,
    data_dir = NULL){
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
  method <- match.arg(method)
  if (is.null(p_arcsec)){
    if (is.null(p_kyr)){
      stop("Give p_kyr or p_arcsec")
    }
    p_arcsec <- 360 * 3600 / (p_kyr * 1000)
  }
  if (method == "angular_momentum"){
    if (is.null(obliquity_deg)){
      obliquity_deg <- .AC$obliquity0_deg
    }
    out <- vector("list", length(p_arcsec))
    for (j in seq_along(p_arcsec)){
      pa <- p_arcsec[[j]]
      ff <- function(a){
        precession_frequency(a, NULL, obliquity_deg)$p_arcsec - pa
      }
      lo <- 25
      hi <- 80    # above 60.14 R_E means "farther than today"

      if (ff(lo) * ff(hi) > 0){
        out[[j]] <- data.frame(a_RE = NA, a_km = NA, lod_h = NA,
            obliquity_deg = obliquity_deg,
            p_arcsec = pa, p_kyr = 360 * 3600 / pa / 1000)
      } else {
        a <- stats::uniroot(ff, c(lo, hi), tol = 1e-10)$root
        out[[j]] <- data.frame(a_RE = a, a_km = a * .AC$RE / 1000, lod_h = lod_from_moon_distance(a),
            obliquity_deg = obliquity_deg, p_arcsec = pa,
            p_kyr = 360 * 3600 / pa / 1000)
      }
    }
    names(out) <- names(p_arcsec)
    return(do.call(rbind, out))
  }
  tab <- .astrogeo22_table(data_dir)
  if (!is.null(age_range)){
    tab <- tab[tab$age_ma >= age_range[1] & tab$age_ma <= age_range[2], ]
  }
  res <- vector("list", length(p_arcsec))
  for (j in seq_along(p_arcsec)){
    pa <- p_arcsec[[j]]
    d <- tab$pre - pa
    cross <- which(d[-1] * d[-length(d)] <= 0)
    if (!length(cross)){
      # outside the model range: fall back to the nearest end of the curve
      i <- which.min(abs(d))
      a <- astrogeo22(tab$age_ma[i], data_dir = data_dir)
      a$p_arcsec_requested <- pa
      a$p_kyr_requested <- 360 * 3600 / pa / 1000
      a$outside_model_range <- TRUE
      res[[j]] <- a
      next
    }
    ages <- numeric(length(cross))
    for (m in seq_along(cross)){
      i <- cross[[m]]
      if (d[i + 1] == d[i]){
        ages[m] <- tab$age_ma[i]
      } else {
        ages[m] <- tab$age_ma[i] - d[i] * (tab$age_ma[i + 1] - tab$age_ma[i]) / (d[i + 1] - d[i])
      }
    }
    names(ages) <- names(cross)
    a <- astrogeo22(ages, data_dir = data_dir)
    a$p_arcsec_requested <- pa
    a$p_kyr_requested <- 360 * 3600 / pa / 1000
    a$outside_model_range <- FALSE
    res[[j]] <- a
  }
  names(res) <- names(p_arcsec)
  do.call(rbind, res)
}
