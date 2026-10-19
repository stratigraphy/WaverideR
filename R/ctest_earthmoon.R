#' @title Earth-Moon distance, length of day and obliquity for CTest samples
#'
#' @description Each sampled precession constant is converted with
#' precession_frequency() (Farhat et al. 2022, eq. A.7), either by inverting
#' the AstroGeo22 curve or by assuming Earth-Moon angular-momentum conservation.
#'
#' @param x a ctest object, or a numeric vector of precession constants (yr)
#' @param age_ma age of the record; used for the obliquity of the
#'   angular-momentum inversion and for the comparison row
#' @param method "angular_momentum" (default) or "astrogeo22"
#' @param obliquity_deg obliquity used in the angular-momentum inversion;
#'   NULL takes it from AstroGeo22 at `age_ma`
#' @param probs quantiles to report
#' @param data_dir folder with the AstroGeo22 table
#' @return list with `summary`, `samples`, `method`, `age_ma` and, when an age
#'   is known, `model_at_age`
#' @export
ctest_earthmoon <- function(x, age_ma = NULL, method = c("angular_momentum", "astrogeo22"),
                            obliquity_deg = NULL, probs = c(0.025, 0.5, 0.975),
                            data_dir = NULL){
  # ---- internal helper functions (local to this function) ----

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

  # --- Astronomical periods at a given age (calculate_orbit9.m) --------------
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

  # ---- main code ----
  method <- match.arg(method)
  if (is.numeric(x)){
    k <- x
  } else {
    k <- x$samples$k
  }
  if (is.null(age_ma) && !is.numeric(x)){
    age_ma <- mean(x$age_ma)
  }
  if (is.null(obliquity_deg) && !is.null(age_ma)){
    obliquity_deg <- astrogeo22(age_ma, data_dir = data_dir)$obliquity_deg
  }
  if (is.null(obliquity_deg)){
    obliquity_deg <- .AC$obliquity0_deg
  }
  p_arcsec <- 360 * 3600 / k
  if (method == "angular_momentum"){
    # invert once on a grid, then map the samples (fast and smooth)
    grid_a <- seq(30, 60.2, by = 0.01)
    grid_p <- numeric(length(grid_a))
    for (j in seq_along(grid_a)){
      grid_p[j] <- precession_frequency(grid_a[j], NULL, obliquity_deg)$p_arcsec
    }
    a_RE <- stats::approx(grid_p, grid_a, xout = p_arcsec, rule = 1)$y
    lod <- lod_from_moon_distance(a_RE)
    obl <- rep(obliquity_deg, length(a_RE))
    model_age <- rep(NA_real_, length(a_RE))
  } else {
    tab <- .astrogeo22_table(data_dir)
    keep <- tab$age_ma <= 1000
    a_RE <- stats::approx(tab$pre[keep], tab$a[keep], xout = p_arcsec, rule = 1)$y
    lod <- stats::approx(tab$pre[keep], tab$lod[keep], xout = p_arcsec, rule = 1)$y
    obl <- stats::approx(tab$pre[keep], tab$obli[keep], xout = p_arcsec, rule = 1)$y
    model_age <- stats::approx(tab$pre[keep], tab$age_ma[keep], xout = p_arcsec, rule = 1)$y
  }
  q <- function(v){
    c(mean = mean(v, na.rm = TRUE), stats::quantile(v, probs, na.rm = TRUE))
  }
  summ <- rbind(
    `k, precession constant (kyr)` = q(k / 1000),
    `precession frequency (arcsec/yr)` = q(p_arcsec),
    `Earth-Moon distance (R_E)` = q(a_RE),
    `Earth-Moon distance (km)` = q(a_RE * .AC$RE / 1000),
    `solar day (h)` = q(lod),
    `sidereal day (h)` = q(lod_hours(earth_omega(lod, "solar"), "sidereal")),
    `days per year` = q(days_per_year(lod)),
    `sidereal months per year` = q(months_per_year(a_RE, "sidereal")),
    `synodic months per year` = q(months_per_year(a_RE, "synodic")),
    `obliquity (deg)` = q(obl))
  if (method == "astrogeo22"){
    summ <- rbind(summ, `equivalent AstroGeo22 age (Ma)` = q(model_age))
  }
  out <- list(summary = as.data.frame(summ), method = method, age_ma = age_ma,
              samples = data.frame(k = k, p_arcsec = p_arcsec, a_RE = a_RE,
                                   a_km = a_RE * .AC$RE / 1000,
                                   lod_solar_h = lod,
                                   lod_sidereal_h = lod_hours(earth_omega(lod, "solar"), "sidereal"),
                                   days_per_year = days_per_year(lod),
                                   sidereal_months_per_year = months_per_year(a_RE, "sidereal"),
                                   synodic_months_per_year = months_per_year(a_RE, "synodic"),
                                   obliquity_deg = obl, model_age_ma = model_age))
  if (!is.null(age_ma)){
    out$model_at_age <- astrogeo22(age_ma, data_dir = data_dir)
  }
  out
}
