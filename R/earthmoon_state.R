#' Full Earth-Moon state from any one observable
#'
#' Give exactly one of: `p_arcsec`/`p_kyr`/`k_yr` (precession), `a_RE`/`a_km`
#' (Earth-Moon distance), `lod_h` (solar day), `months_per_year`, or an
#' observed `cycle_kyr` together with its `term` (e.g. the obliquity cycle).
#'
#' @param method how the Earth's spin and the lunar distance are linked when
#'   inverting a precession constant: "angular_momentum" (Earth-Moon angular
#'   momentum conserved) or "astrogeo22" (read off the Farhat et al. 2022 curve)
#' @param obliquity_deg obliquity to use; defaults to the AstroGeo22 value at
#'   `age_ma`, else the present value
#' @param p_arcsec Precession constant in arcsec/yr.
#' @param p_kyr Precession period in kyr.
#' @param k_yr Precession constant in arcsec/yr.
#' @param a_RE Earth-Moon distance in Earth radii.
#' @param a_km Earth-Moon distance in km.
#' @param lod_h Length of day in hours.
#' @param months_per_year Number of lunar months per year.
#' @param months_kind Type of month, `"sidereal"` or `"synodic"`.
#' @param cycle_kyr Period of a cycle in kyr.
#' @param term Frequency term, e.g. `"k+s3"`.
#' @param age_ma Age(s) in Ma.
#' @param data_dir Directory with the data files (`NULL` = the data installed with the package).
#' @return one-row-per-input data.frame with the precession constant, the
#'   Earth-Moon distance, the solar and sidereal day, days per year, sidereal
#'   and synodic months per year and the obliquity
#' @export
earthmoon_state <- function(p_arcsec = NULL, p_kyr = NULL, k_yr = NULL,
                            a_RE = NULL, a_km = NULL, lod_h = NULL,
                            months_per_year = NULL, months_kind = "sidereal",
                            cycle_kyr = NULL, term = "k+s3",
                            method = c("angular_momentum", "astrogeo22"),
                            obliquity_deg = NULL, age_ma = NULL, data_dir = NULL){
  # ---- internal helper functions (local to this function) ----

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

  .AC <- astro_constants()

  # ---- main code ----
  method <- match.arg(method)
  # vectorised over whichever observable was supplied
  n <- max(lengths(list(p_arcsec, p_kyr, k_yr, a_RE, a_km, lod_h, months_per_year, cycle_kyr)), 1)
  if (n > 1){
    pick <- function(v, i){
      if (is.null(v)){
        NULL
      } else {
        v[[min(i, length(v))]]
      }
    }
    rows <- vector("list", n)
    for (i in seq_len(n)){
      rows[[i]] <- earthmoon_state(pick(p_arcsec, i), pick(p_kyr, i), pick(k_yr, i), pick(a_RE, i),
                                   pick(a_km, i), pick(lod_h, i), pick(months_per_year, i), months_kind,
                                   pick(cycle_kyr, i), term, method,
                                   pick(obliquity_deg, i), pick(age_ma, i), data_dir)
    }
    return(do.call(rbind, rows))
  }
  if (is.null(obliquity_deg)){
    if (!is.null(age_ma)){
      obliquity_deg <- astrogeo22(age_ma, data_dir = data_dir)$obliquity_deg
    } else {
      obliquity_deg <- .AC$obliquity0_deg
    }
  }
  if (!is.null(cycle_kyr)){
    k_yr <- cycles_to_k(cycle_kyr, term, data_dir = data_dir)$k_yr
  }
  if (!is.null(k_yr)){
    p_arcsec <- 360 * 3600 / k_yr
  }
  if (!is.null(p_kyr)){
    p_arcsec <- 360 * 3600 / (p_kyr * 1000)
  }
  if (!is.null(a_km)){
    a_RE <- a_km * 1000 / .AC$RE
  }
  if (!is.null(months_per_year)){
    a_RE <- moon_distance_from_months(months_per_year, months_kind)
  }
  if (!is.null(lod_h) && is.null(a_RE)){
    a_RE <- moon_distance_from_lod(lod_h)
  }

  model_age <- NA_real_
  if (is.null(a_RE)){
    if (is.null(p_arcsec)){
      stop("Give one observable (precession, distance, day length or months per year)")
    }
    inv <- earthmoon_from_precession(p_arcsec = p_arcsec, method = method,
                                     obliquity_deg = obliquity_deg, data_dir = data_dir)
    a_RE <- inv$a_RE[1]
    lod_h <- inv$lod_h[1]
    if (method == "astrogeo22"){
      obliquity_deg <- inv$obliquity_deg[1]
      model_age <- inv$age_ma[1]
    }
  }
  if (is.null(lod_h)){
    lod_h <- lod_from_moon_distance(a_RE)
  }
  if (is.null(p_arcsec)){
    p_arcsec <- precession_frequency(a_RE, lod_h, obliquity_deg)$p_arcsec
  }

  data.frame(p_arcsec = p_arcsec,
             k_kyr = 360 * 3600 / p_arcsec / 1000,
             a_RE = a_RE, a_km = a_RE * .AC$RE / 1000,
             lod_solar_h = lod_h,
             lod_sidereal_h = lod_hours(earth_omega(lod_h, "solar"), "sidereal"),
             days_per_year = days_per_year(lod_h),
             sidereal_months_per_year = months_per_year(a_RE, "sidereal"),
             synodic_months_per_year = months_per_year(a_RE, "synodic"),
             obliquity_deg = obliquity_deg,
             model_age_ma = model_age)
}
