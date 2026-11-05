# One-stop comparison table

#' Compare a measured precession constant with published models
#'
#' @param age_ma Age of the record in Ma (single numeric value, 0 to 4429), at which the published models
#' (AstroGeo22 and the La2004 extrapolation) are evaluated. Required, no default.
#' @param k_yr Measured precession constant of the record, given as its period in years (numeric, e.g. about 25700
#' for today); when given, rows are added with the Earth-Moon state inverted from it (by angular momentum and on the
#' AstroGeo22 curve, see \code{\link{earthmoon_from_precession}}). \code{NULL} gives only the published models.
#' \code{Default=NULL}.
#' @param obliquity_deg Obliquity in degrees (numeric) used for the angular-momentum inversion of \code{k_yr}; only
#' used when \code{k_yr} is given. \code{NULL} uses the AstroGeo22 obliquity at \code{age_ma}.
#' \code{Default=NULL}.
#' @param data_dir Folder that contains the astronomical data tables \code{AstroGeo22.txt} and
#' \code{La04gisi.csv} (character path); \code{NULL} uses the tables installed with the package (or the folder
#' set with \code{options(cocoR.data = ...)}). \code{Default=NULL}.
#' @return data.frame with one row per source: precession period and frequency,
#'   Earth-Moon distance, LOD, obliquity, plus the k+g and k+s periods
#' @export
compare_published <- function(age_ma, k_yr = NULL, obliquity_deg = NULL,
    data_dir = NULL){
  # Internal helper functions (local to this function)

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

  # Main code
  ag <- astrogeo22(age_ma, data_dir = data_dir)
  la <- la2004_extrap(age_ma)
  row <- function(src, p_kyr, a_RE, lod_h, obl){
    op <- orbital_periods(k_yr = p_kyr * 1000, data_dir = data_dir)
    if (is.na(lod_h)){
      lod_sid <- NA_real_
      dpy <- NA_real_
    } else {
      lod_sid <- lod_hours(earth_omega(lod_h, "solar"), "sidereal")
      dpy <- days_per_year(lod_h)
    }
    if (is.na(a_RE)){
      smpy <- NA_real_
    } else {
      smpy <- months_per_year(a_RE, "sidereal")
    }
    data.frame(source = src, age_ma = age_ma,
        p_kyr = p_kyr, p_arcsec = 360 * 3600 / (p_kyr * 1000),
        a_RE = a_RE, a_km = a_RE * .AC$RE / 1000,
        lod_solar_h = lod_h,
        lod_sidereal_h = lod_sid,
        days_per_year = dpy,
        sidereal_months_per_year = smpy,
        obliquity_deg = obl,
        `k+s3_kyr` = op$period_kyr[op$term == "k+s3"],
        `k+g5_kyr` = op$period_kyr[op$term == "k+g5"],
        `k+g2_kyr` = op$period_kyr[op$term == "k+g2"],
        `k+g4_kyr` = op$period_kyr[op$term == "k+g4"],
        check.names = FALSE)
  }
  out <- rbind(
      row("AstroGeo22 (Farhat et al. 2022)", ag$p_kyr, ag$a_RE, ag$lod_h, ag$obliquity_deg),
      row("La2004 extrapolation (eqs 39-40)", la$p_kyr, la$a_RE, la$lod_h, NA))
  if (!is.null(k_yr)){
    if (is.null(obliquity_deg)){
      obliquity_deg <- ag$obliquity_deg
    }
    inv_am <- earthmoon_from_precession(k_yr / 1000, method = "angular_momentum",
        obliquity_deg = obliquity_deg)
    inv_ag <- earthmoon_from_precession(k_yr / 1000, method = "astrogeo22",
        data_dir = data_dir)
    row_am <- row("measured k (angular momentum)", k_yr / 1000, inv_am$a_RE, inv_am$lod_h, obliquity_deg)
    row_ag <- NULL
    if (!is.null(inv_ag) && nrow(inv_ag)){
      row_ag <- row(sprintf("measured k (on AstroGeo22 curve, model age %s Ma)",
          paste(round(inv_ag$age_ma), collapse = "/")),
          k_yr / 1000, inv_ag$a_RE[1], inv_ag$lod_h[1], inv_ag$obliquity_deg[1])
    }
    out <- rbind(row_am, row_ag, out)
  }
  rownames(out) <- NULL
  out
}
