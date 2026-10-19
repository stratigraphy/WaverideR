#' Luni-solar precession frequency (Farhat et al. 2022, eq. A.7)
#'
#' @param a_RE Earth-Moon distance in Earth radii
#' @param lod_h length of day; if NULL it is derived from `a_RE` by
#'   angular-momentum conservation
#' @param obliquity_deg Earth's obliquity
#' @param Ed dynamical ellipticity at the present spin rate
#' @param lod_kind Type of day, `"solar"` or `"sidereal"`.
#' @return list: p in arcsec/yr and the precession period in kyr
#' @export
precession_frequency <- function(a_RE, lod_h = NULL, obliquity_deg = .AC$obliquity0_deg,
                                 Ed = .AC$Ed0, lod_kind = "solar"){
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

  omega0 <- function(){
    earth_omega(.AC$lod0_sid_h, "sidereal")
  }

  .AC <- astro_constants()

  # ---- main code ----
  if (is.null(lod_h)){
    omega <- omega_from_moon_distance(a_RE)
  } else {
    omega <- earth_omega(lod_h, lod_kind)
  }
  aM <- a_RE * .AC$RE
  p_rad_s <- 1.5 * (.AC$GMS / .AC$aE^3 + .AC$GMM / aM^3) * Ed *
    omega / omega0()^2 * cos(obliquity_deg * pi / 180)
  p_arcsec <- p_rad_s * .AC$arcsec_per_rad * .AC$yr
  list(p_arcsec = p_arcsec, p_kyr = 360 * 3600 / p_arcsec / 1000,
       lod_h = lod_hours(omega, lod_kind), omega = omega)
}
