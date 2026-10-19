#' Lunar mean motion and months per year from the Earth-Moon distance
#' (Kepler's third law, n^2 a^3 = G(M_E + M_M), as in Farhat et al. 2022 Table D.2)
#' @param a_RE Earth-Moon distance in Earth radii.
#' @param kind Type of month, `"sidereal"` or `"synodic"`.
#' @export
months_per_year <- function(a_RE, kind = c("sidereal", "synodic")){
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
  kind <- match.arg(kind)
  n <- sqrt(.AC$GMEM / (a_RE * .AC$RE)^3)        # rad/s
  sid <- n * .AC$yr / (2 * pi)
  if (kind == "sidereal"){
    sid
  } else {
    sid - 1                                      # synodic = sidereal - 1 per year
  }
}
