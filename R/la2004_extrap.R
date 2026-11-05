# Older published models for comparison

#' La2004 polynomial extrapolation (Laskar et al. 2004, eqs. 39-40)
#' Strictly valid only over a few hundred Myr; shown for comparison.
#' @param age_ma Age or ages in Ma (numeric vector, >= 0; positive into the past); one output row per age. The
#' polynomials are only reliable over the last few hundred Myr. Required (no default).
#' @export
la2004_extrap <- function(age_ma){
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
  # In La2004 time is negative into the past, so T = -age.
  Tg <- -age_ma / 1000                     # Gyr (negative = past)
  a <- 60.142611 + 6.100887 * Tg - 2.709407 * Tg^2 + 1.366779 * Tg^3 - 1.484062 * Tg^4
  p <- 50.475838 - 26.368583 * Tg + 21.890862 * Tg^2
  data.frame(age_ma = age_ma, a_RE = a, a_km = a * .AC$RE / 1000,
      p_arcsec = p, p_kyr = 360 * 3600 / p / 1000,
      lod_h = lod_from_moon_distance(a))
}
