#' Earth's polar moment of inertia at spin rate omega (Farhat et al. 2022, A.3)
#' @param omega Angular velocity (spin rate) of the Earth, in rad/s (numeric, scalar or vector, > 0), e.g. from
#' \code{\link{earth_omega}}. Present-day value is about 7.292e-5 rad/s. Required (no default).
#' @export
earth_moment_of_inertia <- function(omega){
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

  omega0 <- function(){
    earth_omega(.AC$lod0_sid_h, "sidereal")
  }

  .AC <- astro_constants()

  # Main code
  C0 <- .AC$inertia_factor * .AC$ME * .AC$RE^2
  C0 + 2 * .AC$k2f * .AC$RE^5 / (9 * .AC$G) * (omega^2 - omega0()^2)
}
