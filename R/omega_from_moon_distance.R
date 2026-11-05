#' Earth's spin rate implied by an Earth-Moon distance, assuming the
#' Earth-Moon angular momentum is conserved (solar tidal drain neglected).
#' @param a_RE Earth-Moon distance (semi-major axis of the lunar orbit) in Earth radii (numeric vector, > 0; the
#' present-day value is about 60.14). Required (no default).
#' @param L_total Total angular momentum of the Earth-Moon system in kg m^2/s (numeric, > 0), assumed to be
#' conserved. \code{NULL} uses the present-day value (computed from today's Earth-Moon distance and length of day).
#' \code{Default=NULL}.
#' @export
omega_from_moon_distance <- function(a_RE, L_total = NULL){
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
  if (is.null(L_total)){
    L_total <- earthmoon_angular_momentum(.AC$a0_RE, omega0())
  }
  beta <- .AC$ME * .AC$MM / (.AC$ME + .AC$MM)
  out <- numeric(length(a_RE))
  for (i in seq_along(a_RE)){
    a <- a_RE[[i]]
    if (!is.finite(a) || a <= 0){
      out[i] <- NA_real_
      next
    }
    Lspin <- L_total - beta * sqrt(.AC$GMEM * a * .AC$RE)
    if (Lspin <= 0){
      out[i] <- NA_real_
      next
    }
    out[i] <- stats::uniroot(function(w){
      earth_moment_of_inertia(w) * w - Lspin
    }, c(1e-6, 1e-3), tol = 1e-14)$root
  }
  names(out) <- names(a_RE)
  out
}
