## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

#' Full Earth-Moon state from a precession constant
#' @param k_kyr Precession constant expressed as the period of the axial precession, in kyr (numeric, scalar or
#' vector, > 0; about 25.7 today), i.e. a precession rate of \code{360 * 3600 / (k_kyr * 1000)} arcsec/yr. The
#' Earth-Moon distance giving this rate is searched between 20 Earth radii and the present distance (conserving
#' the Earth-Moon angular momentum); values without a solution give \code{NA}. Required (no default).
#' @param obl_deg Obliquity of the Earth used in the precession equation (Farhat et al. 2022, eq. A.7), in degrees
#' (numeric). \code{Default=AC$obl0_deg}, i.e. the present-day obliquity, 23.2545 degrees.
#' @export
em_state <- function(k_kyr, obl_deg = AC$obl0_deg){
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

  freq_arcsec <- function(period_kyr){
    ARCSEC_PER_CYCLE / (period_kyr * 1000)
  }

  ## polar moment of inertia at spin rate omega (their eq. A.3)
  earth_C <- function(omega){
    C0 <- AC$inertia_factor * AC$ME * AC$RE^2
    C0 + 2 * AC$k2f * AC$RE^5 / (9 * AC$G) * (omega^2 - omega0()^2)
  }

  ## total Earth-Moon angular momentum: Earth spin + lunar orbit
  earthmoon_L <- function(a_RE, omega){
    beta <- AC$ME * AC$MM / (AC$ME + AC$MM)
    earth_C(omega) * omega + beta * sqrt(AC$GMEM * a_RE * AC$RE)
  }

  ## spin rate implied by an Earth-Moon distance, angular momentum conserved
  omega_from_a <- function(a_RE, L_total = L_TODAY){
    beta <- AC$ME * AC$MM / (AC$ME + AC$MM)
    out <- numeric(length(a_RE))
    for (i in seq_along(a_RE)){
      a <- a_RE[[i]]
      if (!is.finite(a) || a <= 0){
        out[i] <- NA_real_
        next
      }
      Lspin <- L_total - beta * sqrt(AC$GMEM * a * AC$RE)
      if (Lspin <= 0){
        out[i] <- NA_real_
        next
      }
      out[i] <- uniroot(function(w){
        earth_C(w) * w - Lspin
      }, c(1e-6, 1e-3), tol = 1e-14)$root
    }
    names(out) <- names(a_RE)
    out
  }

  ## luni-solar precession frequency (their eq. A.7), arcsec/yr
  precession_frequency_omega <- function(a_RE, omega = NULL, obl_deg = AC$obl0_deg,
      Ed = AC$Ed0){
    if (is.null(omega)){
      omega <- omega_from_a(a_RE)
    }
    aM <- a_RE * AC$RE
    p_rad_s <- 1.5 * (AC$GMS / AC$aE^3 + AC$GMM / aM^3) * Ed *
        omega / omega0()^2 * cos(obl_deg * pi / 180)
    p_rad_s * AC$arcsec_per_rad * AC$yr
  }

  .AC <- astro_constants()

  AC <- local({
    AU  <- 1.495978707e11                 # m
    day <- 86400                          # s
    G   <- 6.67430e-11                    # m^3 kg^-1 s^-2
    yr  <- 365.25 * day                   # Julian year, s
    GMS  <- 2.959122082853813556248e-4 * AU^3 / day^2      # G * M_sun
    GMEM <- 8.997011395221144381906e-10 * AU^3 / day^2     # G * (M_E + M_M)
    massratio <- 81.30056789872074318737                   # M_E / M_M
    GMM <- GMEM / (1 + massratio)
    GME <- GMEM - GMM
    list(AU = AU, day = day, G = G, yr = yr, GMS = GMS, GMEM = GMEM,
        GMM = GMM, GME = GME, ME = GME / G, MM = GMM / G,
        RE = 6378.1366e3,                       # Earth radius, m
        aE = AU,                                # Earth semi-major axis, m
        inertia_factor = 0.3306947357075918999972,  # C0 / (M_E R_E^2)
        k2f = 0.93,                             # fluid Love number
        Ed0 = 0.003243,                         # dynamical ellipticity, today
        a0_RE = 60.142611,                      # Earth-Moon distance today, R_E
        lod0_sid_h = 23.934468,                 # sidereal day today, h
        obl0_deg = 23.2545,                     # obliquity today, degrees
        p0_arcsec = 50.467718,                  # precession frequency today
        arcsec_per_rad = 180 * 3600 / pi)
  })

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  L_TODAY <- earthmoon_L(AC$a0_RE, omega0())

  # Main code
  p_target <- freq_arcsec(k_kyr)
  out <- vector("list", length(p_target))
  for (j in seq_along(p_target)){
    pt <- p_target[[j]]
    f <- function(a){
      precession_frequency_omega(a, obl_deg = obl_deg) - pt
    }
    a <- tryCatch(uniroot(f, c(20, AC$a0_RE), tol = 1e-10)$root, error = function(e){
      NA_real_
    })
    if (!is.finite(a)){
      out[[j]] <- rep(NA_real_, 7)
      next
    }
    om <- omega_from_a(a)
    lod <- lod_hours(om, "solar")
    out[[j]] <- c(a, a * AC$RE / 1000, lod, lod_hours(om, "sidereal"),
        days_per_year(lod), months_per_year(a, "sidereal"),
        months_per_year(a, "synodic"))
  }
  names(out) <- names(p_target)
  m <- do.call(rbind, out)
  data.frame(k_kyr = k_kyr, p_arcsec = p_target, a_RE = m[, 1], a_km = m[, 2],
      lod_solar_h = m[, 3], lod_sidereal_h = m[, 4],
      days_per_year = m[, 5], sidereal_months_per_year = m[, 6],
      synodic_months_per_year = m[, 7], obliquity_deg = obl_deg)
}
