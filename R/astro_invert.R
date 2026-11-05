## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

#' Joint inversion of several measured terms, with Monte Carlo uncertainty
#'
#' @param windows Measured period ranges of the k-bearing (obliquity/precession) terms, in kyr (named list;
#'   each element a numeric vector \code{c(low, high)}). Allowed names: \code{"k+s3"} (obliquity), \code{"k+g5"},
#'   \code{"k+g2"}, \code{"k+g4"} and \code{"k+g3"} (precession), e.g. \code{list("k+s3" = c(29.0, 30.4))}.
#'   In each draw one period is sampled uniformly in every window, converted into a precession frequency, and
#'   the frequencies of all windows are averaged. Required, no default.
#' @param nsim Number of Monte Carlo draws (integer, > 0); draws whose k falls outside \code{k_range} are
#'   discarded afterwards. \code{Default=5000}.
#' @param obl_deg Obliquity in degrees used in the precession relation when k is converted into the Earth-Moon
#'   state (numeric); see \code{\link{astro_obliquity}} to obtain it for a given age.
#'   \code{Default=AC$obl0_deg}, i.e. today's obliquity, 23.2545 degrees.
#' @param k_range Accepted range of the precession constant k in kyr (numeric vector \code{c(min, max)}); draws
#'   outside it are rejected. \code{Default=c(10, 34)}.
#' @param seed Seed of the random number generator (integer), passed to \code{set.seed()} so that the
#'   Monte Carlo results can be reproduced. \code{Default=1}.
#' @param gs Fundamental secular frequencies in arcsec/yr (list with numeric vectors \code{g} and \code{s} of
#'   length 9, as in Laskar et al. 2004); g2-g5 and s3 are perturbed with the Laskar et al. (2011)
#'   uncertainties in each draw. \code{Default=LA2004}, the La2004 values stored inside the function.
#' @return list(k = summary of k, state = summary of the Earth-Moon state,
#'              draws = the raw draws)
#' @export
astro_invert <- function(windows, nsim = 5000, obl_deg = AC$obl0_deg,
    k_range = c(10, 34), seed = 1, gs = LA2004){
  # Internal helper functions (local to this function)

  freq_arcsec <- function(period_kyr){
    ARCSEC_PER_CYCLE / (period_kyr * 1000)
  }

  ## the secular offset of each k-bearing term, so that f_term = p + offset
  K_TERMS <- function(gs = LA2004){
    c("k+s3" = unname(gs$s[3]), "k+g5" = unname(gs$g[5]),
        "k+g2" = unname(gs$g[2]), "k+g4" = unname(gs$g[4]),
        "k+g3" = unname(gs$g[3]))
  }

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

  # Secular frequencies
  ## La2004, Laskar et al. (2004) table 3, in arcsec/yr.  These are the values
  ## Acycle bundles as La04gisi.csv and uses in calculate_orbit9.m.

  LA2004 <- list(
      g = c(g1 = 5.59, g2 = 7.452, g3 = 17.368, g4 = 17.916, g5 = 4.257452,
          g6 = 28.245, g7 = 3.087951, g8 = 0.673021, g9 = -0.34994),
      s = c(s1 = -5.59, s2 = -7.05, s3 = -18.85, s4 = -17.755, s5 = 0.00000013,
          s6 = -26.347855, s7 = -2.9925259, s8 = -0.691736, s9 = -0.34998))

  ## 1-sigma uncertainties on the secular frequencies, arcsec/yr
  ## (Laskar et al. 2011, their table 6; the same values Zeeden et al. 2023 use).
  LA_SIGMA <- c(g2 = 0.019, g3 = 0.2, g4 = 0.2, g5 = 0.000030, s3 = 0.066)

  # Main code
  set.seed(seed)
  terms <- names(windows)
  stopifnot(all(terms %in% names(K_TERMS(gs))))
  ks <- numeric(nsim)
  for (i in seq_len(nsim)){
    ## sample the secular frequencies within their uncertainties
    g <- gs$g
    s <- gs$s
    s[3]  <- s[3]  + rnorm(1, 0, LA_SIGMA[["s3"]])
    g[2]  <- g[2]  + rnorm(1, 0, LA_SIGMA[["g2"]])
    g[3]  <- g[3]  + rnorm(1, 0, LA_SIGMA[["g3"]])
    g[4]  <- g[4]  + rnorm(1, 0, LA_SIGMA[["g4"]])
    g[5]  <- g[5]  + rnorm(1, 0, LA_SIGMA[["g5"]])
    off <- c("k+s3" = unname(s[3]), "k+g5" = unname(g[5]),
        "k+g2" = unname(g[2]), "k+g4" = unname(g[4]), "k+g3" = unname(g[3]))
    ## sample one period uniformly inside each measured window and solve for p
    p_i <- numeric(length(terms))
    for (j in seq_along(terms)){
      w <- windows[[terms[j]]]
      p_i[j] <- freq_arcsec(runif(1, min(w), max(w))) - off[[terms[j]]]
    }
    names(p_i) <- terms
    ks[i] <- ARCSEC_PER_CYCLE / mean(p_i) / 1000
  }
  ks <- ks[is.finite(ks) & ks >= k_range[1] & ks <= k_range[2]]
  st <- em_state(ks, obl_deg = obl_deg)
  q <- function(x){
    c(mean = mean(x, na.rm = TRUE),
        `2.5%` = unname(quantile(x, 0.025, na.rm = TRUE)),
        `97.5%` = unname(quantile(x, 0.975, na.rm = TRUE)))
  }
  state_cols <- c("a_RE", "a_km", "lod_solar_h", "days_per_year",
      "synodic_months_per_year")
  state <- matrix(NA_real_, length(state_cols), 3,
      dimnames = list(state_cols, c("mean", "2.5%", "97.5%")))
  for (j in seq_along(state_cols)){
    state[j, ] <- q(st[[state_cols[j]]])
  }
  list(k = q(ks),
      state = state,
      draws = st, windows = windows, obl_deg = obl_deg)
}
