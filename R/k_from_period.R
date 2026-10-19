## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

## ---- 5. inversion: measured periods -> k -----------------------------------

#' Precession constant implied by one measured cycle period
#' @param period_kyr measured period, kyr
#' @param term one of k+s3, k+g5, k+g2, k+g4, k+g3
#' @param gs Vector of fundamental g frequencies (`NULL` = LA2004 values).
#' @export
k_from_period <- function(period_kyr, term = "k+s3", gs = LA2004){
  # ---- internal helper functions (local to this function) ----

  freq_arcsec <- function(period_kyr){
    ARCSEC_PER_CYCLE / (period_kyr * 1000)
  }

  ## the secular offset of each k-bearing term, so that f_term = p + offset
  K_TERMS <- function(gs = LA2004){
    c("k+s3" = unname(gs$s[3]), "k+g5" = unname(gs$g[5]),
      "k+g2" = unname(gs$g[2]), "k+g4" = unname(gs$g[4]),
      "k+g3" = unname(gs$g[3]))
  }

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  ## ---- 2. secular frequencies ------------------------------------------------
  ## La2004, Laskar et al. (2004) table 3, in arcsec/yr.  These are the values
  ## Acycle bundles as La04gisi.csv and uses in calculate_orbit9.m.

  LA2004 <- list(
    g = c(g1 = 5.59, g2 = 7.452, g3 = 17.368, g4 = 17.916, g5 = 4.257452,
          g6 = 28.245, g7 = 3.087951, g8 = 0.673021, g9 = -0.34994),
    s = c(s1 = -5.59, s2 = -7.05, s3 = -18.85, s4 = -17.755, s5 = 0.00000013,
          s6 = -26.347855, s7 = -2.9925259, s8 = -0.691736, s9 = -0.34998))

  # ---- main code ----
  off <- K_TERMS(gs)[[term]]
  p <- freq_arcsec(period_kyr) - off
  ARCSEC_PER_CYCLE / p / 1000
}
