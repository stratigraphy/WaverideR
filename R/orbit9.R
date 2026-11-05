## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

# The nine principal terms (Acycle calculate_orbit9 equivalent)

#' Nine principal astronomical terms for a given precession constant
#' @param k_kyr Precession constant k expressed as the period of the axial precession, in kyr (numeric, > 0; about
#' 25.7 kyr today). Required (no default).
#' @param gs Secular frequencies in arcsec/yr: a list with elements \code{g} and \code{s}, each a numeric vector of
#' the nine fundamental frequencies g1-g9 and s1-s9. The default \code{LA2004} is the La2004 table (Laskar et al.,
#' 2004) defined inside the function. \code{Default=LA2004}.
#' @return data frame: term, group, frequency (arcsec/yr), period (kyr)
#' @export
orbit9 <- function(k_kyr, gs = LA2004){
  # Internal helper functions (local to this function)

  period_kyr  <- function(freq_arcsec){
    ARCSEC_PER_CYCLE / abs(freq_arcsec) / 1000
  }

  freq_arcsec <- function(period_kyr){
    ARCSEC_PER_CYCLE / (period_kyr * 1000)
  }

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  # Secular frequencies
  ## La2004, Laskar et al. (2004) table 3, in arcsec/yr.  These are the values
  ## Acycle bundles as La04gisi.csv and uses in calculate_orbit9.m.

  LA2004 <- list(
      g = c(g1 = 5.59, g2 = 7.452, g3 = 17.368, g4 = 17.916, g5 = 4.257452,
          g6 = 28.245, g7 = 3.087951, g8 = 0.673021, g9 = -0.34994),
      s = c(s1 = -5.59, s2 = -7.05, s3 = -18.85, s4 = -17.755, s5 = 0.00000013,
          s6 = -26.347855, s7 = -2.9925259, s8 = -0.691736, s9 = -0.34998))

  # Main code
  g <- gs$g
  s <- gs$s
  p <- freq_arcsec(k_kyr)
  f <- c(g[2] - g[5], g[3] - g[2], g[4] - g[2], g[3] - g[5], g[4] - g[5],
      p + s[3], p + g[5], p + g[2], p + g[4])
  nm <- c("g2-g5", "g3-g2", "g4-g2", "g3-g5", "g4-g5",
      "k+s3", "k+g5", "k+g2", "k+g4")
  grp <- c("long eccentricity", rep("short eccentricity", 4),
      "obliquity", rep("precession", 3))
  data.frame(term = nm, group = grp, freq_arcsec_yr = unname(f),
      period_kyr = unname(period_kyr(f)), stringsAsFactors = FALSE)
}
