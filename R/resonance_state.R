# Chaotic resonance transitions in deep time
#
# The Earth-Mars secular resonance argument is
#       theta = (s4 - s3) - 2 (g4 - g3)
# which is ~0 today: the 1.2 Myr obliquity modulation is exactly twice as fast
# as the 2.4 Myr eccentricity modulation (a 2:1 lock, stable for at least the
# past 45 Myr; Laskar et al. 2004, 2011a). The Solar System is chaotic, so this
# lock can break: Ma, C., Meyers, S.R., Sageman, B.B. (2017) Nature 542, 468
# and (2019) EPSL 513, 1-11 report a switch from a 2.4 Myr to a 1.2 Myr
# expression between ~90 and 86.5 Ma in the Western Interior Seaway, and
# Olsen et al. (2019) Science 365, eaba0937 map the same behaviour in the
# Triassic-Jurassic. Astronomical solutions cannot be extrapolated across such
# a transition, which is the hard limit on deep-time astrochronology.
#
# These helpers answer two practical questions:
#   resonance_state()   given a measured 2.4 Myr and 1.2 Myr pair, which
#                       resonance state is the record in?
#   chaos_scenario()    if g3 and g4 (and s3, s4) shift, what happens to every
#                       other cycle? Some terms move a lot, some not at all.
#   term_sensitivity()  the same thing as partial derivatives.

#' Which secular resonance state does a record show?
#'
#' @param ecc_kyr Measured period of the long eccentricity modulation (the "2.4 Myr" cycle, g4-g3), in kyr (numeric,
#'   > 0, e.g. \code{2400}). Required, no default.
#' @param obl_kyr Measured period of the obliquity modulation (the "1.2 Myr" cycle, s4-s3), in kyr (numeric, > 0, e.g.
#'   \code{1200}). Required, no default.
#' @param max_ratio Largest integer p and q considered when searching the simple ratio p:q closest to the measured
#'   frequency ratio f(obliquity modulation) / f(eccentricity modulation) (integer, >= 1). \code{Default=4}.
#' @param data_dir Directory that contains the data file \code{La04gisi.csv} with the La2004 secular frequencies used
#'   for the modern values (character path). \code{NULL} uses the data installed with the package (or
#'   \code{options(cocoR.data)}). \code{Default=NULL}.
#' @return list: the implied frequencies, their ratio, the nearest simple
#'   resonance, the resonance argument (arcsec/yr and as a period), and the
#'   modern values for comparison
#' @export
resonance_state <- function(ecc_kyr, obl_kyr, max_ratio = 4, data_dir = NULL){
  # Internal helper functions (local to this function)

  coco_data_dir <- function(){
    if (!is.null(.coco_env$src)){
      src_data <- file.path(dirname(.coco_env$src), "..", "data")
    } else {
      src_data <- NULL
    }
    cand <- c(getOption("cocoR.data"), system.file("extdata", package = "WaverideR"), "data", file.path("cocoR", "data"), src_data)
    for (d in cand){
      if (!is.null(d) && file.exists(file.path(d, "AstroGeo22.txt"))){
        return(d)
      }
    }
    stop("Cannot find AstroGeo22.txt / La04gisi.csv. Reinstall WaverideR or set options(cocoR.data = '/path/to/data').")
  }

  .la04_gisi <- function(data_dir = NULL){
    if (is.null(data_dir)){
      data_dir <- coco_data_dir()
    }
    gs <- readLines(file.path(data_dir, "La04gisi.csv"), warn = FALSE, encoding = "UTF-8")
    gs <- gs[!grepl("%", gs) & nzchar(trimws(gs))]
    gisi <- as.numeric(strsplit(gs[1], ",")[[1]])[1:18]
    list(g = gisi[1:9], s = gisi[10:18])
  }

  # Astronomical periods at a given age (calculate_orbit9.m)
  # Nine principal periods from the La2004 secular frequencies g_i, s_i
  # (La04gisi.csv) and the precession constant k(t) of AstroGeo22
  # (AstroGeo22.txt, column 11, arcsec/yr), both bundled in Acycle
  # code/package/secular_periods and copied to acycle2R/data.
  #   1 g2-g5 (long ecc.)   2 g3-g2   3 g4-g2   4 g3-g5   5 g4-g5 (short ecc.)
  #   6 k+s3 (obliquity)    7 k+g5    8 k+g2    9 k+g4 (precession)
  # The age is rounded to the nearest integer Ma (AstroGeo22 has 1-Myr steps,
  # 0-4429 Ma). g and s are held at their present La2004 values.
  .coco_env <- new.env()

  .coco_env$src <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e){
    NULL
  })

  # Main code
  gs <- .la04_gisi(data_dir)
  f_ecc <- 360 * 3600 / (ecc_kyr * 1000)
  f_obl <- 360 * 3600 / (obl_kyr * 1000)
  ratio <- f_obl / f_ecc
  cands <- expand.grid(p = 1:max_ratio, q = 1:max_ratio)
  cands <- cands[abs(cands$p / cands$q - ratio) == min(abs(cands$p / cands$q - ratio)), ][1, ]
  theta <- f_obl - (cands$p / cands$q) * f_ecc
  modern_ecc <- gs$g[4] - gs$g[3]
  modern_obl <- gs$s[4] - gs$s[3]
  if (abs(theta) > 0){
    arg_period <- 360 * 3600 / abs(theta) / 1000
  } else {
    arg_period <- Inf
  }
  list(ecc_kyr = ecc_kyr, obl_kyr = obl_kyr,
      f_ecc_arcsec = f_ecc, f_obl_arcsec = f_obl,
      ratio = ratio, nearest = sprintf("%d:%d", cands$p, cands$q),
      argument_arcsec = theta,
      argument_period_kyr = arg_period,
      modern_ratio = modern_obl / modern_ecc,
      modern_ecc_kyr = 360 * 3600 / abs(modern_ecc) / 1000,
      modern_obl_kyr = 360 * 3600 / abs(modern_obl) / 1000,
      modern_state = "2:1 (g4-g3 : s4-s3 lock, La2004)",
      matches_modern = abs(ratio - modern_obl / modern_ecc) < 0.1)
}
