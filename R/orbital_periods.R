# Forward: precession constant + secular frequencies -> periods

#' All principal astronomical periods for a given precession constant
#'
#' @param k_yr Precession constant k expressed as the period of the axial precession, in years (numeric, > 0), e.g.
#' 20700 for the Devonian estimate of Zeeden et al. (2023); ignored when \code{p_arcsec} is given. One of
#' \code{k_yr} and \code{p_arcsec} is required. \code{Default=NULL}.
#' @param p_arcsec Precession constant as a frequency in arcsec/yr (numeric, > 0; about 50.5 today); if given it is
#' used instead of \code{k_yr}. \code{Default=NULL}.
#' @param g Secular frequencies g1-g9 of the planetary perihelia in arcsec/yr (numeric vector of length 9).
#' \code{NULL} uses the La2004 values of \code{La04gisi.csv}. \code{Default=NULL}.
#' @param s Secular frequencies s1-s9 of the planetary nodes in arcsec/yr (numeric vector of length 9).
#' \code{NULL} uses the La2004 values of \code{La04gisi.csv}. \code{Default=NULL}.
#' @param terms Which terms are returned (character): \code{"nine"} (the Acycle ORBIT9 set: g2-g5, g3-g2, g4-g2,
#' g3-g5, g4-g5, k+s3, k+g5, k+g2, k+g4) or \code{"all"} (every gi-gj eccentricity term and every k+si and k+gi term
#' for i = 1 to 5). \code{Default=c("nine", "all")}, i.e. \code{"nine"}.
#' @param data_dir Folder that contains the data tables \code{AstroGeo22.txt} and \code{La04gisi.csv} (character
#' path). \code{NULL} uses \code{getOption("cocoR.data")} or the copies installed with WaverideR.
#' \code{Default=NULL}.
#' @return data.frame: term, group, frequency (arcsec/yr), period (kyr)
#' @export
orbital_periods <- function(k_yr = NULL, p_arcsec = NULL, g = NULL, s = NULL,
    terms = c("nine", "all"), data_dir = NULL){
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
  terms <- match.arg(terms)
  if (is.null(p_arcsec)){
    if (is.null(k_yr)){
      stop("Give k_yr or p_arcsec")
    }
    p_arcsec <- 360 * 3600 / k_yr
  }
  if (is.null(g) || is.null(s)){
    gs <- .la04_gisi(data_dir)
    if (is.null(g)){
      g <- gs$g
    }
    if (is.null(s)){
      s <- gs$s
    }
  }
  if (terms == "nine"){
    freq <- c(g[2] - g[5], g[3] - g[2], g[4] - g[2], g[3] - g[5], g[4] - g[5],
        p_arcsec + s[3], p_arcsec + g[5], p_arcsec + g[2], p_arcsec + g[4])
    nm <- c("g2-g5", "g3-g2", "g4-g2", "g3-g5", "g4-g5", "k+s3", "k+g5", "k+g2", "k+g4")
    grp <- c("long eccentricity", rep("short eccentricity", 4), "obliquity",
        rep("precession", 3))
  } else {
    ecc <- utils::combn(5, 2)
    nm <- freq <- grp <- NULL
    for (j in seq_len(ncol(ecc))){
      i1 <- ecc[2, j]
      i2 <- ecc[1, j]
      nm <- c(nm, sprintf("g%d-g%d", i1, i2))
      freq <- c(freq, g[i1] - g[i2])
      grp <- c(grp, "eccentricity")
    }
    for (i in seq_len(5)){
      nm <- c(nm, sprintf("k+s%d", i))
      freq <- c(freq, p_arcsec + s[i])
      grp <- c(grp, "obliquity")
    }
    for (i in seq_len(5)){
      nm <- c(nm, sprintf("k+g%d", i))
      freq <- c(freq, p_arcsec + g[i])
      grp <- c(grp, "precession")
    }
  }
  out <- data.frame(term = nm, group = grp, freq_arcsec_yr = freq,
      period_kyr = 360 * 3600 / abs(freq) / 1000)
  out[order(-out$period_kyr), ]
}
