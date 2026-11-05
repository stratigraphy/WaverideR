#' @title Nine astronomical periods at a given age
#' @description Nine principal orbital periods (calculate_orbit9.m in Acycle):
#'   La2004 secular frequencies g_i, s_i and the AstroGeo22 precession constant
#'   k(t) at the age rounded to the nearest Ma.
#' @param age_ma Age at which the nine periods are calculated, in Ma (single finite numeric value, 0 to 4429);
#' it is rounded to the nearest whole Ma because AstroGeo22 is tabulated in 1-Myr steps. Required, no default.
#' @param data_dir Folder that contains the astronomical data tables \code{AstroGeo22.txt} and
#' \code{La04gisi.csv} (character path); \code{NULL} uses the tables installed with the package (or the folder
#' set with \code{options(cocoR.data = ...)}). \code{Default=NULL}.
#' @return data.frame with columns term, group, freq_arcsec_yr, period_yr, period_kyr.
#' @export
calculate_orbit9 <- function(age_ma, data_dir = NULL){
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
  if (!is.numeric(age_ma) || length(age_ma) != 1 || !is.finite(age_ma)){
    stop("age_ma must be a finite number")
  }
  if (is.null(data_dir)){
    data_dir <- coco_data_dir()
  }
  if (is.null(.coco_env$astro)){
    .coco_env$astro <- utils::read.table(file.path(data_dir, "AstroGeo22.txt"), comment.char = "%")
    gs <- readLines(file.path(data_dir, "La04gisi.csv"), warn = FALSE, encoding = "UTF-8")
    gs <- gs[!grepl("%", gs) & nzchar(trimws(gs))]
    .coco_env$gisi <- as.numeric(strsplit(gs[1], ",")[[1]])[1:18]
  }
  astro <- .coco_env$astro
  age_ma <- round(age_ma)
  d <- abs(astro[, 1] - age_ma / 1000)
  row <- which.min(d)
  if (d[row] > 1e-8){
    stop(sprintf("No AstroGeo22 row for %d Ma (valid 0-%d Ma)", age_ma, round(1000 * max(astro[, 1]))))
  }
  k <- astro[row, 11]
  g <- .coco_env$gisi[1:9]
  s <- .coco_env$gisi[10:18]
  freq <- c(g[2] - g[5], g[3] - g[2], g[4] - g[2], g[3] - g[5], g[4] - g[5],
      k + s[3], k + g[5], k + g[2], k + g[4])
  period <- 360 * 3600 / abs(freq)
  data.frame(term = c("g2-g5", "g3-g2", "g4-g2", "g3-g5", "g4-g5", "k+s3", "k+g5", "k+g2", "k+g4"),
      group = c("long eccentricity", rep("short eccentricity", 4), "obliquity", rep("precession", 3)),
      freq_arcsec_yr = freq, period_yr = period, period_kyr = period / 1000)
}
