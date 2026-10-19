#' @title Build the CTest input table from interpreted cycle wavelengths
#'
#' @param age c(min, max) age of the record in Ma
#' @param e405,e125,e95 c(min, max) wavelength (m) of the 405-kyr, ~125-kyr and
#'   ~95-kyr eccentricity cycles in the depth-domain spectrum
#' @param obl c(min, max) wavelength (m) of the obliquity cycle (k+s3)
#' @param p1,p2,p3 c(min, max) wavelength (m) of the precession cycles, in the
#'   row order of the published spreadsheet, i.e. by DECREASING period:
#'   p1 = longest (k+g5, ~23.7 kyr today), p2 = k+g2 (~22.4 kyr),
#'   p3 = shortest (k+g4 and k+g3 together, ~19 kyr).
#'   Note: the comments inside the published functions label these rows the
#'   other way round, but the code tests them in the order used here (row 6
#'   against k+g5, row 7 against k+g2, row 8 against k+g4 and k+g3).
#' @return the 8x2 matrix used by ctest(), in the published layout
#'   (column 1 = lower frequency = longer wavelength)
#' @export
ctest_input <- function(age, e405, e125, e95, obl = c(NA, NA),
                        p1 = c(NA, NA), p2 = c(NA, NA), p3 = c(NA, NA)){
  # frequencies, low then high
  wl <- function(x){
    c(1 / max(x), 1 / min(x))
  }
  m <- rbind(age, wl(e405), wl(e125), wl(e95), wl(obl), wl(p1), wl(p2), wl(p3))
  dimnames(m) <- list(c("age_Ma", "e405", "e125", "e95", "obliquity",
                        "p1", "p2", "p3"), c("min", "max"))
  m
}
