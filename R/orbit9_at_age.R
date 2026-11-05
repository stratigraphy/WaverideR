#' @title Nine orbital periods (kyr) at age_ma
#' @description The nine periods used by Acycle's COCO/eCOCO GUI (see calculate_orbit9()).
#' @param age_ma Age at which the nine periods are calculated, in Ma (single finite numeric value, 0 to 4429); it is
#' rounded to the nearest whole Ma (see \code{\link{calculate_orbit9}}). Required (no default).
#' @param data_dir Folder that contains the data tables \code{AstroGeo22.txt} and \code{La04gisi.csv} (character
#' path). \code{NULL} uses \code{getOption("cocoR.data")} or the copies installed with WaverideR.
#' \code{Default=NULL}.
#' @return numeric vector of nine periods (kyr).
#' @export
orbit9_at_age <- function(age_ma, data_dir = NULL){
  calculate_orbit9(age_ma, data_dir)$period_kyr
}
