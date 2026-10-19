#' @title Nine orbital periods (kyr) at age_ma
#' @description The nine periods used by Acycle's COCO/eCOCO GUI (see calculate_orbit9()).
#' @param age_ma age in Ma.
#' @param data_dir data folder \code{Default=NULL}.
#' @return numeric vector of nine periods (kyr).
#' @export
orbit9_at_age <- function(age_ma, data_dir = NULL){
  calculate_orbit9(age_ma, data_dir)$period_kyr
}
