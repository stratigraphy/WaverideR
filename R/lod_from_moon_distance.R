#' Length of day implied by an Earth-Moon distance (angular momentum conserved)
#' @param a_RE Earth-Moon distance (semi-major axis of the lunar orbit) in Earth radii (numeric vector, > 0; the
#' present-day value is about 60.14). Required (no default).
#' @param kind Type of day that is returned (character, passed to \code{\link{lod_hours}}): \code{"solar"}
#' (mean solar day) or \code{"sidereal"} (one rotation relative to the stars). \code{Default="solar"}.
#' @param L_total Total angular momentum of the Earth-Moon system in kg m^2/s (numeric, > 0), assumed to be
#' conserved. \code{NULL} uses the present-day value (computed from today's Earth-Moon distance and length of day).
#' \code{Default=NULL}.
#' @export
lod_from_moon_distance <- function(a_RE, kind = "solar", L_total = NULL){
  lod_hours(omega_from_moon_distance(a_RE, L_total), kind)
}
