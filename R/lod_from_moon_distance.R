#' Length of day implied by an Earth-Moon distance (angular momentum conserved)
#' @param a_RE Earth-Moon distance in Earth radii.
#' @param kind Type of day, `"solar"` or `"sidereal"`.
#' @param L_total Total angular momentum (`NULL` = present-day Earth-Moon value).
#' @export
lod_from_moon_distance <- function(a_RE, kind = "solar", L_total = NULL){
  lod_hours(omega_from_moon_distance(a_RE, L_total), kind)
}
