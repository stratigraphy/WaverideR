#' @title Which precession term is that peak? One k per hypothesis.
#'
#' @description With only one or two precession cycles the assignment matters
#' more than the measurement. This returns the precession constant and
#' Earth-Moon state implied by reading a single cycle as each of the candidate
#' terms.
#'
#' @param period_kyr the observed cycle period (kyr)
#' @param age_ma approximate age, used for the obliquity and the model row
#' @param terms candidate terms
#' @param method "angular_momentum" or "astrogeo22"
#' @param data_dir folder with the AstroGeo22 table
#' @return data.frame, one row per candidate term
#' @export
k_hypotheses <- function(period_kyr, age_ma,
                         terms = c("k+g5", "k+g2", "k+g4", "k+g3", "k+s3"),
                         method = c("angular_momentum", "astrogeo22"), data_dir = NULL){
  method <- match.arg(method)
  rows <- vector("list", length(terms))
  names(rows) <- names(terms)
  for (j in seq_along(terms)){
    tm <- terms[[j]]
    st <- earthmoon_state(cycle_kyr = period_kyr, term = tm, age_ma = age_ma,
                          method = method, data_dir = data_dir)
    rows[[j]] <- data.frame(term = tm, cycle_kyr = period_kyr, k_kyr = st$k_kyr,
                            p_arcsec = st$p_arcsec, a_RE = st$a_RE,
                            lod_solar_h = st$lod_solar_h,
                            days_per_year = st$days_per_year,
                            moon_beyond_today = st$a_RE > 60.142611)
  }
  out <- do.call(rbind, rows)
  ref <- astrogeo22(age_ma, data_dir = data_dir)
  out$dk_vs_model_kyr <- out$k_kyr - ref$p_kyr
  out
}
