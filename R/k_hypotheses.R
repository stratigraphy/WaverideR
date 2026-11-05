#' @title Which precession term is that peak? One k per hypothesis.
#'
#' @description With only one or two precession cycles the assignment matters
#' more than the measurement. This returns the precession constant and
#' Earth-Moon state implied by reading a single cycle as each of the candidate
#' terms.
#'
#' @param period_kyr Observed period of the precession (or obliquity) cycle in kyr (numeric, one value > 0).
#' Required (no default).
#' @param age_ma Approximate age of the record in Ma (numeric, one value); used for the obliquity of the
#' angular-momentum inversion and for the AstroGeo22 model value that the column \code{dk_vs_model_kyr} is compared
#' with. Required (no default).
#' @param terms Candidate astronomical terms the cycle may represent (character vector; each must contain k), e.g.
#' \code{"k+g5"}, \code{"k+g2"}, \code{"k+g4"}, \code{"k+g3"}, \code{"k+s3"}; one output row per term (names of a
#' named vector become row names). \code{Default=c("k+g5", "k+g2", "k+g4", "k+g3", "k+s3")}.
#' @param method How the precession constant is converted to Earth-Moon distance and length of day (character,
#' passed to \code{\link{earthmoon_state}}): \code{"angular_momentum"} (conservation of Earth-Moon angular
#' momentum) or \code{"astrogeo22"} (inversion of the AstroGeo22 curve of Farhat et al., 2022).
#' \code{Default=c("angular_momentum", "astrogeo22")}, i.e. \code{"angular_momentum"}.
#' @param data_dir Folder that contains the data tables \code{AstroGeo22.txt} and \code{La04gisi.csv} (character
#' path). \code{NULL} uses \code{getOption("cocoR.data")} or the copies installed with WaverideR.
#' \code{Default=NULL}.
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
