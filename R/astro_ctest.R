## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

#' Consistency-test inversion (the estimator used in the paper)
#'
#' After Zeeden et al. (2023): draw a precession constant and a set of secular
#' frequencies, predict every term, and keep the draw only if all measured
#' windows are satisfied at once.  Two details matter.  First, the secular
#' frequencies are sampled as periods with the Laskar et al. (2011)
#' uncertainties scaled by sqrt(age / 100 Myr).  Second, because a record tuned
#' to a nominal 405 kyr inherits the error in the true g2-g5 period, every
#' measured period is rescaled by (sampled 405-kyr period) / 405 before it is
#' compared with the prediction; this is what allows the tuned time scale
#' itself to stretch slightly.
#'
#' @param windows named list, term -> c(low, high) in kyr, as in astro_invert()
#' @param age_ma age of the record, or c(min, max)
#' @param nsim number of accepted draws to collect
#' @param k_range prior range for the precession constant, kyr
#' @param anchor_kyr the period assumed when the record was tuned
#' @param seed Random seed for reproducibility.
#' @param block Number of random draws per block.
#' @param max_attempts Maximum number of draws before the search is aborted.
#' @param verbose Print progress messages.
#' @export
astro_ctest <- function(windows, age_ma, nsim = 3000, k_range = c(10, 34),
                        anchor_kyr = 405, seed = 1, block = 2e5,
                        max_attempts = 5e7, verbose = TRUE){
  set.seed(seed)
  if (length(age_ma) == 1){
    age_ma <- age_ma + c(-5, 5)
  }
  terms <- names(windows)
  ## La2004 secular frequencies as periods in years, and their relative
  ## 1-sigma uncertainties (Laskar et al. 2011)
  P <- c(g2 = 173913, g3 = 74620, g4 = 72338, g5 = 304407, s3 = 68753)
  drel <- c(g2 = 0.019 / 7.453, g3 = 0.2 / 17.368, g4 = 0.2 / 17.916,
            g5 = 0.000030 / 4.257482, s3 = 0.066 / 18.848)
  keep <- list()
  got <- 0
  tried <- 0
  while (got < nsim){
    tried <- tried + block
    if (tried > max_attempts){
      stop("no consistent draws; widen the windows, the age range or k_range")
    }
    t_yr <- runif(block, age_ma[1] * 1e6, age_ma[2] * 1e6)
    k_yr <- runif(block, k_range[1] * 1000, k_range[2] * 1000)
    sc <- sqrt(t_yr / 1e8)
    G2 <- rnorm(block, P[["g2"]], P[["g2"]] * drel[["g2"]] * sc)
    G3 <- rnorm(block, P[["g3"]], P[["g3"]] * drel[["g3"]] * sc)
    G4 <- rnorm(block, P[["g4"]], P[["g4"]] * drel[["g4"]] * sc)
    G5 <- rnorm(block, P[["g5"]], P[["g5"]] * drel[["g5"]] * sc)
    S3 <- rnorm(block, P[["s3"]], drel[["s3"]] * sc)
    per <- list("e405" = 1 / (1 / G2 - 1 / G5) / 1000,
                "k+s3" = 1 / (1 / k_yr - 1 / S3) / 1000,
                "k+g5" = 1 / (1 / k_yr + 1 / G5) / 1000,
                "k+g2" = 1 / (1 / k_yr + 1 / G2) / 1000,
                "k+g4" = 1 / (1 / k_yr + 1 / G4) / 1000,
                "k+g3" = 1 / (1 / k_yr + 1 / G3) / 1000)
    scale <- per$e405 / anchor_kyr       # the tuned axis may stretch
    ok <- rep(TRUE, block)
    for (tm in terms){
      w <- sort(windows[[tm]])
      obs <- per[[tm]] / scale
      ok <- ok & is.finite(obs) & obs >= w[1] & obs <= w[2]
    }
    idx <- which(ok)
    if (!length(idx)){
      next
    }
    idx <- idx[seq_len(min(length(idx), nsim - got))]
    keep[[length(keep) + 1]] <- data.frame(
      k_kyr = k_yr[idx] / 1000, age_ma = t_yr[idx] / 1e6, scale = scale[idx],
      obl_ks3 = per$`k+s3`[idx], prec_kg5 = per$`k+g5`[idx],
      prec_kg2 = per$`k+g2`[idx], prec_kg4 = per$`k+g4`[idx],
      prec_kg3 = per$`k+g3`[idx], ecc405 = per$e405[idx])
    got <- got + length(idx)
  }
  smp <- do.call(rbind, keep)
  if (verbose){
    message(sprintf("astro_ctest: %d accepted of %d draws (%.3f%%)",
                    nrow(smp), tried, 100 * nrow(smp) / tried))
  }
  st <- em_state(smp$k_kyr, obl_deg = astro_obliquity(mean(age_ma)))
  q <- function(x){
    c(mean = mean(x, na.rm = TRUE),
      `2.5%` = unname(quantile(x, 0.025, na.rm = TRUE)),
      `97.5%` = unname(quantile(x, 0.975, na.rm = TRUE)))
  }
  state_cols <- c("a_RE", "a_km", "lod_solar_h", "days_per_year",
                  "synodic_months_per_year")
  state <- matrix(NA_real_, length(state_cols), 3,
                  dimnames = list(state_cols, c("mean", "2.5%", "97.5%")))
  for (j in seq_along(state_cols)){
    state[j, ] <- q(st[[state_cols[j]]])
  }
  list(k = q(smp$k_kyr),
       state = state,
       draws = cbind(smp, st[, -1]), windows = windows, samples = smp,
       obl_deg = astro_obliquity(mean(age_ma)))
}
