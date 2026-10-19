#' @title Best sedimentation rate per eCOCO window
#' @description Row of maximum finite score in every window
#'   (bestPerWindowSummary in Acycle ecoco.m).
#' @param sr sedimentation-rate grid (cm/kyr)
#' @param depth window centres (m)
#' @param rho,p,norbit,score rate-by-window matrices
#' @return data.frame with depth, sr, rho, p_global, norbit, score, sr_low, sr_high
#' @export
ecoco_best_per_window <- function(sr, depth, rho, p, norbit, score){
  # ---- internal helper functions (local to this function) ----

  .best_finite_index <- function(v){
    ok <- which(is.finite(v))
    if (!length(ok)){
      return(integer(0))
    }
    ok[which.max(v[ok])]
  }

  # ---- main code ----
  nw <- length(depth)
  out <- data.frame(depth = depth, sr = NA_real_, rho = NA_real_, p_global = NA_real_,
                    norbit = NA_real_, score = NA_real_, sr_low = NA_real_, sr_high = NA_real_)
  for (j in seq_len(nw)){
    b <- .best_finite_index(score[, j])
    if (!length(b)){
      next
    }
    out[j, 2:6] <- c(sr[b], rho[b, j], p[b, j], norbit[b, j], score[b, j])
  }
  out
}
