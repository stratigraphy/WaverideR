#' @title Best sedimentation rate per eCOCO window
#' @description Row of maximum finite score in every window
#'   (bestPerWindowSummary in Acycle ecoco.m).
#' @param sr Sedimentation-rate grid of the score map, in cm/kyr (numeric vector, one value per row of the
#' matrices). Required, no default.
#' @param depth Depths of the window centres, in metres (numeric vector, one value per column of the matrices).
#' Required, no default.
#' @param rho Correlation coefficients (numeric matrix, rates in rows, windows in columns), e.g. element \code{rho}
#' of the output of \code{\link{ecoco_adaptive}}; the value at the best rate is returned. Required, no default.
#' @param p Global Monte Carlo p-values (numeric matrix of the same size as \code{rho}), e.g. element
#' \code{p_global}; the value at the best rate is returned. Required, no default.
#' @param norbit Number of contributing orbital periods (numeric matrix of the same size as \code{rho}). Required,
#' no default.
#' @param score Score used to choose the best rate in each window, the rate with the highest finite score (numeric
#' matrix of the same size as \code{rho}), e.g. element \code{score} (pCOCO times the fraction of contributing
#' periods). Windows without a finite score give \code{NA}. Required, no default.
#' @return data.frame with depth, sr, rho, p_global, norbit, score, sr_low, sr_high
#' @export
ecoco_best_per_window <- function(sr, depth, rho, p, norbit, score){
  # Internal helper functions (local to this function)

  .best_finite_index <- function(v){
    ok <- which(is.finite(v))
    if (!length(ok)){
      return(integer(0))
    }
    ok[which.max(v[ok])]
  }

  # Main code
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
