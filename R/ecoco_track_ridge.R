#' @title Track a smooth sedimentation-rate ridge through an eCOCO score map
#' @description Within-window min-max normalised score, quadratic jump penalty
#'   (penalty 0.35, jump scale max(2, 4*sr step)), Viterbi back-tracking
#'   (trackEcocoRidge in Acycle ecoco.m).
#' @param score Score map through which the ridge is tracked (numeric matrix, rates in rows, windows in columns),
#' e.g. element \code{score} of the output of \code{\link{ecoco_adaptive}}; it is min-max normalised within each
#' window before tracking. Required, no default.
#' @param sr Sedimentation-rate grid of the score map, in cm/kyr (numeric vector, one value per row of
#' \code{score}). Required, no default.
#' @param depth Depths of the window centres, in metres (numeric vector, one value per column of \code{score}).
#' Required, no default.
#' @param rho Correlation coefficients (numeric matrix of the same size as \code{score}); the values along the
#' ridge are returned. Required, no default.
#' @param p Global Monte Carlo p-values (numeric matrix of the same size as \code{score}); the values along the
#' ridge are returned. Required, no default.
#' @param norbit Number of contributing orbital periods (numeric matrix of the same size as \code{score}); the values
#' along the ridge are returned. Required, no default.
#' @param fallback Result returned for windows the ridge does not reach, the output of
#' \code{\link{ecoco_best_per_window}} (data frame); \code{NULL} computes it from the other arguments.
#' \code{Default=NULL}.
#' @param jump_penalty Weight of the quadratic penalty on rate jumps between neighbouring windows (numeric, >= 0):
#' each step costs \code{jump_penalty * (|sr jump| / max(2, 4 * sr step))^2}. Larger values give a smoother
#' ridge, \code{0} picks the best rate in every window independently. \code{Default=0.35}.
#' @return data.frame with depth, sr, rho, p_global, norbit, score, sr_low, sr_high
#' @export
ecoco_track_ridge <- function(score, sr, depth, rho, p, norbit, fallback = NULL,
    jump_penalty = 0.35){
  # Internal helper functions (local to this function)

  .local_sr_range <- function(score_col, sr, row, rel = 0.9){
    if (row < 1 || row > length(score_col) || !is.finite(score_col[row])){
      return(c(NA, NA))
    }
    thr <- rel * score_col[row]
    ok <- is.finite(score_col) & score_col >= thr
    ok[row] <- TRUE
    lo <- row
    while (lo > 1 && ok[lo - 1]){
      lo <- lo - 1
    }
    hi <- row
    while (hi < length(ok) && ok[hi + 1]){
      hi <- hi + 1
    }
    range(sr[c(lo, hi)])
  }

  # Main code
  if (is.null(fallback)){
    fallback <- ecoco_best_per_window(sr, depth, rho, p, norbit, score)
  }
  out <- fallback
  nSr <- nrow(score)
  nWin <- ncol(score)
  if (nSr == 0 || nWin == 0 || length(sr) != nSr){
    return(out)
  }
  sn <- matrix(NA_real_, nSr, nWin)
  for (j in seq_len(nWin)){
    cs <- score[, j]
    ok <- is.finite(cs)
    if (!any(ok)){
      next
    }
    mn <- min(cs[ok])
    mx <- max(cs[ok])
    if (mx > mn){
      sn[ok, j] <- (cs[ok] - mn) / (mx - mn)
    } else {
      sn[ok, j] <- 1
    }
  }
  if (!any(is.finite(sn))){
    return(out)
  }
  srStep <- stats::median(abs(diff(sr)))
  if (!is.finite(srStep) || srStep <= 0){
    srStep <- 1
  }
  jumpScale <- max(2, 4 * srStep)
  dp <- matrix(-Inf, nSr, nWin)
  back <- matrix(NA_integer_, nSr, nWin)
  v1 <- is.finite(sn[, 1])
  dp[v1, 1] <- sn[v1, 1]
  if (nWin >= 2){
    for (j in 2:nWin){
      now <- which(is.finite(sn[, j]))
      if (!length(now)){
        next
      }
      prev <- which(is.finite(dp[, j - 1]))
      if (!length(prev)){
        dp[now, j] <- sn[now, j]
        next
      }
      for (r in now){
        tr <- dp[prev, j - 1] - jump_penalty * ((abs(sr[r] - sr[prev]) / jumpScale)^2)
        b <- which.max(tr)
        dp[r, j] <- sn[r, j] + tr[b]
        back[r, j] <- prev[b]
      }
    }
  }
  col_finite <- logical(nWin)
  for (j in seq_len(nWin)){
    col_finite[j] <- any(is.finite(dp[, j]))
  }
  anyf <- which(col_finite)
  if (!length(anyf)){
    return(out)
  }
  last <- max(anyf)
  path <- rep(NA_integer_, nWin)
  path[last] <- which.max(dp[, last])
  if (last >= 2){
    for (j in last:2){
      pv <- back[path[j], j]
      if (is.na(pv)){
        break
      }
      path[j - 1] <- pv
    }
  }
  for (j in seq_len(nWin)){
    r <- path[j]
    if (is.na(r)){
      next
    }
    rg <- .local_sr_range(score[, j], sr, r, 0.9)
    out[j, ] <- list(depth[j], sr[r], rho[r, j], p[r, j], norbit[r, j], score[r, j], rg[1], rg[2])
  }
  out
}
