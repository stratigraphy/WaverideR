#' @title Smooth sedimentation rate curve with uncertainty from a flmw score matrix
#'
#' @description Turns the sedimentation rate x depth score matrix of \code{\link{flmw}} into a continuous
#' sedimentation rate curve. The sedimentation rate is not allowed to jump around from one depth to the next: the
#' function finds the most probable path through the score matrix (Viterbi algorithm) and also the probability of every
#' sedimentation rate at every depth (forward-backward algorithm of a hidden Markov model). The probabilities give
#' credible intervals, and prior bounds or known hiatuses can be supplied.
#'
#' @param x Output of \code{\link{flmw}} (any \code{method}), or any list with the elements \code{score_mat} (matrix,
#' trial sedimentation rates in rows, depths in columns), \code{depth} and \code{y_axis} (trial sedimentation rates).
#' Required (no default).
#' @param jump_sd Smoothness of the curve: standard deviation of the relative (natural-log) change in sedimentation rate
#' over 1 metre of depth, scaled with the square root of the depth step (random walk) (numeric, > 0). Small values give
#' a stiff curve, large values let the rate change quickly. \code{Default=0.2}.
#' @param sharpness How strongly the score is trusted (numeric, > 0): the likelihood of a sedimentation rate is
#' \code{exp(sharpness * score)} with scores between 0 and 1. Lower values give wider credible intervals.
#' \code{Default=3}.
#' @param normalize Divide the scores by their 99th percentile (and cap them at 1) so that \code{sharpness} means the
#' same for every \code{flmw} method (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param sedrate_range Prior bounds \code{c(min, max)} of the sedimentation rate in cm/kyr (numeric vector of
#' length 2): trial rates outside this range are excluded. \code{NULL} uses all trial rates. \code{Default=NULL}.
#' @param hiatus Depths of hiatuses or other known breaks, in the depth unit of \code{x} (numeric vector). At these
#' depths the sedimentation rate may jump freely. \code{NULL} means no breaks. \code{Default=NULL}.
#' @param probs Cumulative probabilities of the credible interval limits that are returned for every depth (numeric
#' vector, values between 0 and 1). \code{Default=c(0.025, 0.16, 0.84, 0.975)} (95 and 68 percent intervals).
#' @param genplot Plot the posterior probability of the sedimentation rate against depth with the most probable path
#' (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param x_lab Label of the x-axis (depth) of the plot (character). \code{Default="depth"}.
#' @param y_lab Label of the y-axis (sedimentation rate) of the plot (character). \code{Default="sedrate"}.
#'
#' @return A list with \code{depth}, \code{y_axis} (trial sedimentation rates), \code{path} (most probable
#' sedimentation rate curve), \code{mean} (posterior mean), \code{interval} (matrix with the credible interval limits
#' for every depth), \code{posterior} (matrix of probabilities, rates x depths) and \code{score_mat}.
#'
#' @examples
#' # synthetic score matrix with a gradual change in sedimentation rate
#' depth <- seq(0, 20, by = 0.1)
#' rates <- seq(1, 6, by = 0.25)
#' true_rate <- 2 + 3 * depth / 20
#' score <- outer(rates, true_rate, function(r, t) exp(-0.5 * ((r - t) / 0.4)^2))
#' set.seed(1)
#' score <- pmin(pmax(score + matrix(rnorm(length(score), 0, 0.2), nrow(score)), 0), 1)
#' fake <- list(score_mat = score, depth = depth, y_axis = rates)
#' res <- flmw_path(fake)
#' head(res$path)
#'
#' @param normalise British spelling of \code{normalize} (logical, TRUE or FALSE); if given (not \code{NULL}) it
#' overrides \code{normalize}. \code{Default=NULL}.
#' @export
#' @importFrom graphics image lines
#' @importFrom grDevices hcl.colors
flmw_path <- function(x,
    jump_sd = 0.2,
    sharpness = 3,
    normalize = TRUE,
    sedrate_range = NULL,
    hiatus = NULL,
    probs = c(0.025, 0.16, 0.84, 0.975),
    genplot = FALSE,
    x_lab = "depth",
    y_lab = "sedrate",
    normalise = NULL) {
  if (!is.null(normalise)) {
    normalize <- normalise
  }  # British spelling synonym
  sc <- x$score_mat
  if (is.null(sc)) {
    stop("x has no score_mat: use the result of flmw()")
  }
  depth <- as.numeric(x$depth)
  sr <- as.numeric(x$y_axis)
  nsr <- length(sr)
  nd <- length(depth)
  sc[!is.finite(sc)] <- 0
  if (normalize) {                       # scores of different methods have different scales
    q <- stats::quantile(sc[sc > 0], 0.99, names = FALSE)
    if (is.finite(q) && q > 0) {
      sc <- pmin(sc / q, 1)
    }
  }
  loge <- sharpness * sc
  if (!is.null(sedrate_range)) {
    out <- sr < min(sedrate_range) | sr > max(sedrate_range)
    if (all(out)) {
      stop("sedrate_range excludes all trial sedimentation rates")
    }
    loge[out, ] <- -Inf
  }
  dz <- stats::median(abs(diff(depth)))
  sd_step <- jump_sd * sqrt(dz)
  lr <- log(sr)
  logT <- -0.5 * (outer(lr, lr, "-") / sd_step)^2        # [from, to]
  Tm <- exp(logT)
  Tm <- Tm / rowSums(Tm)
  logTm <- log(Tm)
  flat <- matrix(1 / nsr, nsr, nsr)
  logflat <- log(flat)
  hi_idx <- integer(0)
  if (length(hiatus)) {
    hi_raw <- integer(length(hiatus))
    for (ih in seq_along(hiatus)) {
      i <- which(depth >= hiatus[[ih]])[1]
      hi_raw[ih] <- if (is.na(i) || i < 2) { NA_integer_ } else { i }
    }
    hi_idx <- unique(hi_raw)
    hi_idx <- hi_idx[!is.na(hi_idx)]
  }
  # Viterbi (most probable path)
  delta <- matrix(-Inf, nsr, nd)
  psi <- matrix(1L, nsr, nd)
  delta[, 1] <- loge[, 1]
  for (t in seq_len(nd)[-1]) {
    lt <- if (t %in% hi_idx) { logflat } else { logTm }
    cand <- delta[, t - 1] + lt
    arg <- max.col(t(cand), ties.method = "first")
    psi[, t] <- arg
    delta[, t] <- cand[cbind(arg, seq_len(nsr))] + loge[, t]
  }
  path_i <- integer(nd)
  path_i[nd] <- which.max(delta[, nd])
  for (t in rev(seq_len(nd - 1))) {
    path_i[t] <- psi[path_i[t + 1], t + 1]
  }
  # forward-backward (probabilities)
  col_max <- numeric(ncol(loge))
  for (j in seq_len(ncol(loge))) {
    col_max[j] <- max(loge[, j])
  }
  e <- exp(sweep(loge, 2, col_max, "-"))
  norm1 <- function(v) {
    s <- sum(v)
    if (!is.finite(s) || s <= 0) {
      rep(1 / length(v), length(v))
    } else {
      v / s
    }
  }
  alpha <- matrix(0, nsr, nd)
  alpha[, 1] <- norm1(e[, 1])
  for (t in seq_len(nd)[-1]) {
    Tt <- if (t %in% hi_idx) { flat } else { Tm }
    alpha[, t] <- norm1(as.numeric(crossprod(Tt, alpha[, t - 1])) * e[, t])
  }
  beta <- matrix(1, nsr, nd)
  for (t in rev(seq_len(nd - 1))) {
    Tt <- if ((t + 1) %in% hi_idx) { flat } else { Tm }
    beta[, t] <- norm1(as.numeric(Tt %*% (e[, t + 1] * beta[, t + 1])))
  }
  post <- alpha * beta
  for (j in seq_len(ncol(post))) {
    post[, j] <- norm1(post[, j])
  }
  post <- matrix(post, nsr, nd)
  cs <- post
  for (j in seq_len(ncol(post))) {
    cs[, j] <- cumsum(post[, j])
  }
  cs <- matrix(cs, nsr, nd)
  interval <- matrix(NA_real_, nd, length(probs))
  for (ip in seq_along(probs)) {
    interval[, ip] <- sr[pmin(colSums(cs < probs[[ip]] - 1e-12) + 1, nsr)]
  }
  interval <- matrix(interval, nrow = nd, dimnames = list(NULL, paste0("p", format(probs, trim = TRUE))))
  res <- list(depth = depth,
      y_axis = sr,
      path = sr[path_i],
      mean = colSums(post * sr),
      interval = interval,
      posterior = post,
      score_mat = x$score_mat)
  if (genplot) {
    image(x = depth, y = sr, z = t(post), col = .wr_hcl_cols(50, "viridis"),
        xlab = x_lab, ylab = y_lab, useRaster = TRUE)
    lines(depth, res$path, lwd = 2)
  }
  res
}
