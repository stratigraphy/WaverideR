#' @title Smooth sedimentation rate curve with uncertainty from a flmw score matrix
#'
#' @description Turns the sedimentation rate x depth score matrix of \code{\link{flmw}} into a continuous
#' sedimentation rate curve. The sedimentation rate is not allowed to jump around from one depth to the next: the
#' function finds the most probable path through the score matrix (Viterbi algorithm) and also the probability of every
#' sedimentation rate at every depth (forward-backward algorithm of a hidden Markov model). The probabilities give
#' credible intervals, and prior bounds or known hiatuses can be supplied.
#'
#' @param x Result of \code{\link{flmw}} (any \code{method}).
#' @param jump_sd Smoothness: standard deviation of the relative (log) change in sedimentation rate over 1 metre of
#' depth, scaled with the square root of the depth step (random walk). Small values give a stiff curve
#' \code{Default=0.2}.
#' @param sharpness How strongly the score is trusted. The likelihood of a sedimentation rate is
#' \code{exp(sharpness * score)} with scores between 0 and 1. Lower values give wider intervals \code{Default=3}.
#' @param normalize Divide the scores by their 99th percentile (and cap at 1) so that \code{sharpness} means the same
#' for every method \code{Default=TRUE}.
#' @param sedrate_range Optional vector \code{c(min, max)} (cm/kyr): sedimentation rates outside this range are
#' excluded.
#' @param hiatus Optional depths of hiatuses or other known breaks. At these depths the sedimentation rate may change
#' freely.
#' @param probs Probabilities for the credible interval limits \code{Default=c(0.025, 0.16, 0.84, 0.975)}.
#' @param genplot Plot the posterior probability with the path \code{Default=FALSE}.
#' @param x_lab Label of the x-axis \code{Default="depth"}.
#' @param y_lab Label of the y-axis \code{Default="sedrate"}.
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
#' @param normalise British spelling of \code{normalize}; if given it overrides \code{normalize}.
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
  if (!is.null(normalise)) normalize <- normalise  # British spelling synonym
  sc <- x$score_mat
  if (is.null(sc)) stop("x has no score_mat: use the result of flmw()")
  depth <- as.numeric(x$depth)
  sr <- as.numeric(x$y_axis)
  nsr <- length(sr)
  nd <- length(depth)
  sc[!is.finite(sc)] <- 0
  if (normalize) {                       # scores of different methods have different scales
    q <- stats::quantile(sc[sc > 0], 0.99, names = FALSE)
    if (is.finite(q) && q > 0) sc <- pmin(sc / q, 1)
  }
  loge <- sharpness * sc
  if (!is.null(sedrate_range)) {
    out <- sr < min(sedrate_range) | sr > max(sedrate_range)
    if (all(out)) stop("sedrate_range excludes all trial sedimentation rates")
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
    hi_idx <- unique(vapply(hiatus, function(h) {
      i <- which(depth >= h)[1]
      if (is.na(i) || i < 2) NA_integer_ else i
    }, integer(1)))
    hi_idx <- hi_idx[!is.na(hi_idx)]
  }
  # Viterbi (most probable path)
  delta <- matrix(-Inf, nsr, nd)
  psi <- matrix(1L, nsr, nd)
  delta[, 1] <- loge[, 1]
  for (t in seq_len(nd)[-1]) {
    lt <- if (t %in% hi_idx) logflat else logTm
    cand <- delta[, t - 1] + lt
    arg <- max.col(t(cand), ties.method = "first")
    psi[, t] <- arg
    delta[, t] <- cand[cbind(arg, seq_len(nsr))] + loge[, t]
  }
  path_i <- integer(nd)
  path_i[nd] <- which.max(delta[, nd])
  for (t in rev(seq_len(nd - 1))) path_i[t] <- psi[path_i[t + 1], t + 1]
  # forward-backward (probabilities)
  e <- exp(sweep(loge, 2, apply(loge, 2, max), "-"))
  norm1 <- function(v) {
    s <- sum(v)
    if (!is.finite(s) || s <= 0) rep(1 / length(v), length(v)) else v / s
  }
  alpha <- matrix(0, nsr, nd)
  alpha[, 1] <- norm1(e[, 1])
  for (t in seq_len(nd)[-1]) {
    Tt <- if (t %in% hi_idx) flat else Tm
    alpha[, t] <- norm1(as.numeric(crossprod(Tt, alpha[, t - 1])) * e[, t])
  }
  beta <- matrix(1, nsr, nd)
  for (t in rev(seq_len(nd - 1))) {
    Tt <- if ((t + 1) %in% hi_idx) flat else Tm
    beta[, t] <- norm1(as.numeric(Tt %*% (e[, t + 1] * beta[, t + 1])))
  }
  post <- alpha * beta
  post <- apply(post, 2, norm1)
  post <- matrix(post, nsr, nd)
  cs <- apply(post, 2, cumsum)
  cs <- matrix(cs, nsr, nd)
  interval <- vapply(probs, function(p) sr[pmin(colSums(cs < p - 1e-12) + 1, nsr)], numeric(nd))
  interval <- matrix(interval, nrow = nd, dimnames = list(NULL, paste0("p", format(probs, trim = TRUE))))
  res <- list(depth = depth,
              y_axis = sr,
              path = sr[path_i],
              mean = colSums(post * sr),
              interval = interval,
              posterior = post,
              score_mat = x$score_mat)
  if (genplot) {
    image(x = depth, y = sr, z = t(post), col = hcl.colors(50, "YlGnBu", rev = TRUE),
          xlab = x_lab, ylab = y_lab, useRaster = TRUE)
    lines(depth, res$path, lwd = 2)
  }
  res
}
