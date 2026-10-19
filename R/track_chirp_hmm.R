#' Chirp-steered tracking of the sedimentation rate (rate and rate gradient as the hidden state)
#'
#' Hidden Markov tracking on the chirp-matched template scores of \code{\link{analyse_chirp_eha}} (with
#' \code{criterion = "template"}): the hidden state is the pair (sedimentation rate, relative rate gradient). The
#' evidence of a window is the match of the astronomical template with the spectrum de-chirped by that gradient, at
#' that rate; the transition moves the rate by the gradient (log r' = log r + g step) and lets the gradient drift
#' slowly. The gradient therefore steers the track along the direction the de-chirped evidence points to, and a
#' change of the rate has to be announced by the chirp of the cycles before the track follows it: a wiggle of the
#' rate that the spectra support in both its level and its slope is followed, one that the slopes contradict is not.
#'
#' @param ce Result of \code{\link{analyse_chirp_eha}} with \code{criterion = "template"} or of
#' \code{\link{analyse_chirp_wavelet}}; or any scalogram (\code{\link{analyse_wavelet}},
#' \code{\link{analyse_superlet_sst}}, ...): the evidence of \code{\link{track_astro_hmm}} is then tracked with
#' the gradient as a hidden state (arguments of \code{track_astro_hmm} such as \code{astro_periods} and
#' \code{sedrate} in \code{...}, gradient grid \code{chirp}); the evidence then does not depend on the gradient,
#' only the continuity of the rate's slope steers the track.
#' @param chirp Gradient grid for a plain scalogram.
#' @param ... Arguments of \code{\link{track_astro_hmm}} for a plain scalogram.
#' @param weight Weight of the template score of one window against the transition prior \code{Default=10}. On CIP
#' Case 3 (window 20 m, step 0.6 m) the result was the same for 5-20; below 3 the track no longer follows the
#' evidence.
#' @param sd_rate Standard deviation (octaves per window step) of the rate beyond the move given by the gradient
#' \code{Default=0.01}.
#' @param sd_gradient Standard deviation of the change of the gradient per window step, in units of the gradient grid
#' spacing \code{Default=0.7}.
#' @param debias Remove, for every gradient, the mean score over all rates and windows \code{Default=TRUE} (the
#' warped windows of a strong gradient have a slightly different noise continuum, which otherwise favours or
#' penalises that gradient everywhere).
#' @param jump_prob Probability per window step of a jump of the rate to any value (an abrupt change: a hiatus,
#' an erosion surface, a facies change), after which the gradient restarts at zero \code{Default=0}. A jump is
#' not a gradient: without this the gradient-steered track can only reach a new rate through a ramp. Any value above
#' zero lets the track jump wherever the evidence gains more than the (log) probability: on the synthetic records
#' with three jumps \code{1e-3} halved the age errors, on CIP Case 3 (no jumps) it introduced a few spurious jumps
#' (correlation 0.97 to 0.88), so it is off by default and should be switched on for records with hiatuses or
#' abrupt facies changes.
#' @param gradient_prior Standard deviation (in steps of the gradient grid) of a zero-centred prior on the gradient,
#' applied at every window \code{Default=Inf} (none; a prior of 3 grid steps lowered the Case 3 result from 0.97 to
#' 0.55 and did not help the records with jumps).
#' @param scalogram Optional scalogram (\code{\link{analyse_wavelet}} etc.) whose \code{\link{track_astro_hmm}}
#' evidence is added as complementary, local evidence (same \code{astro_periods} in \code{scalogram_args}): the
#' windows of the chirp-matched spectra need 15-25 cycles and straddle abrupt changes of the rate, the short wavelets
#' see them; the sum has the sharp lines of the one and the locality of the other.
#' @param scalogram_weight Weight of that evidence relative to the template scores (both scaled by their standard
#' deviation) \code{Default=1}.
#' @param scalogram_args Arguments of \code{\link{track_astro_hmm}} for the scalogram evidence.
#' @param n_samples Posterior paths for the credible band \code{Default=100}.
#' @param level Credible level \code{Default=0.95}.
#' @param seed Seed \code{Default=1}.
#' @return A \code{wr_track} data frame (depth, sedrate, gradient, sedrate_lower, sedrate_upper, confidence10, time,
#' time_lower, time_upper) with the attribute \code{posterior} (depth, sedrate grid, posterior over the rate) and
#' \code{gradient_grid}.
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 3000, by = 1)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 19) + 0.8 * sin(2 * pi * t / 41)
#' sr <- 2 + 0.7 * sin(2 * pi * t / 600)
#' zt <- c(0, cumsum(sr[-1] * diff(t)))
#' z <- seq(0, max(zt), by = 1)
#' x <- approx(zt, x_t, z)$y + rnorm(length(z), sd = 1.5)
#' ce <- analyse_chirp_eha(data.frame(z, x), window = 250, step = 10, chirp = seq(-0.004,
#'     0.004, by = 0.0005),
#'                         lowerPeriod = 10, upperPeriod = 150, criterion = "template",
#'                         astro_periods = c(41, 23.7, 19), sedrate = c(0.8, 5))
#' tr <- track_chirp_hmm(ce)
#' plot(tr$depth, tr$sedrate, type = "l"); lines(z, approx(zt, sr, z)$y, col = 2)
#' }
#' @export
track_chirp_hmm <- function(ce, weight = 10, sd_rate = 0.01, sd_gradient = 1.5, jump_prob = 0, gradient_prior = Inf, scalogram = NULL,
                            scalogram_weight = 1, scalogram_args = list(), n_samples = 100,
                            level = 0.95, seed = 1, debias = TRUE, chirp = seq(-0.05, 0.05, by = 0.005), ...) {
  if (inherits(ce, "analyze.wavelet") && is.null(ce$template_score_all)) {
    # any scalogram: the evidence of track_astro_hmm() (rates x positions) is used with the gradient as a hidden
    # state that steers the rate (the evidence itself does not depend on the gradient)
    old <- getOption("WaverideR.track_evidence"); options(WaverideR.track_evidence = TRUE)
    E0 <- track_astro_hmm(ce, ...)
    options(WaverideR.track_evidence = old)
    nr0 <- nrow(E0$logE); nt0 <- ncol(E0$logE)
    TS <- array(E0$logE / max(1e-9, stats::sd(E0$logE)), c(nr0, 1, nt0))[, rep(1, length(chirp)), , drop = FALSE]
    ce <- list(template_score_all = TS, sedrate_grid = 2^E0$ls, chirp_grid = chirp, positions = E0$zc,
               dt = stats::median(diff(E0$zc)))
  }
  TS <- ce$template_score_all
  if (is.null(TS)) stop("ce must come from analyse_chirp_eha(criterion = 'template') or analyse_chirp_wavelet()")
  lr <- log2(ce$sedrate_grid); dj <- lr[2] - lr[1]
  g <- ce$chirp_grid; dg <- if (length(g) > 1) g[2] - g[1] else 1
  nr <- length(lr); ng <- length(g); nt <- dim(TS)[3]
  step <- ce$dt
  if (!is.null(scalogram)) {
    # complementary local evidence from a scalogram (track_astro_hmm's evidence: short wavelets see an abrupt change
    # of the rate where a window straddles it), interpolated to the rates and positions of the windows and added to
    # every gradient layer
    old <- getOption("WaverideR.track_evidence"); options(WaverideR.track_evidence = TRUE)
    E0 <- do.call(track_astro_hmm, c(list(scalogram, sedrate = range(ce$sedrate_grid)), scalogram_args))
    options(WaverideR.track_evidence = old)
    zc <- if (!is.null(ce$positions)) ce$positions else ce$axis.1
    Ez <- t(apply(E0$logE, 1, function(v) stats::approx(E0$zc, v, zc, rule = 2)$y))
    Er <- apply(Ez, 2, function(v) stats::approx(E0$ls, v, lr, rule = 2)$y)
    Er <- Er / max(1e-9, stats::sd(Er)) * stats::sd(TS, na.rm = TRUE) * scalogram_weight
    for (gi in seq_len(ng)) TS[, gi, ] <- TS[, gi, ] + Er
  }
  E <- TS; E[!is.finite(E)] <- min(E, na.rm = TRUE)
  if (debias && ng > 1) {
    # the de-chirped spectra of a warped window are not exactly comparable between gradients (the warp changes the
    # continuum of noise); the mean score of every gradient over all rates and windows is removed so that no
    # gradient is favoured a priori
    for (gi in seq_len(ng)) E[, gi, ] <- E[, gi, ] - mean(E[, gi, ])
  }
  logE <- weight * E
  if (is.finite(gradient_prior) && gradient_prior > 0 && ng > 1) {
    # a stationary prior on the gradient: zero-centred, so that the track only takes a slope when the evidence
    # asks for it (without it, sliding along the aliases of the template through the gradient states was cheap)
    gpen <- -(seq_len(ng) - which.min(abs(g)))^2 / (2 * gradient_prior^2)
    for (gi in seq_len(ng)) logE[, gi, ] <- logE[, gi, ] + gpen[gi]
  }
  # transition in rate for every gradient: shift in bins (log2 r' = log2 r + g step / ln 2), fractional, plus a
  # Gaussian spread sd_rate / dj bins
  shift <- g * step / log(2) / dj
  sr <- max(sd_rate / dj, 0.3)
  kr <- max(1, ceiling(3 * sr)); offs <- -kr:kr
  # gradient kernel
  sg <- max(sd_gradient, 0.3); kg <- max(1, ceiling(3 * sg)); goffs <- -kg:kg
  Kg <- exp(-goffs^2 / (2 * sg^2)); Kg <- Kg / sum(Kg)
  move <- function(V, gi) {
    # V: vector over rates (log prob) at gradient gi; returns the best-predecessor max over the rate move
    out <- rep(-Inf, nr); arg <- rep(NA_integer_, nr)
    for (o in offs) {
      d <- shift[gi] + o
      lo <- floor(d); fr <- d - lo
      w <- -o^2 / (2 * sr^2)
      for (part in list(c(lo, 1 - fr), c(lo + 1, fr))) {
        sft <- part[1]; pw <- part[2]
        if (pw <= 0) next
        src <- seq_len(nr) - sft
        ok <- src >= 1 & src <= nr
        cand <- rep(-Inf, nr); cand[ok] <- V[src[ok]] + w + log(pw)
        better <- cand > out
        out[better] <- cand[better]; arg[better] <- src[better]
      }
    }
    list(best = out, arg = arg)
  }
  # Viterbi over (rate, gradient); a jump goes from the best previous state to any rate with gradient zero
  delta <- matrix(logE[, , 1], nr, ng)
  ptr_r <- vector("list", nt); ptr_g <- vector("list", nt); jump_from <- vector("list", nt)
  lKg <- log(Kg)
  g0 <- which.min(abs(g)); ljump <- if (jump_prob > 0) log(jump_prob) - log(nr) else -Inf
  lstay <- if (jump_prob > 0) log1p(-jump_prob) else 0
  for (t in 2:nt) {
    # step 1: gradient drift (max over previous gradient)
    D1 <- matrix(-Inf, nr, ng); G1 <- matrix(NA_integer_, nr, ng)
    for (gi in seq_len(ng)) for (m in seq_along(goffs)) {
      gp <- gi - goffs[m]
      if (gp < 1 || gp > ng) next
      cand <- delta[, gp] + lKg[m]
      better <- cand > D1[, gi]
      D1[better, gi] <- cand[better]; G1[better, gi] <- gp
    }
    # step 2: rate move by the (new) gradient
    D2 <- matrix(-Inf, nr, ng); R2 <- matrix(NA_integer_, nr, ng)
    for (gi in seq_len(ng)) { mv <- move(D1[, gi], gi); D2[, gi] <- mv$best + lstay; R2[, gi] <- mv$arg }
    if (is.finite(ljump)) {
      bi <- which(delta == max(delta), arr.ind = TRUE)[1, ]
      jv <- delta[bi[1], bi[2]] + ljump
      better <- jv > D2[, g0]
      if (any(better)) { D2[better, g0] <- jv; R2[better, g0] <- -bi[1]; jump_from[[t]] <- bi }
    }
    delta <- D2 + matrix(logE[, , t], nr, ng)
    ptr_r[[t]] <- R2; ptr_g[[t]] <- G1
  }
  path_r <- integer(nt); path_g <- integer(nt)
  i <- which(delta == max(delta), arr.ind = TRUE)[1, ]
  path_r[nt] <- i[1]; path_g[nt] <- i[2]
  for (t in nt:2) {
    r <- path_r[t]; gi <- path_g[t]
    rp <- ptr_r[[t]][r, gi]
    if (rp < 0) { bi <- jump_from[[t]]; path_r[t - 1] <- bi[1]; path_g[t - 1] <- bi[2]; next }
    gp <- ptr_g[[t]][rp, gi]
    path_r[t - 1] <- rp; path_g[t - 1] <- gp
  }
  # posterior over the rate by forward-backward in log space (gradient marginalised), for the band and confidence
  lse <- function(a) { m <- max(a); if (!is.finite(m)) return(-Inf); m + log(sum(exp(a - m))) }
  fwd <- array(-Inf, c(nr, ng, nt)); fwd[, , 1] <- logE[, , 1]
  logE_t <- function(t) matrix(logE[, , t], nr, ng)
  spread <- function(A) {
    # A: nr x ng log probabilities; apply gradient drift then rate move (sum instead of max)
    B <- matrix(-Inf, nr, ng)
    for (gi in seq_len(ng)) {
      # log-sum over the gradient kernel
      vals <- sapply(seq_along(goffs), function(m) { gp <- gi - goffs[m]; if (gp < 1 || gp > ng) rep(-Inf, nr) else A[, gp] + lKg[m] })
      acc <- apply(matrix(vals, nr), 1, lse)
      # rate move: sum over the offsets
      tot <- rep(-Inf, nr)
      for (o in offs) {
        d <- shift[gi] + o; lo <- floor(d); fr <- d - lo; w <- -o^2 / (2 * sr^2)
        for (part in list(c(lo, 1 - fr), c(lo + 1, fr))) {
          sft <- part[1]; pw <- part[2]; if (pw <= 0) next
          src <- seq_len(nr) - sft; ok <- src >= 1 & src <= nr
          cand <- rep(-Inf, nr); cand[ok] <- acc[src[ok]] + w + log(pw)
          m <- pmax(tot, cand)
          tot <- ifelse(is.finite(m), m + log1p(exp(-abs(tot - cand))), m)
        }
      }
      B[, gi] <- tot + lstay
    }
    if (is.finite(ljump)) {
      jv <- lse(as.vector(A)) + ljump
      m <- pmax(B[, g0], jv); B[, g0] <- ifelse(is.finite(m), m + log1p(exp(-abs(B[, g0] - jv))), m)
    }
    B
  }
  for (t in 2:nt) {
    f <- spread(matrix(fwd[, , t - 1], nr, ng)) + logE_t(t)
    fwd[, , t] <- f - max(f)                                 # normalised per window (avoids overflow)
  }
  post <- matrix(0, nr, nt)
  for (t in seq_len(nt)) { v <- apply(matrix(fwd[, , t], nr, ng), 1, lse); v <- exp(v - max(v)); post[, t] <- v / sum(v) }
  # (filtering posterior; the backward pass is omitted for speed: the band is from the filtered posterior)
  rate <- ce$sedrate_grid[path_r]
  grad <- g[path_g]
  z <- if (!is.null(ce$positions)) ce$positions else ce$axis.1
  lo_i <- apply(post, 2, function(p) which(cumsum(p) >= (1 - level) / 2)[1])
  hi_i <- apply(post, 2, function(p) which(cumsum(p) >= 1 - (1 - level) / 2)[1])
  conf10 <- sapply(seq_len(nt), function(t) sum(post[abs(lr - lr[path_r[t]]) <= log2(1.1), t]))
  tm <- c(0, cumsum(diff(z) / rate[-1]))
  out <- data.frame(depth = z, sedrate = rate, gradient = grad, sedrate_lower = ce$sedrate_grid[lo_i],
                    sedrate_upper = ce$sedrate_grid[hi_i], confidence10 = conf10, time = tm,
                    time_lower = c(0, cumsum(diff(z) / ce$sedrate_grid[hi_i][-1])),
                    time_upper = c(0, cumsum(diff(z) / ce$sedrate_grid[lo_i][-1])))
  attr(out, "posterior") <- list(depth = z, sedrate = ce$sedrate_grid, post = post, evidence = apply(logE, c(1, 3), max))
  attr(out, "gradient_grid") <- g
  class(out) <- c("wr_track", "data.frame")
  out
}
