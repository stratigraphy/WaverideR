#' @title Wavelet synchrosqueezing transform
#'
#' @description Computes the continuous wavelet transform of a record and sharpens it with
#' synchrosqueezing (Daubechies, Lu and Wu 2011; Thakur et al. 2013). A wavelet scalogram smears the
#' energy of a cycle over neighbouring periods because the wavelet has a finite bandwidth. Synchrosqueezing
#' estimates the instantaneous frequency of every wavelet coefficient from the phase of the transform
#' (the time derivative of the coefficient divided by the coefficient) and moves the coefficient to
#' that frequency. The result is a time-frequency representation in which a cycle appears as a thin ridge
#' at its instantaneous period, so that closely spaced or slowly changing cycles are better separated
#' than in the scalogram, and the transform stays invertible: summing the squeezed coefficients gives
#' back the record, so a cycle can be reconstructed from a band of the squeezed transform.
#'
#' The result has the same structure as the output of \code{\link{analyse_wavelet}} (class
#' \code{analyze.wavelet}), so it can be plotted with \code{\link{plot_wavelet}} and used by the tracking and extraction
#' functions (e.g. \code{\link{track_period_wavelet}}, \code{\link{extract_signal}}) in place of the scalogram.
#'
#' @param data Input is a time series with the first column being depth or time and the second column being a proxy.
#' @param dj Spacing of the scales in octaves (as in \code{\link{analyse_wavelet}}); also the spacing of the squeezed
#' period axis.
#' @param lowerPeriod,upperPeriod Shortest and longest period of the transform.
#' @param omega_nr Central angular frequency of the Morlet wavelet (see \code{\link{analyse_wavelet}}); several values (e.g. \code{c(6, 10, 16, 24)}) select adaptive synchrosqueezing, see \code{adaptive_window}. Synchrosqueezing separates cycles that are further apart than the bandwidth of the wavelet, so a larger value (10 to 12) separates close cycles better; the interference of two close cycles otherwise makes the ridges wave.
#' @param wavelet,wavelet_order Mother wavelet and its order, see \code{\link{analyse_wavelet}}. Synchrosqueezing needs an
#' analytic wavelet of \code{\link{analyse_wavelet}} (not the real wavelets dog, morl, hhat, poisson, haar, db, sym and battle).
#' @param n_iter Number of synchrosqueezing iterations (multisynchrosqueezing, Yu et al. 2019) \code{Default=1}.
#' @param threshold Coefficients with a modulus below this value are not squeezed, because the phase of a
#' coefficient that is dominated by noise gives a meaningless frequency, and squeezing them scatters noise
#' over the whole spectrum. \code{"riskshrink"} (default) estimates the noise level from the finest scales as
#' \code{sqrt(2 log(N)) * 1.4826 * MAD(|W|)} (Donoho and Johnstone 1994; as in the Synchrosqueezing Toolbox of Brevdo
#' and Thakur and in \code{ssqueezepy}); a number is used directly, in standard deviations of the (standardised) record,
#' since the transform is normalised so that white noise gives the same modulus at every scale; \code{"none"}
#' uses a negligible value (sqrt of the machine epsilon), which squeezes everything.
#' @param sst_order \code{1} (default) for the classic synchrosqueezing, \code{2} for the second-order
#' synchrosqueezing transform (Oberlin, Meignen and Perrier 2015; Oberlin and Meignen 2017). The first-order
#' estimate of the instantaneous frequency is the frequency at the centre of gravity (in time) of the part of the
#' signal that a wavelet sees, which lags behind when the period changes quickly (a changing sedimentation rate makes
#' the cycle a chirp). The second-order estimate measures the rate of change of the frequency (the chirp rate) from
#' additional transforms and corrects for it, so that ridges follow fast changes of the sedimentation rate without
#' giving up frequency resolution.
#' @param type \code{"sst"} (default) squeezes every coefficient to its instantaneous period; \code{"set"} is the
#' synchroextracting transform (Yu et al. 2017), which keeps only the coefficients whose instantaneous period equals
#' the period of their own scale (the ridge itself) and gives the thinnest ridges, at the cost of discarding the energy
#' away from the ridge (it cannot be used to reconstruct the full signal).
#' @param adaptive_window Only used when \code{omega_nr} has several values (adaptive synchrosqueezing): the transform is
#' computed for every \code{omega_nr} and at every position the one with the most concentrated energy (lowest Renyi
#' entropy of order 3 within a moving window of this length, in the units of the depth/time axis) is used. A short
#' wavelet is then chosen where the period changes quickly and a long wavelet where close cycles need to be separated.
#' \code{Default=NULL}: three times the geometric mean of \code{lowerPeriod} and \code{upperPeriod}.
#' @param adaptive_on Transform on which the entropy of the adaptive selection is measured: \code{"sst"} (default),
#' the squeezed transform, which picks short wavelets where the period changes quickly, but in very noisy records also
#' prefers short wavelets (noise is squeezed into sparse spikes); or \code{"cwt"}, the wavelet transform, which is
#' robust against noise but nearly always prefers the longest wavelet.
#' @param verbose Print messages.
#'
#' @details Every coefficient W(s, t) is moved to the period that its instantaneous frequency gives. With
#' the Torrence and Compo (1998) normalisation the reconstruction sum W / sqrt(s) is conserved: the squeezed
#' coefficients are stored such that \code{Re(Wave) / sqrt(Scale)} summed over the periods equals the same sum of
#' the scalogram, which is what \code{\link{extract_signal}} and related functions use. Coefficients whose instantaneous
#' period falls outside \code{lowerPeriod} - \code{upperPeriod} are discarded, so very strong cycles at the edge of
#' the range lose part of their energy there.
#'
#' The ridge of a stationary cycle is sharp; for a cycle with a changing period the ridge follows the
#' period (the instantaneous frequency). The power of the squeezed transform is not proportional to the
#' squared amplitude of the cycle in the same way as in a scalogram, because energy is concentrated.
#'
#' @return A list of class \code{analyze.wavelet} like \code{\link{analyse_wavelet}} (\code{Wave}, \code{Power},
#' \code{Ampl}, \code{Phase}, \code{Period}, \code{Scale}, \code{coi.1}, \code{coi.2}, \code{axis.1}, \code{axis.2}, ...)
#' for the squeezed transform, with the extra elements \code{synchrosqueezed = TRUE}, \code{Wave_cwt} and
#' \code{Power_cwt} (the unsqueezed transform) and \code{inst_period} (matrix with the instantaneous period of every coefficient).
#'
#' @references
#' Daubechies, I., Lu, J., Wu, H.-T., 2011. Synchrosqueezed wavelet transforms: an empirical mode
#' decomposition-like tool. Applied and Computational Harmonic Analysis 30, 243-261.
#'
#' Thakur, G., Brevdo, E., Fuckar, N.S., Wu, H.-T., 2013. The synchrosqueezing algorithm for time-varying
#' spectral analysis: robustness properties and new paleoclimate applications. Signal Processing 93, 1079-1094.
#'
#' Torrence, C., Compo, G.P., 1998. A practical guide to wavelet analysis. Bulletin of the American
#' Meteorological Society 79, 61-78.
#'
#' Donoho, D.L., Johnstone, I.M., 1994. Ideal spatial adaptation by wavelet shrinkage. Biometrika 81, 425-455.
#'
#' Oberlin, T., Meignen, S., Perrier, V., 2015. Second-order synchrosqueezing transform or invertible reassignment?
#' Towards ideal time-frequency representations. IEEE Transactions on Signal Processing 63, 1335-1344.
#'
#' Oberlin, T., Meignen, S., 2017. The second-order wavelet synchrosqueezing transform. 2017 IEEE International
#' Conference on Acoustics, Speech and Signal Processing (ICASSP), 3994-3998.
#'
#' Yu, G., Wang, Z., Zhao, P., 2019. Multisynchrosqueezing transform. IEEE Transactions on Industrial Electronics 66,
#' 5441-5455.
#'
#' Auger, F., Flandrin, P., 1995. Improving the readability of time-frequency and time-scale representations by the
#' reassignment method. IEEE Transactions on Signal Processing 43, 1068-1089.
#'
#' Yu, G., Yu, M., Xu, C., 2017. Synchroextracting transform. IEEE Transactions on Industrial Electronics 64, 8042-8054.
#'
#' Sheu, Y.-L., Hsu, L.-Y., Chou, P.-T., Wu, H.-T., 2017. Entropy-based time-varying window width selection for
#' nonlinear-type time-frequency analysis. International Journal of Data Science and Analytics 3, 231-245.
#'
#' @examples
#' t <- seq(0, 1000, by = 2)
#' x <- sin(2 * pi * t / 45) + sin(2 * pi * t / 60) + rnorm(length(t), sd = 0.2)
#' cwt <- analyse_wavelet(data.frame(t, x), lowerPeriod = 20, upperPeriod = 150, dj = 1 / 40)
#' sst <- analyse_wavelet_sst(data.frame(t, x), lowerPeriod = 20, upperPeriod = 150, dj = 1 / 40)
#' op <- par(mfrow = c(1, 2))
#' plot_wavelet(cwt)
#' plot_wavelet(sst)
#' par(op)
#' @importFrom stats approx sd
#' @export
analyse_wavelet_sst <- function(data = NULL,
                                dj = 1 / 100,
                                lowerPeriod = 2,
                                upperPeriod = 1024,
                                omega_nr = 12,
                                wavelet = "morlet",
                                wavelet_order = NULL,
                                threshold = "riskshrink",
                                sst_order = 1,
                                type = c("sst", "set", "reassign"),
                                n_iter = 1,
                                adaptive_window = NULL,
                                adaptive_on = c("sst", "cwt"),
                                verbose = FALSE) {
  settings <- mget(setdiff(names(formals(sys.function())), c("data", "verbose")))
  type <- match.arg(type)
  adaptive_on <- match.arg(adaptive_on)
  if (!(sst_order %in% c(1, 2))) stop("sst_order must be 1 or 2")
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dt <- stats::median(diff(dat[, 1]))
  xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
  xi <- stats::approx(dat[, 1], dat[, 2], xout)$y
  if (verbose){
    cat("dataset interpolated to: ", round(dt, 10), "\n")
  }
  n <- length(xi)
  xs <- (xi - mean(xi)) / stats::sd(xi)

  runs <- lapply(omega_nr, function(om) {
    wdef <- .wr_wavelet(wavelet, wavelet_order, om)
    if (!wdef$analytic){
      stop("synchrosqueezing needs an analytic wavelet: use an analytic wavelet such as \"morlet\", \"morlet0\", ",
           "\"paul\", \"morse\", \"bump\", \"cmor\", \"cgau\", \"loggabor\" or \"gammatone\"; real wavelets ",
           "(dog, morl, hhat, poisson, haar, db, sym, battle) cannot be synchrosqueezed")
    }
    # common period axis for all omega_nr: the scales follow from the Fourier factor of each wavelet
    J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
    periods <- lowerPeriod * 2^((0:J) * dj)
    scales <- periods / wdef$fourier_factor
    r <- .wr_sst_core(xs, dt, wdef, scales, dj, sst_order, type, threshold, n_iter)
    r$wdef <- wdef
    r$scales <- scales
    r$periods <- periods
    r
  })
  periods <- runs[[1]]$periods
  ns <- length(periods)
  sel <- rep(1L, n)
  if (length(runs) > 1) {
    # adaptive: per position the omega_nr with the lowest Renyi entropy (order 3) in a moving window
    if (is.null(adaptive_window)) adaptive_window <- 3 * sqrt(lowerPeriod * upperPeriod)
    L <- max(3, 2 * floor(adaptive_window / dt / 2) + 1)
    L <- min(L, 2 * floor((n - 1) / 2) + 1)
    k <- rep(1, L)
    H <- sapply(runs, function(r) {
      P <- if (adaptive_on == "cwt") Mod(r$wave)^2 / r$scales else Mod(r$tx)^2 / r$scales
      s1 <- as.numeric(stats::filter(colSums(P), k, sides = 2))
      s3 <- as.numeric(stats::filter(colSums(P^3), k, sides = 2))
      h <- -0.5 * log(s3 / s1^3)
      h[!is.finite(h)] <- NA
      h
    })
    H <- apply(H, 2, function(h) {   # fill the edges of the window with the nearest value
      ok <- which(!is.na(h))
      if (!length(ok)) return(rep(0, n))
      stats::approx(ok, h[ok], seq_len(n), rule = 2)$y
    })
    sel <- max.col(-H, ties.method = "first")
  }
  pick <- function(el) {
    out <- runs[[1]][[el]]
    if (length(runs) > 1) for (i in 2:length(runs)) {
      j <- which(sel == i)
      if (length(j)) out[, j] <- runs[[i]][[el]][, j]
    }
    out
  }
  tx <- pick("tx")
  wave <- pick("wave")
  p_inst <- pick("p_inst")
  scales <- runs[[1]]$scales
  wdef <- runs[[1]]$wdef
  # the squeezed coefficients are stored on the scales of the first wavelet (the sum of W / sqrt(s) is kept)
  if (length(runs) > 1) for (i in 2:length(runs)) {
    j <- which(sel == i)
    if (length(j)) {
      tx[, j] <- tx[, j] * sqrt(scales / runs[[i]]$scales)
      wave[, j] <- wave[, j] * sqrt(scales / runs[[i]]$scales)
    }
  }
  wave_sst <- tx
  Power <- Mod(wave_sst)^2 / matrix(rep(scales, n), nrow = ns)
  Ampl <- Mod(wave_sst) / matrix(rep(sqrt(scales), n), nrow = ns)
  fourier.factor <- wdef$fourier_factor
  lp <- log2(periods)

  x_axis <- xout
  coi <- fourier.factor * .wr_coi_factor(wdef) * dt *
    c(1e-5, 1:((n + 1) / 2 - 1), rev(1:(n / 2 - 1)), 1e-5)
  coi.x <- c(x_axis[c(1, 1)] - dt * 0.5, x_axis, x_axis[c(n, n)] + dt * 0.5)
  logyint <- lp[2] - lp[1]
  yl <- c(lp[ns] + 0.5 * logyint, lp[1] - 0.5 * logyint)
  coi.y <- c(yl, log2(coi), rev(yl))
  Power_cwt <- Mod(wave)^2 / matrix(rep(scales, n), nrow = ns)

  out <- list(Wave = wave_sst, Phase = Arg(wave_sst), Ampl = Ampl, Power = Power,
              dt = dt, dj = dj, Power.avg = rowMeans(Power), Period = periods, Scale = scales,
              coi.1 = coi.x, coi.2 = coi.y, nc = n, nr = ns, axis.1 = x_axis, axis.2 = lp,
              omega_nr = omega_nr[1], wavelet = wdef$name, wavelet_order = wdef$order,
              fourier_factor = fourier.factor, psi0 = wdef$psi0,
              x = xout, y = xi, Power.avg.pval = NULL, Power.pval = NULL,
              synchrosqueezed = TRUE, sst_order = sst_order, type = type,
              threshold = sapply(runs, function(r) r$gam),
              omega_nr_all = omega_nr, omega_selected = omega_nr[sel],
              Wave_cwt = wave, Power_cwt = Power_cwt, inst_period = p_inst)
  out$settings <- settings
  out$transform <- "analyse_wavelet_sst"
  class(out) <- "analyze.wavelet"
  invisible(out)
}

# Core of the (second-order) synchrosqueezing / synchroextracting transform of a standardised,
# regularly sampled record xs. Returns the squeezed coefficients tx on the given scales (stored such
# that the sum over scales of Re(tx) / sqrt(s) equals that of the wavelet transform), the wavelet
# transform, the instantaneous periods and the threshold.
.wr_sst_core <- function(xs, dt, wdef, scales, dj, sst_order = 1, type = "sst", threshold = "riskshrink",
                         n_iter = 1) {
  n <- length(xs)
  ns <- length(scales)
  pad <- 2^(trunc(log2(n) + 0.5) + 1) - n
  N <- n + pad
  omega.k <- 1:floor(N / 2) * (2 * pi) / (N * dt)
  omega.k <- c(0, omega.k, -omega.k[floor((N - 1) / 2):1])
  fx <- fft(c(xs, rep(0, pad)))
  inv <- function(v) (fft(fx * v, inverse = TRUE) / N)[1:n]
  periods <- scales * wdef$fourier_factor
  lp <- log2(periods)

  wave <- matrix(0 + 0i, ns, n)
  for (j in seq_len(ns)) wave[j, ] <- inv(wdef$daughter(scales[j], omega.k, dt))
  amod <- Mod(wave)
  if (is.character(threshold)){
    threshold <- match.arg(threshold, c("riskshrink", "none"))
    if (threshold == "riskshrink"){
      nfine <- min(ns, max(1, round(1 / dj)))
      fine <- as.numeric(amod[seq_len(nfine), ])
      gam <- sqrt(2 * log(n)) * 1.4826 * stats::mad(fine, constant = 1)
      # a nearly noise-free series makes the finest scales carry signal rather than noise, and the universal
      # threshold would then discard every coefficient: keep at least a tenth of them
      gam <- min(gam, stats::quantile(amod, 0.9, names = FALSE))
    } else {
      gam <- sqrt(.Machine$double.eps)
    }
  } else {
    gam <- threshold
  }
  tx <- matrix(0 + 0i, ns, n)
  p_inst <- matrix(NA_real_, ns, n)
  bmat <- matrix(NA_integer_, ns, n)
  if (type == "reassign") {
    tbin <- matrix(NA_integer_, ns, n)
    reas <- matrix(0, ns, n)
  }
  for (j in seq_len(ns)) {
    dau <- wdef$daughter(scales[j], omega.k, dt)
    w0 <- wave[j, ]
    valid <- amod[j, ] > gam
    wv <- ifelse(valid, w0, NA)
    w1 <- inv(dau * 1i * omega.k)
    om <- Im(w1 / wv)
    om1 <- om
    if (sst_order == 2 || type == "reassign") {
      # derivative of the daughter with respect to omega: transform with the time weighted wavelet
      h <- 1e-4 / scales[j]
      ddau <- (wdef$daughter(scales[j], omega.k + h, dt) - wdef$daughter(scales[j], omega.k - h, dt)) / (2 * h)
      w2 <- inv(dau * -(omega.k^2))
      wt <- inv(1i * ddau)
      wt1 <- inv(1i * ddau * 1i * omega.k)
      r1 <- w1 / wv
      d_om <- Im(w2 / wv - r1^2)                       # d omega_hat / d b
      tau <- Re(wt / wv)                               # b - t_hat
      d_t <- 1 - Re(wt1 / wv - wt * w1 / wv^2)          # d t_hat / d b
      q <- ifelse(abs(d_t) > 1e-6, d_om / d_t, 0)
      om2 <- om + q * tau
      # use the corrected estimate where it is defined and positive
      if (sst_order == 2) om <- ifelse(is.finite(om2) & om2 > 0, om2, om)
      if (type == "reassign") {
        # reassignment in depth/time as well: the centre of gravity t_hat = b - Re(W_tau / W)
        tb <- round((seq_len(n) - 1) - tau / dt) + 1
        tb[!is.finite(tb) | tb < 1 | tb > n] <- NA
        tbin[j, ] <- tb
      }
    }
    pj <- ifelse(is.finite(om) & om > 0, 2 * pi / om, NA_real_)
    p_inst[j, ] <- pj
    b <- round((log2(pj) - lp[1]) / dj) + 1
    b[!is.finite(b) | b < 1 | b > ns] <- NA
    if (type == "set") {
      # synchroextracting: keep only the coefficients on the ridge, i.e. those whose (first-order)
      # instantaneous period lies within half a bin of the period of their own scale
      p1 <- ifelse(is.finite(om1) & om1 > 0, 2 * pi / om1, NA_real_)
      on_ridge <- !is.na(p1) & abs(log2(p1) - lp[j]) <= dj / 2
      b <- ifelse(on_ridge, j, NA)
    }
    bmat[j, ] <- b
  }
  if (type == "sst" && n_iter > 1) {
    # multisynchrosqueezing (Yu et al. 2019): apply the frequency reassignment repeatedly, i.e. move a
    # coefficient to the instantaneous period estimated at the period it was moved to
    cols <- rep(seq_len(n), each = ns)
    for (it in 2:n_iter) {
      cur <- as.vector(bmat)
      nxt <- cur
      okk <- !is.na(cur)
      nxt[okk] <- bmat[cbind(cur[okk], cols[okk])]
      nxt[okk & is.na(nxt)] <- cur[okk & is.na(nxt)]
      bmat <- matrix(nxt, ns, n)
    }
  }
  for (j in seq_len(ns)) {
    b <- bmat[j, ]
    ok <- which(!is.na(b))
    if (!length(ok)) next
    if (type == "reassign") {
      tb <- tbin[j, ok]
      k <- !is.na(tb)
      if (!any(k)) next
      idx <- cbind(b[ok][k], tb[k])
      # accumulate power (|W|^2 / s): the reassigned scalogram is not invertible
      v <- Mod(wave[j, ok][k])^2 / scales[j]
      agg <- stats::aggregate(v, by = list(idx[, 1], idx[, 2]), FUN = sum)
      reas[cbind(agg[[1]], agg[[2]])] <- reas[cbind(agg[[1]], agg[[2]])] + agg[[3]]
    } else {
      idx <- cbind(b[ok], ok)
      tx[idx] <- tx[idx] + wave[j, ok] / sqrt(scales[j])
    }
  }
  if (type == "reassign") {
    # store as coefficients with the power of the reassigned scalogram (phase zero)
    tx <- sqrt(reas) + 0i
    return(list(tx = tx * sqrt(scales), wave = wave, p_inst = p_inst, gam = gam))
  }
  list(tx = tx * sqrt(scales), wave = wave, p_inst = p_inst, gam = gam)
}
