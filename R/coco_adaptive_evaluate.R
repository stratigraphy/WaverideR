# Core evaluator (cocoAdaptiveEvaluate.m)
#' @title Evaluate Adaptive COCO for a batch of spectra (cocoAdaptiveEvaluate.m)
#'
#' @param data_power One-sided periodograms of the data to evaluate (numeric matrix, one spectrum per column, with
#' \code{floor(max(pad, n_target)/2) + 1} rows on the frequency grid \code{data_frequency}); negative values are set
#' to zero. Required, no default.
#' @param data_frequency Frequencies of the rows of \code{data_power}, in cycles per metre (numeric vector, evenly
#' spaced from 0). Required, no default.
#' @param n_target Number of samples of the analysed record or window (integer); the target sinusoids are
#' synthesised on this many points. Required, no default.
#' @param spacing Depth sampling interval of the analysed record, in metres (numeric, > 0). Required, no default.
#' @param pad Length (number of points) to which the record was zero-padded for the periodogram (integer); the
#' target spectra are padded to \code{max(pad, n_target)}. Required, no default.
#' @param orbit9 Astronomical target periods in kyr (numeric vector; nine distinct values are needed for
#' \code{amplitude_mode="four-group-area"} and for the grouped weights of \code{"fixed"}). Required, no default.
#' @param sr_grid Sedimentation rates to evaluate, in cm/kyr (numeric vector, > 0). Required, no default.
#' @param method Correlation coefficient between the data and target spectra (character): \code{"Pearson"} or
#' \code{"Spearman"}. \code{Default="Pearson"}.
#' @param max_frequency Highest frequency of the target that enters the correlation, in cycles/kyr (numeric, must be
#' at least \code{max(1/orbit9)}); \code{NULL} uses \code{1.2 * max(1/orbit9)}. \code{Default=NULL}.
#' @param target_model How the target spectrum is synthesised (character): \code{"coherent-nine"} (spectrum of the
#' sum of the sinusoids, with their interference) or \code{"phase-averaged"} (sum of the individual spectra,
#' averaged over sine and cosine phase). \code{Default=c("coherent-nine", "phase-averaged")}, i.e.
#' \code{"coherent-nine"}.
#' @param amplitude_mode How the target amplitudes are set (character): \code{"four-group-area"} (one amplitude
#' per group of long eccentricity, short eccentricity, obliquity and precession, fitted to the band areas; Acycle
#' default for Adaptive COCO/eCOCO; needs \code{target_model="coherent-nine"}), \code{"per-orbit"} (each amplitude
#' fitted to the spectral peak at its own frequency, Acycle adaptive9a) or \code{"fixed"} (weights 1.0/0.8/0.6 from
#' \code{\link{coco_fixed_target_weights}}). \code{Default=c("four-group-area", "per-orbit", "fixed")}, i.e.
#' \code{"four-group-area"}.
#' @param observed_p Also compute the parametric p-value of the correlation for the first column of
#' \code{data_power} at every rate (logical, TRUE or FALSE); used for the observed spectrum, not for the Monte Carlo
#' spectra. \code{FALSE} returns \code{NA} p-values. \code{Default=FALSE}.
#' @export
coco_adaptive_evaluate <- function(data_power, data_frequency, n_target, spacing,
    pad, orbit9, sr_grid,
    method = "Pearson", max_frequency = NULL,
    target_model = c("coherent-nine", "phase-averaged"),
    amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
    observed_p = FALSE){
  # Internal helper functions (local to this function)

  detrend_linear <- function(x){
    # Least-squares removal of a straight line from every column (MATLAB detrend(x,1)).
    xm <- as.matrix(x)
    n <- nrow(xm)
    if (n < 2){
      if (is.matrix(x)){
        return(xm * 0)
      } else {
        return(as.vector(xm) * 0)
      }
    }
    t <- seq_len(n)
    X <- cbind(1, t)
    q <- qr(X)
    res <- qr.resid(q, xm)
    if (is.matrix(x)){
      res
    } else {
      as.vector(res)
    }
  }

  # Correlation helpers
  cor_with_p <- function(x, y, method = c("Pearson", "Spearman")){
    method <- match.arg(method)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    n <- length(x)
    if (n < 3){
      return(c(rho = NA_real_, p = NA_real_))
    }
    if (stats::sd(x) == 0 || stats::sd(y) == 0){
      return(c(rho = NA_real_, p = NA_real_))
    }
    if (method == "Pearson"){
      cm <- "pearson"
    } else {
      cm <- "spearman"
    }
    r <- stats::cor(x, y, method = cm)
    tt <- r * sqrt((n - 2) / max(1e-300, 1 - r^2))
    p <- 2 * stats::pt(-abs(tt), n - 2)
    c(rho = r, p = p)
  }

  # Physical groups of a 9-period target by descending rank (cocoOrbitGroups.m)
  coco_orbit_groups <- function(orbit9){
    if (length(orbit9) != 9 || length(unique(orbit9)) != 9){
      stop("orbit9 must contain nine distinct periods")
    }
    o <- order(orbit9, decreasing = TRUE)
    idx <- integer(9)
    idx[o[1]] <- 1L
    idx[o[2:5]] <- 2L
    idx[o[6]] <- 3L
    idx[o[7:9]] <- 4L
    idx
  }

  .column_correlation <- function(x, y, method){
    x <- as.matrix(x)
    y <- as.matrix(y)
    if (method == "Spearman"){
      xr <- x
      for (j in seq_len(ncol(x))){
        xr[, j] <- rank(x[, j])
      }
      yr <- y
      for (j in seq_len(ncol(y))){
        yr[, j] <- rank(y[, j])
      }
      x <- as.matrix(xr)
      y <- as.matrix(yr)
    }
    n <- nrow(x)
    xc <- sweep(x, 2, colMeans(x))
    yc <- sweep(y, 2, colMeans(y))
    ssx <- colSums(xc^2)
    ssy <- colSums(yc^2)
    maxx <- numeric(ncol(x))
    for (j in seq_len(ncol(x))){
      maxx[j] <- max(abs(x[, j]))
    }
    maxy <- numeric(ncol(y))
    for (j in seq_len(ncol(y))){
      maxy[j] <- max(abs(y[, j]))
    }
    tolx <- n * (256 * .eps * pmax(maxx, .Machine$double.xmin))^2
    toly <- n * (256 * .eps * pmax(maxy, .Machine$double.xmin))^2
    r <- colSums(xc * yc) / sqrt(ssx * ssy)
    r[ssx <= tolx | ssy <= toly] <- 0
    r[!is.finite(r)] <- NA
    r
  }

  .nnls4 <- function(M, E){
    # Exact nonnegative least squares min ||M q - E||, q >= 0, by enumerating
    # the 15 active sets (cocoNonnegativeLeakageSolve.m). E is 4-by-k.
    E <- as.matrix(E)
    k <- ncol(E)
    cs <- numeric(k)
    for (j in seq_len(k)){
      cs[j] <- max(E[, j])
    }
    cs[!is.finite(cs) | cs <= 0] <- 1
    S <- sweep(E, 2, cs, "/")
    G <- matrix(0, 4, k)
    best <- colSums(S^2)
    for (code in 1:15){
      act <- as.logical(bitwAnd(code, c(1L, 2L, 4L, 8L)))
      D <- M[, act, drop = FALSE]
      q <- qr(D)
      if (q$rank < sum(act)){
        next
      }
      cand <- matrix(qr.coef(q, S), nrow = sum(act))
      cmax <- numeric(k)
      for (j in seq_len(k)){
        cmax[j] <- max(abs(cand[, j]))
      }
      tol <- 256 * .eps * pmax(1, cmax)
      above <- cand >= -matrix(tol, nrow(cand), k, byrow = TRUE)
      feas <- logical(k)
      for (j in seq_len(k)){
        feas[j] <- all(above[, j])
      }
      if (!any(feas)){
        next
      }
      cand <- pmax(cand, 0)
      res <- colSums((D %*% cand - S)^2)
      better <- feas & is.finite(res) & res < best
      if (any(better)){
        G[, better] <- 0
        G[act, better] <- cand[, better]
        best[better] <- res[better]
      }
    }
    sweep(G, 2, cs, "*")
  }

  # Small helpers
  .eps <- .Machine$double.eps

  # Main code
  target_model <- match.arg(target_model)
  amplitude_mode <- match.arg(amplitude_mode)
  if (amplitude_mode == "four-group-area" && target_model != "coherent-nine"){
    stop("four-group-area requires target_model = 'coherent-nine'")
  }
  P <- as.matrix(data_power)
  P[P < 0] <- 0
  orbit9 <- as.vector(orbit9)
  if (is.null(max_frequency)){
    max_frequency <- 1.2 * max(1 / orbit9)
  }
  if (max_frequency < max(1 / orbit9) * (1 - 1e-12)){
    stop("max_frequency must include the highest orbital frequency")
  }
  if (amplitude_mode == "four-group-area"){
    group <- coco_orbit_groups(orbit9)
  } else {
    group <- NULL
  }
  if (amplitude_mode == "fixed"){
    fixed_amp <- coco_fixed_target_weights(orbit9)
  } else {
    fixed_amp <- NULL
  }
  nfft <- max(ceiling(pad), n_target)
  nf <- floor(nfft / 2) + 1
  if (nrow(P) != nf){
    stop("data_power rows must equal floor(pad/2)+1")
  }
  rbw_spatial <- 1 / (n_target * spacing)       # enbw(rectwin(n), 1/dz)
  df_spatial <- mean(diff(data_frequency))
  w1 <- rep(2, nf)
  w1[1] <- 1
  if (nfft %% 2 == 0){
    w1[nf] <- 1
  }
  nRate <- length(sr_grid)
  k <- ncol(P)
  no <- length(orbit9)
  rho <- matrix(NA_real_, nRate, k)
  pval <- rep(NA_real_, nRate)
  nMissing <- rep(NA_real_, nRate)
  center <- 1 / orbit9
  tphase0 <- (0:(n_target - 1))

  for (m in seq_len(nRate)){
    sr <- sr_grid[m]
    Fs <- sr / (100 * spacing)
    tt <- tphase0 / Fs
    phase <- 2 * pi * outer(tt, 1 / orbit9)
    freq <- (0:(nf - 1)) * Fs / nfft
    rbw_t <- Fs / n_target
    # rows needed: everything up to max_frequency plus one resolution band
    need <- which(freq <= max_frequency + rbw_t + 64 * .eps)
    padz <- function(M){
      rbind(M, matrix(0, nfft - n_target, ncol(M)))
    }
    Sfull <- mvfft(padz(detrend_linear(sin(phase))))
    S <- Sfull[need, , drop = FALSE]
    scale <- 1 / (Fs * n_target)
    wn <- w1[need]
    sineT <- (Mod(S)^2) * wn * scale
    phaseT <- NULL
    if (target_model == "phase-averaged" || amplitude_mode == "four-group-area"){
      C <- mvfft(padz(detrend_linear(cos(phase))))[need, , drop = FALSE]
      phaseT <- 0.5 * (Mod(S)^2 + Mod(C)^2) * wn * scale
    }
    if (target_model == "coherent-nine"){
      template <- sineT
    } else {
      template <- phaseT
    }
    fr <- freq[need]
    unitPeak <- numeric(no)
    basisRes <- logical(no)
    for (j in 1:no){
      ib <- abs(fr - center[j]) <= rbw_t
      if (center[j] < Fs / 2 && any(ib)){
        unitPeak[j] <- max(template[ib, j])
        basisRes[j] <- is.finite(unitPeak[j]) && unitPeak[j] > 0
      }
    }
    tRay <- rbw_spatial * sr / 100
    tNyq <- (1 / (2 * spacing)) * sr / 100
    eligible <- center >= tRay & center < tNyq & center <= max_frequency & basisRes
    tf <- data_frequency[need] * sr / 100       # identical to fr
    fmask <- fr <= max_frequency + 64 * .eps * max(1, max_frequency)
    if (sum(fmask) < 3){
      stop(sprintf("Fewer than three frequency bins below max_frequency at %.4g cm/kyr", sr))
    }
    rbw_band <- rbw_spatial * sr / 100
    Pn <- P[need, , drop = FALSE]

    if (amplitude_mode == "four-group-area"){
      cand <- matrix(FALSE, length(tf), no)
      for (j in 1:no){
        if (eligible[j]){
          cand[, j] <- abs(tf - center[j]) <= rbw_band &
              tf <= max_frequency + 64 * .eps * max(1, max_frequency)
        }
      }
      active <- eligible & colSums(cand) > 0
      gT <- matrix(0, length(fr), 4)
      for (j in which(active)){
        gT[, group[j]] <- gT[, group[j]] + phaseT[, j]
      }
      bandIdx <- vector("list", 4)
      activeG <- logical(4)
      for (g in 1:4){
        mem <- group == g & active
        if (any(mem)){
          bandIdx[[g]] <- which(rowSums(cand[, mem, drop = FALSE]) > 0)
        } else {
          bandIdx[[g]] <- integer(0)
        }
        activeG[g] <- any(mem) && length(bandIdx[[g]]) > 0
      }
      valid <- TRUE
      upper <- min(max_frequency, freq[nf])
      lo <- pmax(0, center - rbw_band)
      hi <- pmin(upper, center + rbw_band)
      ai <- which(active)
      tol_o <- 128 * .eps * max(c(1, abs(lo), abs(hi), upper))
      if (length(ai) > 1){
        for (a in 1:(length(ai) - 1)){
          for (b in (a + 1):length(ai)){
            i1 <- ai[a]
            i2 <- ai[b]
            if (group[i1] != group[i2] && (min(hi[i1], hi[i2]) - max(lo[i1], lo[i2])) > tol_o){
              valid <- FALSE
            }
          }
        }
      }
      L <- matrix(0, 4, 4)
      nDf <- Fs / nfft
      for (g in 1:4){
        if (length(bandIdx[[g]])){
          L[g, ] <- colSums(gT[bandIdx[[g]], , drop = FALSE]) * nDf
        }
      }
      if (!any(activeG)){
        valid <- FALSE
      } else {
        A <- L[activeG, activeG, drop = FALSE]
        if (any(!is.finite(A)) || any(diag(A) <= 0)){
          valid <- FALSE
        } else {
          rc <- rcond(A)
          if (!is.finite(rc) || rc < 1e-10){
            valid <- FALSE
          }
        }
      }
      nMissing[m] <- sum(!active)
      if (!valid){
        next
      }
      Ms <- L
      for (g in which(!activeG)){
        Ms[g, ] <- 0
        Ms[, g] <- 0
        Ms[g, g] <- 1
      }
      E <- matrix(0, 4, k)
      for (g in 1:4){
        if (length(bandIdx[[g]])){
          E[g, ] <- colSums(Pn[bandIdx[[g]], , drop = FALSE]) * df_spatial
        }
      }
      E[!activeG, ] <- 0
      GP <- .nnls4(Ms, E)
      GP[!activeG, ] <- 0
      mx <- numeric(ncol(GP))
      for (j in seq_len(ncol(GP))){
        mx[j] <- max(GP[, j])
      }
      GP[GP <= matrix(mx * 1e-16, 4, k, byrow = TRUE)] <- 0
      GA <- sqrt(GP)
      Bg <- matrix(0i, length(fr), 4)
      for (g in 1:4){
        mem <- group == g & active
        if (any(mem)){
          Bg[, g] <- rowSums(S[, mem, drop = FALSE])
        }
      }
      Bm <- Bg[fmask, , drop = FALSE]
      tp <- (Mod(Bm %*% GA)^2) * wn[fmask] * scale
    } else {
      cand <- matrix(FALSE, length(tf), no)
      for (j in which(eligible)){
        cand[, j] <- tf >= center[j] - rbw_band & tf <= center[j] + rbw_band & tf <= max_frequency
      }
      multi <- which(rowSums(cand) > 1)
      for (b in multi){
        mem <- which(cand[b, ])
        d <- abs(tf[b] - center[mem])
        tied <- mem[d == min(d)]
        if (length(tied) > 1){
          ch <- tied[which.min(center[tied])]
        } else {
          ch <- tied
        }
        cand[b, mem] <- FALSE
        cand[b, ch] <- TRUE
      }
      bandRes <- colSums(cand) > 0
      if (amplitude_mode == "fixed"){
        active <- eligible
      } else {
        active <- eligible & bandRes
      }
      nMissing[m] <- sum(!active)
      if (amplitude_mode == "fixed"){
        amp <- matrix(fixed_amp * active, no, k)
      } else {
        pk <- matrix(0, no, k)
        for (j in which(active)){
          for (cc in seq_len(k)){
            pk[j, cc] <- max(Pn[cand[, j], cc])
          }
        }
        pk[!is.finite(pk) | pk < 0] <- 0
        amp <- matrix(0, no, k)
        ok <- active & is.finite(unitPeak) & unitPeak > 0
        amp[ok, ] <- sqrt(pk[ok, , drop = FALSE] / unitPeak[ok])
      }
      if (target_model == "coherent-nine"){
        tp <- (Mod(S[fmask, , drop = FALSE] %*% amp)^2) * wn[fmask] * scale
      } else {
        tp <- template[fmask, , drop = FALSE] %*% (amp^2)
      }
    }
    dm <- Pn[fmask, , drop = FALSE]
    rho[m, ] <- .column_correlation(tp, dm, method)
    if (observed_p){
      cp <- cor_with_p(tp[, 1], dm[, 1], method)
      if (is.na(cp["rho"])){
        pval[m] <- 1
      } else {
        pval[m] <- cp["p"]
      }
    }
  }
  list(rho = rho, p = pval, nMissing = nMissing)
}
