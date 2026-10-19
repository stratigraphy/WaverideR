#' @title Plot the result of a bispectral analysis
#'
#' @description The \code{\link{plot_bispectral}} function plots the result of
#' \code{\link{bispectral_analysis}}. The "overview" layout follows Acycle:
#' two aligned 2pi multitaper power spectra (NW = 2, 3 tapers, NFFT = 5N) above
#' maps of the bispectrum magnitude (log10|B|) and the squared bicoherence,
#' with the strongest fraction of each map retained, the significance
#' boundary (if inference was calculated) drawn in black, the strongest
#' (significant) bicoherence peaks annotated and period axes on top.
#' Single maps of "bicoherence", "magnitude", "biphase", "real" or
#' "imaginary" can also be drawn. Rendering is simplified with respect to
#' Acycle (no interpolated display mesh or alpha ramp).
#' Translated from Acycle v3.1 bispectralPlot.m (Mingsong Li).
#'
#' @param result Output of \code{\link{bispectral_analysis}}.
#' @param quantity "overview", "bicoherence", "magnitude", "biphase", "real"
#' or "imaginary". \code{Default="overview"}
#' @param keep_bispectrum Fraction (0,1] of strongest |B| values shown.
#' \code{Default=0.5}
#' @param keep_bicoherence Fraction (0,1] of strongest b^2 values shown.
#' \code{Default=0.5}
#' @param color_grid Number of colours (4-256). \code{Default=32}
#' @param freq_min Lower plotted frequency. \code{Default=0}
#' @param freq_max Upper plotted frequency; \code{NULL} shows the full
#' computed axis. \code{Default=NULL}
#' @param reference_periods Periods drawn as guide lines f1 + f2 = 1/period.
#' \code{Default=NULL}
#' @param frequency_pairs Matrix [f1 f2] of coordinates marked with guides.
#' \code{Default=NULL}
#' @param peak_count Number of bicoherence peaks annotated. \code{Default=5}
#' @param show_period_axes Add a period axis on top of the maps.
#' \code{Default=TRUE}
#' @param show_significance Draw the significance boundary. \code{Default=TRUE}
#' @param title Figure title. \code{Default=NULL}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle bispectral module).
#'
#' @references
#' Li, M., Hinnov, L., Kump, L., 2019. Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences 127,
#' 12-22.
#'
#' @examples
#' \donttest{
#' set.seed(4)
#' t <- 0:1023
#' y <- cos(2 * pi * 0.1 * t) + cos(2 * pi * 0.17 * t) +
#'   0.5 * cos(2 * pi * 0.27 * t) + 0.3 * rnorm(1024)
#' res <- bispectral_analysis(data = cbind(t, y), max_freq_bins = 128)
#' plot_bispectral(res)
#' plot_bispectral(res, quantity = "biphase")
#' }
#' @return
#' Invisibly returns \code{NULL}; the plot is drawn on the current device.
#'
#' @param colour_grid British spelling of \code{color_grid}; if given it overrides \code{color_grid}.
#' @export
#' @importFrom graphics layout
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics contour
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom graphics axis
#' @importFrom graphics box
#' @importFrom graphics abline
#' @importFrom graphics mtext
#' @importFrom graphics segments
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.hold
#' @importFrom grDevices dev.flush

plot_bispectral <- function(result = NULL,
                            quantity = "overview",
                            keep_bispectrum = 0.5,
                            keep_bicoherence = 0.5,
                            color_grid = 32,
                            freq_min = 0,
                            freq_max = NULL,
                            reference_periods = NULL,
                            frequency_pairs = NULL,
                            peak_count = 5,
                            show_period_axes = TRUE,
                            show_significance = TRUE,
                            title = NULL,
  colour_grid = NULL){
  if (!is.null(colour_grid)) color_grid <- colour_grid  # British spelling synonym
  # ---- internal helper functions (local to this function) ----

  # ---------------------------------------------------------------------------
  # small utilities
  # ---------------------------------------------------------------------------

  acy_nextpow2 <- function(n){
    if (n <= 1){
      return(0)
    }
    ceiling(log2(n) - 1e-12)
  }

  acy_dpss <- function(n, nw, k = NULL){
    if (is.null(k)){
      k <- round(2 * nw) - 1
    }
    key <- paste(n, format(nw, digits = 17), k, sep = "_")
    if (exists(key, envir = acy_dpss_cache, inherits = FALSE)){
      return(get(key, envir = acy_dpss_cache))
    }
    w <- nw / n
    i <- 0:(n - 1)
    dg <- ((n - 1 - 2 * i) / 2)^2 * cos(2 * pi * w)
    od <- (1:(n - 1)) * (n - (1:(n - 1))) / 2
    if (n <= 2048){
      Tm <- diag(dg, n, n)
      for (j in 1:(n - 1)){
        Tm[j, j + 1] <- od[j]
        Tm[j + 1, j] <- od[j]
      }
      ev <- eigen(Tm, symmetric = TRUE)
      E <- ev$vectors[, 1:k, drop = FALSE]
    } else {
      E <- acy_dpss_tridiag(dg, od, k)
    }
    tc <- i - (n - 1) / 2
    for (j in 1:k){
      E[, j] <- E[, j] / sqrt(sum(E[, j]^2))
      if (j %% 2 == 1){
        if (sum(E[, j]) < 0){
          E[, j] <- -E[, j]
        }
      } else {
        if (sum(tc * E[, j]) < 0){
          E[, j] <- -E[, j]
        }
      }
    }
    V <- numeric(k)
    for (j in 1:k){
      V[j] <- acy_dpss_concentration(E[, j], w)
    }
    out <- list(E = E, V = V)
    assign(key, out, envir = acy_dpss_cache)
    out
  }

  acy_dpss_concentration <- function(e, w){
    n <- length(e)
    m <- 2^acy_nextpow2(2 * n)
    fe <- stats::fft(c(e, rep(0, m - n)))
    r <- Re(stats::fft(Mod(fe)^2, inverse = TRUE)) / m
    lag <- 1:(n - 1)
    a <- sin(2 * pi * w * lag) / (pi * lag)
    2 * w * r[1] + 2 * sum(a * r[lag + 1])
  }

  acy_dpss_tridiag <- function(dg, od, k){
    n <- length(dg)
    sturm_count <- function(lambda){
      cnt <- 0
      q <- dg[1] - lambda
      if (q < 0){
        cnt <- cnt + 1
      }
      for (j in 2:n){
        if (q == 0){
          q <- 1e-300
        }
        q <- dg[j] - lambda - od[j - 1]^2 / q
        if (q < 0){
          cnt <- cnt + 1
        }
      }
      cnt
    }
    bound <- max(abs(dg)) + 2 * max(abs(od))
    E <- matrix(0, n, k)
    for (j in 1:k){
      target <- n - j
      lo <- -bound
      hi <- bound
      for (it in 1:200){
        mid <- (lo + hi) / 2
        if (sturm_count(mid) > target){
          hi <- mid
        } else {
          lo <- mid
        }
        if (hi - lo <= 1e-14 * max(1, abs(mid))){
          break
        }
      }
      lambda <- (lo + hi) / 2 + 1e-10 * max(1, abs(hi))
      v <- rep(1, n) / sqrt(n)
      for (it in 1:6){
        v <- acy_tridiag_solve(od, dg - lambda, od, v)
        v <- v / sqrt(sum(v^2))
      }
      E[, j] <- v
    }
    E
  }

  acy_tridiag_solve <- function(a, b, c, d){
    n <- length(b)
    cp <- numeric(n)
    dp <- numeric(n)
    cp[1] <- c[1] / b[1]
    dp[1] <- d[1] / b[1]
    for (i in 2:n){
      den <- b[i] - a[i - 1] * cp[i - 1]
      if (den == 0){
        den <- 1e-300
      }
      if (i < n){
        cp[i] <- c[i] / den
      }
      dp[i] <- (d[i] - a[i - 1] * dp[i - 1]) / den
    }
    x <- numeric(n)
    x[n] <- dp[n]
    for (i in (n - 1):1){
      x[i] <- dp[i] - cp[i] * x[i + 1]
    }
    x
  }

  # ---------------------------------------------------------------------------
  # MATLAB periodogram / pmtm (one-sided, adaptive weights)
  # ---------------------------------------------------------------------------

  acy_datawrap <- function(x, nfft){
    n <- length(x)
    if (n <= nfft){
      return(c(x, rep(0, nfft - n)))
    }
    out <- rep(0, nfft)
    nseg <- ceiling(n / nfft)
    xx <- c(x, rep(0, nseg * nfft - n))
    for (s in 1:nseg){
      out <- out + xx[((s - 1) * nfft + 1):(s * nfft)]
    }
    out
  }

  acy_onesided <- function(S, nfft){
    if (nfft %% 2 == 0){
      sel <- 1:(nfft / 2 + 1)
      P <- S[sel]
      if (length(P) > 2){
        P[2:(length(P) - 1)] <- 2 * P[2:(length(P) - 1)]
      }
    } else {
      sel <- 1:((nfft + 1) / 2)
      P <- S[sel]
      if (length(P) > 1){
        P[2:length(P)] <- 2 * P[2:length(P)]
      }
    }
    P
  }

  acy_freq_vector <- function(nfft, fs = NULL){
    if (nfft %% 2 == 0){
      m <- nfft / 2 + 1
    } else {
      m <- (nfft + 1) / 2
    }
    if (is.null(fs)){
      return(2 * pi * (0:(m - 1)) / nfft)
    }
    fs * (0:(m - 1)) / nfft
  }

  acy_pmtm <- function(x, nw = 4, nfft = NULL, fs = NULL, drop_last = TRUE,
                       method = "adapt"){
    x <- as.numeric(x)
    n <- length(x)
    if (is.null(nfft)){
      nfft <- max(256, 2^acy_nextpow2(n))
    }
    if (drop_last == TRUE){
      k <- round(2 * nw) - 1
    } else {
      k <- round(2 * nw)
    }
    if (k < 1){
      k <- 1
    }
    tp <- acy_dpss(n, nw, k)
    E <- tp$E
    V <- tp$V
    Sk <- matrix(0, nfft, k)
    for (j in 1:k){
      Sk[, j] <- Mod(stats::fft(acy_datawrap(E[, j] * x, nfft)))^2
    }
    if (method == "adapt" && k > 1){
      sig2 <- sum(x^2) / n
      S <- (Sk[, 1] + Sk[, 2]) / 2
      S1 <- rep(0, nfft)
      tol <- 0.0005 * sig2 / nfft
      a <- sig2 * (1 - V)
      it <- 0
      while (sum(abs(S - S1) / nfft) > tol && it < 100){
        it <- it + 1
        num <- 0
        den <- 0
        for (j in 1:k){
          b <- S / (S * V[j] + a[j])
          wk <- b^2 * V[j]
          num <- num + wk * Sk[, j]
          den <- den + wk
        }
        S1 <- S
        S <- num / den
      }
    } else if (method == "eigen"){
      S <- rep(0, nfft)
      for (j in 1:k){
        S <- S + V[j] * Sk[, j]
      }
      S <- S / sum(V)
    } else {
      S <- rowMeans(Sk)
    }
    P <- acy_onesided(S, nfft)
    if (is.null(fs)){
      P <- P / (2 * pi)
    } else {
      P <- P / fs
    }
    list(pxx = P, f = acy_freq_vector(nfft, fs), tapers = E, eigenvalues = V)
  }

  # @keywords internal
  # @noRd
  bisp_round <- function(x){
    # MATLAB round (half away from zero)
    sign(x) * floor(abs(x) + 0.5)
  }

  # ---------------------------------------------------------------------------
  # plotting helpers (bispectralPlot.m, simplified)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_mtm_power <- function(result){
    y <- result$processed_data[, 2]
    y <- y - mean(y)
    n <- length(y)
    mt <- acy_pmtm(y, nw = 2, nfft = 5 * n)
    dt <- result$meta$sample_interval
    list(f = mt$f / (2 * pi * dt), p = mt$pxx)
  }

  # @keywords internal
  # @noRd
  bisp_log_power_limits <- function(p, headroom){
    if (length(p) == 0){
      return(c(1e-12, 1))
    }
    lo <- min(p)
    hi <- max(p)
    if (hi <= lo * (1 + 64 * .Machine$double.eps)){
      return(c(lo, lo * 10))
    }
    up <- exp(log(lo) + (log(hi) - log(lo)) / (1 - headroom))
    c(lo, up)
  }

  # @keywords internal
  # @noRd
  bisp_positive_map <- function(n){
    grDevices::colorRampPalette(c("#FFFFFF", "#FAD1B8", "#EB6133", "#A61414"))(n)
  }

  # @keywords internal
  # @noRd
  bisp_diverging_map <- function(n){
    grDevices::colorRampPalette(c("#1A3D9E", "#6BADD6", "#FFFFFF", "#F29461",
                                  "#A61414"))(n)
  }

  # @keywords internal
  # @noRd
  bisp_robust_limits <- function(z){
    v <- sort(z[is.finite(z)])
    if (length(v) == 0){
      return(c(0, 1))
    }
    lo <- v[max(1, bisp_round(0.01 * length(v)))]
    hi <- v[max(1, bisp_round(0.995 * length(v)))]
    if (!(hi > lo)){
      d <- max(1, abs(lo)) * 1e-6
      return(c(lo - d, lo + d))
    }
    c(lo, hi)
  }

  # @keywords internal
  # @noRd
  bisp_symmetric_limits <- function(z){
    v <- sort(abs(z[is.finite(z)]))
    if (length(v) == 0){
      return(c(-1, 1))
    }
    m <- v[max(1, bisp_round(0.995 * length(v)))]
    if (!(is.finite(m) && m > 0)){
      m <- 1
    }
    c(-m, m)
  }

  # @keywords internal
  # @noRd
  bisp_draw_map <- function(result, quantity, keep, ngrid, fmin, fmax,
                            reference_periods, frequency_pairs, peak_count,
                            show_period_axes, show_significance){
    f <- result$frequency
    rescale <- FALSE
    if (quantity %in% c("bicoherence", "bicoherence-squared", "b2")){
      quantity <- "bicoherence"
      z <- result$bicoherence_squared
      metric <- z
      rescale <- TRUE
      lab <- "b^2"
      ttl <- "Magnitude-squared bicoherence, b^2"
      lim <- c(0, 1)
    } else if (quantity %in% c("magnitude", "bispectrum-magnitude")){
      quantity <- "magnitude"
      raw <- result$bispectrum_magnitude
      fp <- raw[is.finite(raw) & raw > 0]
      if (length(fp) == 0){
        floor_value <- .Machine$double.xmin
      } else {
        floor_value <- max(.Machine$double.xmin, min(fp) * 0.1)
      }
      z <- matrix(NA_real_, nrow(raw), ncol(raw))
      ok <- is.finite(raw)
      z[ok] <- log10(pmax(raw[ok], floor_value))
      metric <- z
      rescale <- TRUE
      lab <- "log10|B|"
      ttl <- "Bispectrum magnitude, log10|B|"
      lim <- bisp_robust_limits(z)
    } else if (quantity %in% c("real", "bispectrum-real")){
      z <- Re(result$bispectrum)
      metric <- result$bispectrum_magnitude
      lab <- "Re(B)"
      ttl <- "Real bispectrum"
      lim <- bisp_symmetric_limits(z)
    } else if (quantity %in% c("imaginary", "bispectrum-imaginary")){
      z <- Im(result$bispectrum)
      metric <- result$bispectrum_magnitude
      lab <- "Im(B)"
      ttl <- "Imaginary bispectrum"
      lim <- bisp_symmetric_limits(z)
    } else if (quantity %in% c("biphase", "phase")){
      quantity <- "biphase"
      z <- result$biphase
      metric <- result$bicoherence_squared
      lab <- "Biphase (rad)"
      ttl <- "Biphase"
      lim <- c(-pi, pi)
    } else {
      stop("Unknown quantity for plot_bispectral().")
    }
    cols <- bisp_diverging_map(ngrid)
    f1g <- matrix(f, length(f), length(f), byrow = TRUE)
    f2g <- matrix(f, length(f), length(f))
    vis <- f1g >= fmin & f1g <= fmax & f2g >= fmin & f2g <= fmax
    cutoff <- -Inf
    if (keep < 1){
      vm <- is.finite(z) & is.finite(metric) & vis
      fv <- sort(metric[vm])
      if (length(fv) > 0){
        first <- max(1, min(length(fv), floor((1 - keep) * length(fv)) + 1))
        cutoff <- fv[first]
        if (rescale == TRUE){
          rc <- z[vm & metric >= cutoff]
          lo <- min(rc)
          hi <- max(rc)
          if (!(hi > lo)){
            lo <- hi - max(1, abs(hi)) * 1e-6
          }
          lim <- c(lo, hi)
          cols <- bisp_positive_map(ngrid)
        }
      }
    }
    zd <- z
    hide <- !(is.finite(metric) & metric >= cutoff)
    zd[hide] <- NA
    zd <- pmin(pmax(zd, lim[1]), lim[2])
    graphics::par(mar = c(4.2, 4.5, 4.2, 1.5))
    graphics::image(f, f, t(zd), col = cols, zlim = lim, xlim = c(fmin, fmax),
                    ylim = c(fmin, fmax), xlab = "f1 (cycles/unit)",
                    ylab = "f2 (cycles/unit)", useRaster = FALSE, asp = NA)
    if (quantity %in% c("bicoherence", "magnitude") && any(is.finite(zd))){
      lev <- seq(lim[1], lim[2], length.out = 7)[2:6]
      graphics::contour(f, f, t(zd), levels = lev, add = TRUE, drawlabels = FALSE,
                        col = "#4D3D33", lwd = 0.45)
    }
    has_inf <- !is.null(result$significance) && result$significance$method != "none"
    if (show_significance == TRUE && has_inf && any(result$significant_mask)){
      src <- matrix(NA_real_, length(f), length(f))
      src[result$valid_mask] <- as.numeric(result$significant_mask[result$valid_mask])
      graphics::contour(f, f, t(src), levels = 0.5, add = TRUE, drawlabels = FALSE,
                        col = "black", lwd = 2)
    }
    if (peak_count > 0 && quantity == "bicoherence"){
      b2 <- result$bicoherence_squared
      cand <- result$valid_mask & is.finite(b2) & !hide & vis
      if (has_inf){
        cand <- cand & result$significant_mask
      }
      idx <- which(cand)
      if (length(idx) > 0){
        idx <- idx[order(-b2[idx], method = "radix")]
        chosen <- matrix(0, 0, 2)
        for (i in seq_len(length(idx))){
          r <- (idx[i] - 1) %% length(f) + 1
          cidx <- (idx[i] - 1) %/% length(f) + 1
          p1 <- f[cidx]
          p2 <- f[r]
          np <- c((p1 - fmin) / (fmax - fmin), (p2 - fmin) / (fmax - fmin))
          far <- TRUE
          if (nrow(chosen) > 0){
            dd <- sqrt((chosen[, 1] - np[1])^2 + (chosen[, 2] - np[2])^2)
            far <- all(dd >= 0.18)
          }
          if (far){
            chosen <- rbind(chosen, np)
            graphics::points(p1, p2, pch = 21, bg = "white", cex = 0.8)
            graphics::text(p1, p2, sprintf("  %.3g, %.3g", p1, p2), pos = 4,
                           cex = 0.75, offset = 0.2)
            if (nrow(chosen) >= peak_count){
              break
            }
          }
        }
      }
    }
    if (!is.null(reference_periods)){
      for (i in seq_len(length(reference_periods))){
        s <- 1 / reference_periods[i]
        graphics::abline(a = s, b = -1, lty = 2, col = "grey30")
      }
    }
    if (!is.null(frequency_pairs)){
      frequency_pairs <- matrix(frequency_pairs, ncol = 2)
      for (i in seq_len(nrow(frequency_pairs))){
        graphics::abline(v = frequency_pairs[i, 1], h = frequency_pairs[i, 2],
                         lty = 3, col = "grey20")
      }
    }
    if (show_period_axes == TRUE){
      at <- pretty(c(fmin, fmax))
      at <- at[at > 0 & at >= fmin & at <= fmax]
      if (length(at) > 0){
        graphics::axis(3, at = at, labels = signif(1 / at, 3), cex.axis = 0.8)
        graphics::mtext("Period (unit)", side = 3, line = 1.9, cex = 0.7)
      }
    }
    graphics::box()
    list(title = ttl, label = lab, limits = lim, colors = cols)
  }

  # @keywords internal
  # @noRd
  bisp_colorbar <- function(map){
    graphics::par(mar = c(2.2, 6, 0.5, 2))
    lv <- seq(map$limits[1], map$limits[2], length.out = length(map$colors) + 1)
    mid <- (lv[-1] + lv[-length(lv)]) / 2
    graphics::image(mid, 1, matrix(mid, ncol = 1), col = map$colors,
                    zlim = map$limits, xlab = "", ylab = "", yaxt = "n")
    graphics::mtext(map$label, side = 2, las = 1, line = 0.5, cex = 0.8)
    graphics::box()
  }

  # ---------------------------------------------------------------------------
  # Slepian (DPSS) tapers, MATLAB dpss() conventions (unit energy)
  # ---------------------------------------------------------------------------

  acy_dpss_cache <- new.env(parent = emptyenv())

  # ---- main code ----
  if (is.null(result) || is.null(result$bicoherence_squared)){
    stop("result must be the output of bispectral_analysis()")
  }
  if (!(keep_bispectrum > 0 && keep_bispectrum <= 1 &&
        keep_bicoherence > 0 && keep_bicoherence <= 1)){
    stop("keep fractions must lie in (0,1]")
  }
  if (!(color_grid >= 4 && color_grid <= 256 && color_grid == round(color_grid))){
    stop("color_grid must be an integer from 4 through 256")
  }
  quantity <- tolower(trimws(quantity))
  f <- result$frequency
  fmin <- max(freq_min, 0)
  if (is.null(freq_max)){
    fmax <- f[length(f)]
  } else {
    fmax <- min(freq_max, f[length(f)])
  }
  if (!(fmax > fmin)){
    fmin <- f[1]
    fmax <- f[length(f)]
  }
  if (is.null(title)){
    title <- "Bispectral analysis"
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  grDevices::dev.hold()
  on.exit(grDevices::dev.flush(), add = TRUE)

  if (quantity == "overview"){
    graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = TRUE),
                     heights = c(1.1, 3, 0.55))
    graphics::par(oma = c(0, 0, 2, 0))
    mtm <- bisp_mtm_power(result)
    vis <- mtm$f >= fmin & mtm$f <= fmax & is.finite(mtm$p) & mtm$p > 0
    ylim <- bisp_log_power_limits(mtm$p[vis], 0.2)
    ttl <- c("Complex bispectrum amplitude, |B|",
             "Magnitude-squared bicoherence, b^2")
    for (k in 1:2){
      graphics::par(mar = c(2.2, 4.5, 2, 1.5))
      pp <- mtm$p
      pp[!(is.finite(pp) & pp > 0)] <- NA
      graphics::plot(mtm$f, pp, type = "l", log = "y", xlim = c(fmin, fmax),
                     ylim = ylim, xaxs = "i", xlab = "",
                     ylab = "2pi MTM power", main = ttl[k], font.main = 1,
                     cex.main = 1)
    }
    m1 <- bisp_draw_map(result, "magnitude", keep_bispectrum, color_grid,
                        fmin, fmax, reference_periods, frequency_pairs,
                        peak_count, show_period_axes, show_significance)
    m2 <- bisp_draw_map(result, "bicoherence", keep_bicoherence, color_grid,
                        fmin, fmax, reference_periods, frequency_pairs,
                        peak_count, show_period_axes, show_significance)
    bisp_colorbar(m1)
    bisp_colorbar(m2)
    graphics::mtext(title, side = 3, outer = TRUE, font = 2, line = 0.5)
  } else {
    graphics::layout(matrix(c(1, 2), 2, 1), heights = c(4, 0.6))
    if (quantity %in% c("magnitude", "real", "imaginary")){
      keep <- keep_bispectrum
    } else {
      keep <- keep_bicoherence
    }
    m <- bisp_draw_map(result, quantity, keep, color_grid, fmin, fmax,
                       reference_periods, frequency_pairs, peak_count,
                       show_period_axes, show_significance)
    graphics::mtext(paste0(title, ": ", m$title), side = 3, line = 2.9,
                    font = 2)
    bisp_colorbar(m)
  }
  invisible(NULL)
}
