## Origin: astro_parameters.R (secular frequencies, precession constant k and Earth-Moon state
## from measured obliquity/precession periods). Released under the MIT licence.

# Plotting

#' Six-panel summary of an inversion
#' @param fit Inversion result to plot: the list returned by \code{\link{astro_invert}} (with elements \code{k},
#' \code{state}, \code{draws} and \code{windows}). Required, no default.
#' @param peaks Measured spectral peak periods in kyr (numeric vector). If given, panel (b) shows the
#' \code{\link{k_scan}} of these peaks (number of terms matched by each precession constant k); \code{NULL} leaves
#' panel (b) empty. \code{Default=NULL}.
#' @param age_ma Age of the record in Ma (numeric, single value). If given, the AstroGeo22 (Farhat et al. 2022) value
#' of k at that age is added to panel (c) and the fitted Earth-Moon distance and length of day are plotted at that age
#' in panels (d) and (e); \code{NULL} skips these comparisons. \code{Default=NULL}.
#' @param file Path of an output file (character). A name ending in \code{".png"} writes a PNG file, any other name a
#' PDF file; \code{NULL} draws on the current device. \code{Default=NULL}.
#' @param width Width of the output file in inches (numeric, > 0); only used when \code{file} is given.
#' \code{Default=11}.
#' @param height Height of the output file in inches (numeric, > 0); only used when \code{file} is given.
#' \code{Default=7}.
#' @param res Resolution of the PNG output in dpi (numeric, > 0); only used when \code{file} ends in \code{".png"}.
#' \code{Default=300}.
#' @export
plot_astro <- function(fit, peaks = NULL, age_ma = NULL, file = NULL,
    width = 11, height = 7, res = 300){
  # Internal helper functions (local to this function)

  period_kyr  <- function(freq_arcsec){
    ARCSEC_PER_CYCLE / abs(freq_arcsec) / 1000
  }

  freq_arcsec <- function(period_kyr){
    ARCSEC_PER_CYCLE / (period_kyr * 1000)
  }

  ## the secular offset of each k-bearing term, so that f_term = p + offset
  K_TERMS <- function(gs = LA2004){
    c("k+s3" = unname(gs$s[3]), "k+g5" = unname(gs$g[5]),
        "k+g2" = unname(gs$g[2]), "k+g4" = unname(gs$g[4]),
        "k+g3" = unname(gs$g[3]))
  }

  astrogeo22_nearest <- function(age_ma = NULL, file = NULL, url = ASTROGEO22_URL){
    if (is.null(file)){
      file <- file.path(tempdir(), "AstroGeo22.txt")
      if (!file.exists(file)){
        try(utils::download.file(url, file, quiet = TRUE), silent = TRUE)
      }
    }
    if (!file.exists(file)){
      message("AstroGeo22.txt not available (no local copy and no network); ",
          "comparisons with the tidal model are skipped.")
      return(NULL)
    }
    A <- utils::read.table(file, comment.char = "%")
    out <- data.frame(age_ma = A[[1]] * 1000, a_RE = A[[2]], lod_solar_h = A[[5]],
        obliquity_deg = A[[8]], p_arcsec = A[[11]])
    out$k_kyr <- period_kyr(out$p_arcsec)
    out$days_per_year <- days_per_year(out$lod_solar_h)
    if (is.null(age_ma)){
      return(out)
    }
    idx <- numeric(length(age_ma))
    for (i in seq_along(age_ma)){
      idx[i] <- which.min(abs(out$age_ma - age_ma[[i]]))
    }
    names(idx) <- names(age_ma)
    out[idx, ]
  }

  ARCSEC_PER_CYCLE <- 360 * 3600                 # 1 296 000

  # Secular frequencies
  ## La2004, Laskar et al. (2004) table 3, in arcsec/yr.  These are the values
  ## Acycle bundles as La04gisi.csv and uses in calculate_orbit9.m.

  LA2004 <- list(
      g = c(g1 = 5.59, g2 = 7.452, g3 = 17.368, g4 = 17.916, g5 = 4.257452,
          g6 = 28.245, g7 = 3.087951, g8 = 0.673021, g9 = -0.34994),
      s = c(s1 = -5.59, s2 = -7.05, s3 = -18.85, s4 = -17.755, s5 = 0.00000013,
          s6 = -26.347855, s7 = -2.9925259, s8 = -0.691736, s9 = -0.34998))

  # Reference models

  ## AstroGeo22 (Farhat et al. 2022).  Supply a local copy of AstroGeo22.txt, or
  ## let the function fetch it from the Acycle repository.  Columns used:
  ## 1 time (Ga), 2 a (R_E), 5 length of day (h), 8 obliquity (deg),
  ## 11 precession frequency (arcsec/yr).
  ASTROGEO22_URL <- paste0("https://raw.githubusercontent.com/mingsongli/acycle/",
      "master/code/package/secular_periods/AstroGeo22.txt")

  ## Zeeden et al. (2023), their table 1, combined Devonian estimate (1 sigma)
  ZEEDEN2023 <- list(k_kyr = c(mean = 20.70, sd = 0.22),
      periods = c("k+g5" = 19.35, "k+g2" = 18.47,
          "k+g4" = 16.13, "k+g3" = 16.18, "k+s3" = 29.50),
      sd = c("k+g5" = 0.20, "k+g2" = 0.18, "k+g4" = 0.23,
          "k+g3" = 0.16, "k+s3" = 0.46))

  # Main code
  if (!is.null(file)){
    if (grepl("\\.png$", file)){
      png(file, width * res, height * res, res = res)
    } else {
      pdf(file, width = width, height = height)
    }
    on.exit(dev.off())
  }
  AG <- astrogeo22_nearest()
  COL <- c(obl = "#7570B3", p5 = "#1B9E77", p2 = "#D95F02", p4 = "#E7298A")
  terms <- names(fit$windows)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3.0, 1.0), mgp = c(2.5, 0.6, 0),
      cex.axis = 0.85, cex.lab = 0.95)

  ## (a) term periods as a function of k, with the measured windows
  kk <- seq(15, 30, by = 0.05)
  off <- K_TERMS()[terms]
  plot(NA, xlim = range(kk), ylim = c(10, 45), xlab = "precession constant k (kyr)",
      ylab = "period (kyr)", main = "(a) terms as a function of k")
  for (i in seq_along(terms)){
    lines(kk, period_kyr(freq_arcsec(kk) + off[[i]]), col = COL[i], lwd = 2)
    w <- fit$windows[[terms[i]]]
    rect(fit$k[["2.5%"]], min(w), fit$k[["97.5%"]], max(w),
        col = paste0(COL[i], "40"), border = NA)
  }
  abline(v = fit$k[["mean"]], lty = 2)
  legend("topleft", terms, col = COL[seq_along(terms)], lwd = 2, bty = "n", cex = 0.75)

  ## (b) how many terms a single k can satisfy
  if (!is.null(peaks)){
    sc <- k_scan(peaks, terms = intersect(terms, c("k+s3","k+g5","k+g2","k+g4","k+g3")))
    plot(sc$k_kyr, sc$matched, type = "l", lwd = 1.6, xlab = "k (kyr)",
        ylab = "terms matched", main = "(b) single-k match score",
        ylim = c(0, length(terms)))
    abline(v = fit$k[["mean"]], col = "#B2182B", lty = 2)
  } else {
    plot.new()
  }

  ## (c) the posterior of k
  hist(fit$draws$k_kyr, breaks = 40, col = "grey85", border = "white",
      xlab = "k (kyr)", main = "(c) precession constant")
  abline(v = fit$k[c("2.5%", "mean", "97.5%")], lty = c(3, 1, 3), lwd = c(1, 2, 1))
  if (!is.null(AG) && !is.null(age_ma)){
    abline(v = astrogeo22_nearest(age_ma)$k_kyr, col = "grey40", lty = 2, lwd = 2)
  }
  abline(v = ZEEDEN2023$k_kyr[["mean"]], col = "#B2182B", lty = 3, lwd = 2)
  legend("topleft", c("this study", "AstroGeo22", "Zeeden et al. 2023"),
      col = c("black", "grey40", "#B2182B"), lty = c(1, 2, 3), lwd = 2,
      bty = "n", cex = 0.7)

  ## (d) Earth-Moon distance through time
  if (!is.null(AG)){
    k <- AG$age_ma <= 1000
    plot(AG$age_ma[k], AG$a_RE[k], type = "l", lwd = 2, col = "grey40",
        xlab = "age (Ma)", ylab = expression("Earth-Moon distance ("*R[E]*")"),
        main = "(d) Earth-Moon distance")
    if (!is.null(age_ma)){
      s <- fit$state["a_RE", ]
      arrows(age_ma, s[["2.5%"]], age_ma, s[["97.5%"]], length = 0.03,
          angle = 90, code = 3, col = "#0072B2", lwd = 2)
      points(age_ma, s[["mean"]], pch = 16, col = "#0072B2", cex = 1.2)
    }
  } else {
    plot.new()
  }

  ## (e) length of day through time
  if (!is.null(AG)){
    k <- AG$age_ma <= 1000
    plot(AG$age_ma[k], AG$lod_solar_h[k], type = "l", lwd = 2, col = "grey40",
        xlab = "age (Ma)", ylab = "length of solar day (h)",
        main = "(e) length of day")
    if (!is.null(age_ma)){
      s <- fit$state["lod_solar_h", ]
      arrows(age_ma, s[["2.5%"]], age_ma, s[["97.5%"]], length = 0.03,
          angle = 90, code = 3, col = "#0072B2", lwd = 2)
      points(age_ma, s[["mean"]], pch = 16, col = "#0072B2", cex = 1.2)
    }
  } else {
    plot.new()
  }

  ## (f) measured against predicted periods at the fitted k
  pred <- period_kyr(freq_arcsec(fit$k[["mean"]]) + off)
  obs  <- numeric(length(fit$windows))
  err  <- numeric(length(fit$windows))
  for (i in seq_along(fit$windows)){
    w <- fit$windows[[i]]
    obs[i] <- mean(w)
    err[i] <- diff(range(w)) / 2
  }
  names(obs) <- names(fit$windows)
  names(err) <- names(fit$windows)
  rng  <- range(c(pred, obs)) + c(-2, 2)
  plot(pred, obs, pch = 16, col = COL[seq_along(terms)], cex = 1.3, xlim = rng,
      ylim = rng, xlab = "predicted period (kyr)", ylab = "measured period (kyr)",
      main = "(f) measured against predicted")
  arrows(pred, obs - err, pred, obs + err, length = 0.03, angle = 90, code = 3,
      col = COL[seq_along(terms)])
  abline(0, 1, lty = 2, col = "grey50")
  text(pred, obs, terms, pos = 4, cex = 0.7, col = "grey25")
  invisible(fit)
}
