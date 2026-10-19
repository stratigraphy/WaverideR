#' @title CTest: astronomical parameters from interpreted spectral peaks
#'
#' @param data 8x2 table (matrix, data.frame or tibble) as read from the
#'   published spreadsheet, or built with ctest_input()
#' @param nsim number of successful samples to collect
#' @param mode "PT", "P" or "T" (see file header)
#' @param ptest range of precession constants to sample, in years
#' @param sampleorb "U" uniform within 2 sigma (used in the paper) or "N" normal
#' @param seed optional RNG seed
#' @param max_attempts abort if no nsim successes are found within this many draws
#' @param verbose print progress every 10\%
#' @param engine "vectorized" (default, draws in blocks; fast) or "loop",
#'   which follows the published draw order exactly and reproduces
#'   mmc3/mmc4/mmc5 sample for sample under the same seed
#' @param block number of draws per block for the vectorized engine
#' @param s3_uncertainty "published" keeps the published code's draw for s3,
#'   rnorm(1, s3, ds3*sqrt(t/1e8)), in which the relative uncertainty ds3 is
#'   applied as an absolute one in years, so s3 is effectively fixed.
#'   "scaled" applies it the same way as for g2-g5 (s3 * ds3 * sqrt(t/1e8)),
#'   which gives s3 a realistic spread. Use "published" to reproduce the paper.
#' @return list with `summary` (table as in the paper), `samples`
#'   (all accepted samples), `sampled_range` of g/s, and `acceptance`
#' @export
ctest <- function(data, nsim = 1000, mode = c("PT", "P", "T"),
                  ptest = c(10000, 28000), sampleorb = c("U", "N"),
                  seed = NULL, max_attempts = 5e7, verbose = TRUE,
                  engine = c("vectorized", "loop"), block = 1e5,
                  s3_uncertainty = c("published", "scaled")){
  mode <- match.arg(mode)
  sampleorb <- match.arg(sampleorb)
  engine <- match.arg(engine)
  s3_uncertainty <- match.arg(s3_uncertainty)
  if (!is.null(seed)){
    set.seed(seed)
  }
  d <- function(i, j){
    as.numeric(as.matrix(data)[i, j])
  }

  tmin <- d(1, 1) * 1e6
  tmax <- d(1, 2) * 1e6
  high_405 <- 1 / d(2, 1)
  low_405 <- 1 / d(2, 2)
  high_125 <- 1 / d(3, 1)
  low_125 <- 1 / d(3, 2)
  high_95 <- 1 / d(4, 1)
  low_95 <- 1 / d(4, 2)
  o1_min <- 1 / d(5, 1)
  o1_max <- 1 / d(5, 2)
  p1_min <- 1 / d(6, 1)
  p1_max <- 1 / d(6, 2)
  p2_min <- 1 / d(7, 1)
  p2_max <- 1 / d(7, 2)
  p3_min <- 1 / d(8, 1)
  p3_max <- 1 / d(8, 2)

  useT <- mode %in% c("PT", "T")
  useP <- mode %in% c("PT", "P")
  if (useT && any(!is.finite(c(o1_min, o1_max)))){
    stop("mode '", mode, "' needs an obliquity range")
  }
  if (useP && any(!is.finite(c(p1_min, p1_max, p2_min, p2_max, p3_min, p3_max)))){
    stop("mode '", mode, "' needs the three precession ranges")
  }

  # secular frequencies as periods in yr (Laskar et al. 2011)
  g2 <- 173913
  g3 <- 74620
  g4 <- 72338
  g5 <- 304407
  s3 <- 68753
  # d100 uncertainties of Laskar et al. (2011, Table 6), arcsec/yr, relative
  dg2 <- 0.019 / 7.453
  dg3 <- 0.2 / 17.368
  dg4 <- 0.2 / 17.916
  dg5 <- 0.000030 / 4.257482
  ds3 <- 0.066 / 18.848

  # sedimentation-rate bounds from the eccentricity interpretation (m/kyr)
  SRul <- min(high_405 / 405, high_125 / 110, high_95 / 90)
  SRll <- max(low_405 / 406.18, low_125 / 140, low_95 / 105)

  keep <- matrix(NA_real_, nsim, 19)
  colnames(keep) <- c("SR", "Ecc405", "Ecc124", "Ecc130", "Ecc95", "Ecc99",
                      "Obl", "--", "Prec1", "Prec2", "Prec3", "Prec4", "--",
                      "g2", "g3", "g4", "g5", "s3", "p")
  SE <- matrix(c(g2, g3, g4, g5, s3), 5, 2,
               dimnames = list(c("g2", "g3", "g4", "g5", "s3"), c("min", "max")))
  i <- 1
  attempts <- 0
  report <- max(1, round(nsim / 10))

  accepted_total <- 0
  if (engine == "vectorized"){
    block <- max(1000, min(block, 5e6))
    while (i <= nsim){
      attempts <- attempts + block
      if (attempts > max_attempts){
        stop("No consistent samples found; check the interpreted cycle ranges.")
      }
      SR <- stats::runif(block, floor(min(SRll, SRul) * 1000), ceiling(max(SRul, SRll) * 1000)) / 1000
      t <- stats::runif(block, tmin, tmax)
      p <- stats::runif(block, ptest[1], ptest[2])
      sc <- sqrt(t / 1e8)
      if (sampleorb == "N"){
        G2 <- stats::rnorm(block, g2, g2 * dg2 * sc)
        G3 <- stats::rnorm(block, g3, g3 * dg3 * sc)
        G4 <- stats::rnorm(block, g4, g4 * dg4 * sc)
        G5 <- stats::rnorm(block, g5, g5 * dg5 * sc)
      } else {
        G2 <- stats::runif(block, g2 - g2 * 2 * dg2 * sc, g2 + g2 * 2 * dg2 * sc)
        G3 <- stats::runif(block, g3 - g3 * 2 * dg3 * sc, g3 + g3 * 2 * dg3 * sc)
        G4 <- stats::runif(block, g4 - g4 * 2 * dg4 * sc, g4 + g4 * 2 * dg4 * sc)
        G5 <- stats::runif(block, g5 - g5 * 2 * dg5 * sc, g5 + g5 * 2 * dg5 * sc)
      }
      if (s3_uncertainty == "scaled"){
        sd_s3 <- s3 * ds3 * sc
      } else {
        sd_s3 <- ds3 * sc
      }
      if (useT){
        S3 <- stats::rnorm(block, s3, sd_s3)
      } else {
        S3 <- rep(NA_real_, block)
      }
      Re405 <- 1 / (1 / G2 - 1 / G5) / 1000
      Re95 <- 1 / (1 / G4 - 1 / G5) / 1000
      Re124 <- 1 / (1 / G4 - 1 / G2) / 1000
      Re99 <- 1 / (1 / G3 - 1 / G5) / 1000
      Re130 <- 1 / (1 / G3 - 1 / G2) / 1000
      Rp1 <- 1 / (1 / p + 1 / G5) / 1000
      Rp2 <- 1 / (1 / p + 1 / G2) / 1000
      Rp3 <- 1 / (1 / p + 1 / G4) / 1000
      Rp4 <- 1 / (1 / p + 1 / G3) / 1000
      if (useT){
        Ro3 <- 1 / (1 / p - 1 / S3) / 1000
      } else {
        Ro3 <- rep(NA_real_, block)
      }
      SE[1, ] <- c(min(SE[1, 1], G2), max(SE[1, 2], G2))
      SE[2, ] <- c(min(SE[2, 1], G3), max(SE[2, 2], G3))
      SE[3, ] <- c(min(SE[3, 1], G4), max(SE[3, 2], G4))
      SE[4, ] <- c(min(SE[4, 1], G5), max(SE[4, 2], G5))
      if (useT){
        SE[5, ] <- c(min(SE[5, 1], S3), max(SE[5, 2], S3))
      }
      ok <- low_405 < Re405 * SR & Re405 * SR < high_405 &
        low_95 < Re95 * SR & low_95 < Re99 * SR & Re95 * SR < high_95 & Re99 * SR < high_95 &
        low_125 < Re124 * SR & low_125 < Re130 * SR & Re124 * SR < high_125 & Re130 * SR < high_125
      if (useT){
        ok <- ok & o1_min > Ro3 * SR & Ro3 * SR > o1_max
      }
      if (useP){
        ok <- ok & p3_min > Rp3 * SR & p3_min > Rp4 * SR & Rp3 * SR > p3_max &
          Rp4 * SR > p3_max & p2_min > Rp2 * SR & Rp2 * SR > p2_max &
          p1_min > Rp1 * SR & Rp1 * SR > p1_max
      }
      naccept <- sum(ok)
      accepted_total <- accepted_total + naccept
      idx <- which(ok)
      if (!length(idx)){
        next
      }
      idx <- idx[seq_len(min(length(idx), nsim - i + 1))]
      rows <- i:(i + length(idx) - 1)
      keep[rows, 1] <- SR[idx]
      keep[rows, 2:6] <- cbind(Re405, Re124, Re130, Re95, Re99)[idx, ]
      if (useT){
        keep[rows, 7] <- Ro3[idx]
      }
      if (useP){
        keep[rows, 9:12] <- cbind(Rp1, Rp2, Rp3, Rp4)[idx, ]
      }
      keep[rows, 14:17] <- cbind(G2, G3, G4, G5)[idx, ]
      if (useT){
        keep[rows, 18] <- S3[idx]
      }
      keep[rows, 19] <- p[idx]
      i <- i + length(idx)
      if (verbose){
        message(sprintf("CTest %s: %d of %d successful samples (%.3f%% accepted)",
                        mode, i - 1, nsim, 100 * accepted_total / attempts))
      }
    }
    i <- nsim + 1
  }

  while (i <= nsim){
    attempts <- attempts + 1
    if (attempts > max_attempts){
      stop("No consistent samples found; check the interpreted cycle ranges.")
    }
    SR <- stats::runif(1, floor(min(SRll, SRul) * 1000), ceiling(max(SRul, SRll) * 1000)) / 1000
    t <- stats::runif(1, tmin, tmax)
    p <- stats::runif(1, ptest[1], ptest[2])
    sc <- sqrt(t / 1e8)
    if (sampleorb == "N"){
      G2 <- stats::rnorm(1, g2, g2 * dg2 * sc)
      G3 <- stats::rnorm(1, g3, g3 * dg3 * sc)
      G4 <- stats::rnorm(1, g4, g4 * dg4 * sc)
      G5 <- stats::rnorm(1, g5, g5 * dg5 * sc)
    } else {
      G2 <- stats::runif(1, g2 - g2 * 2 * dg2 * sc, g2 + g2 * 2 * dg2 * sc)
      G3 <- stats::runif(1, g3 - g3 * 2 * dg3 * sc, g3 + g3 * 2 * dg3 * sc)
      G4 <- stats::runif(1, g4 - g4 * 2 * dg4 * sc, g4 + g4 * 2 * dg4 * sc)
      G5 <- stats::runif(1, g5 - g5 * 2 * dg5 * sc, g5 + g5 * 2 * dg5 * sc)
    }
    if (s3_uncertainty == "scaled"){
      sd_s3 <- s3 * ds3 * sc
    } else {
      sd_s3 <- ds3 * sc
    }
    if (useT){
      S3 <- stats::rnorm(1, s3, sd_s3)
    } else {
      S3 <- NA_real_
    }

    # eccentricity, obliquity and precession periods (kyr) for this sample
    Re405 <- 1 / (1 / G2 - 1 / G5) / 1000
    Re95 <- 1 / (1 / G4 - 1 / G5) / 1000
    Re124 <- 1 / (1 / G4 - 1 / G2) / 1000
    Re99 <- 1 / (1 / G3 - 1 / G5) / 1000
    Re130 <- 1 / (1 / G3 - 1 / G2) / 1000
    Rp1 <- 1 / (1 / p + 1 / G5) / 1000     # k+g5, longest precession (row 6)
    Rp2 <- 1 / (1 / p + 1 / G2) / 1000     # k+g2
    Rp3 <- 1 / (1 / p + 1 / G4) / 1000     # k+g4
    Rp4 <- 1 / (1 / p + 1 / G3) / 1000     # k+g3
    if (useT){
      Ro3 <- 1 / (1 / p - 1 / S3) / 1000
    } else {
      Ro3 <- NA_real_
    }

    SE[1, ] <- c(min(SE[1, 1], G2), max(SE[1, 2], G2))
    SE[2, ] <- c(min(SE[2, 1], G3), max(SE[2, 2], G3))
    SE[3, ] <- c(min(SE[3, 1], G4), max(SE[3, 2], G4))
    SE[4, ] <- c(min(SE[4, 1], G5), max(SE[4, 2], G5))
    if (useT){
      SE[5, ] <- c(min(SE[5, 1], S3), max(SE[5, 2], S3))
    }

    ok_ecc <- low_405 < Re405 * SR && Re405 * SR < high_405 &&
      low_95 < Re95 * SR && low_95 < Re99 * SR && Re95 * SR < high_95 && Re99 * SR < high_95 &&
      low_125 < Re124 * SR && low_125 < Re130 * SR && Re124 * SR < high_125 && Re130 * SR < high_125
    ok_obl <- !useT || (o1_min > Ro3 * SR && Ro3 * SR > o1_max)
    ok_pre <- !useP || (p3_min > Rp3 * SR && p3_min > Rp4 * SR && Rp3 * SR > p3_max &&
                          Rp4 * SR > p3_max && p2_min > Rp2 * SR && Rp2 * SR > p2_max &&
                          p1_min > Rp1 * SR && Rp1 * SR > p1_max)
    if (!(ok_ecc && ok_obl && ok_pre)){
      next
    }

    keep[i, 1] <- SR
    keep[i, 2:6] <- c(Re405, Re124, Re130, Re95, Re99)
    if (useT){
      keep[i, 7] <- Ro3
    }
    if (useP){
      keep[i, 9:12] <- c(Rp1, Rp2, Rp3, Rp4)
    }
    keep[i, 14:17] <- c(G2, G3, G4, G5)
    if (useT){
      keep[i, 18] <- S3
    }
    keep[i, 19] <- p
    if (verbose && i %% report == 0){
      message(sprintf("CTest %s: %d of %d successful samples (%.2f%% accepted)",
                      mode, i, nsim, 100 * i / attempts))
    }
    i <- i + 1
  }

  rng <- function(x){
    c(range(x), mean(x), stats::sd(x))
  }
  na4 <- rep(NA_real_, 4)
  row_obl <- na4
  row_s3 <- na4
  if (useT){
    row_obl <- rng(keep[, 7])
    row_s3 <- rng(keep[, 18])
  }
  row_kg5 <- na4
  row_kg2 <- na4
  row_kg4 <- na4
  row_kg3 <- na4
  if (useP){
    row_kg5 <- rng(keep[, 9])
    row_kg2 <- rng(keep[, 10])
    row_kg4 <- rng(keep[, 11])
    row_kg3 <- rng(keep[, 12])
  }
  rows <- list(
    "sedimentation rate (m/kyr)" = rng(keep[, 1]),
    "405 kyr ecc"               = rng(keep[, 2]),
    "125 kyr ecc"               = rng(keep[, 3:4]),
    "95 kyr ecc"                = rng(keep[, 5:6]),
    "obliquity (k+s3)"          = row_obl,
    "precession k+g5"           = row_kg5,
    "precession k+g2"           = row_kg2,
    "precession k+g4"           = row_kg4,
    "precession k+g3"           = row_kg3,
    "g2 (yr)" = rng(keep[, 14]), "g3 (yr)" = rng(keep[, 15]),
    "g4 (yr)" = rng(keep[, 16]), "g5 (yr)" = rng(keep[, 17]),
    "s3 (yr)" = row_s3,
    "precession constant k (yr)" = rng(keep[, 19]))
  summ <- as.data.frame(do.call(rbind, rows))
  names(summ) <- c("lower", "upper", "mean", "sd")
  summ <- cbind(quantity = rownames(summ), summ)
  rownames(summ) <- NULL

  samples <- as.data.frame(keep[, c(1:7, 9:12, 14:19)])
  names(samples) <- c("SR", "ecc405", "ecc124", "ecc130", "ecc95", "ecc99",
                      "obl_ks3", "prec_kg5", "prec_kg2", "prec_kg4", "prec_kg3",
                      "g2", "g3", "g4", "g5", "s3", "k")
  if (engine == "vectorized"){
    acceptance <- accepted_total / attempts
  } else {
    acceptance <- nsim / attempts
  }
  out <- list(summary = summ, samples = samples, raw = keep, sampled_range = SE,
              mode = mode, nsim = nsim, attempts = attempts,
              acceptance = acceptance, ptest = ptest, sampleorb = sampleorb,
              age_ma = c(d(1, 1), d(1, 2)))
  class(out) <- "ctest"
  out
}

#' @title Print a CTest result
#' @param x a ctest object
#' @param ... ignored
#' @return x, invisibly
#' @export
print.ctest <- function(x, ...){
  cat(sprintf("CTest (%s), %d samples, %.3f%% of draws accepted, age %g-%g Ma\n",
              x$mode, x$nsim, 100 * x$acceptance, x$age_ma[1], x$age_ma[2]))
  k <- x$samples$k
  cat(sprintf("Precession constant k = %.0f yr (%.2f kyr), sd %.0f yr, 95%% CI [%.0f, %.0f]\n",
              mean(k), mean(k) / 1000, stats::sd(k),
              stats::quantile(k, 0.025), stats::quantile(k, 0.975)))
  print(x$summary, digits = 5)
  invisible(x)
}

#' @title Plot the CTest posterior distributions
#' @param x a ctest object
#' @param ... ignored
#' @export
plot.ctest <- function(x, ...){
  op <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(op))
  graphics::plot(stats::density(x$samples$k), lwd = 2, main = "Precession constant",
                 xlab = "k (yr)", ylab = "density")
  graphics::abline(v = mean(x$samples$k), col = "red", lty = 2)
  graphics::plot(stats::density(x$samples$SR), lwd = 2, main = "Sedimentation rate",
                 xlab = "m/kyr", ylab = "density")
  if (all(is.finite(x$samples$obl_ks3))){
    graphics::plot(stats::density(x$samples$obl_ks3), lwd = 2,
                   main = "Obliquity (k+s3)", xlab = "kyr", ylab = "density")
  }
  if (all(is.finite(x$samples$prec_kg5))){
    graphics::plot(stats::density(x$samples$prec_kg5), lwd = 2, col = "red",
                   xlim = range(c(x$samples$prec_kg5, x$samples$prec_kg3)),
                   main = "Precession components", xlab = "kyr", ylab = "density")
    graphics::lines(stats::density(x$samples$prec_kg2), col = "black")
    graphics::lines(stats::density(x$samples$prec_kg4), col = "green4")
    graphics::lines(stats::density(x$samples$prec_kg3), col = "blue")
    graphics::legend("topright", c("k+g5", "k+g2", "k+g4", "k+g3"), bty = "n",
                     col = c("red", "black", "green4", "blue"), lty = 1, cex = 0.8)
  }
}
