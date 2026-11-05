## Origin: M. Arts and A.-C. Da Silva, test_tuned (version 3), Pridoli analysis. Released under the MIT licence.

#' @title Precession-constant test for a tuned record
#'
#' @description Tests a tuned (time-domain) record with the precession constant k.
#' Significant multi-taper (MTM) spectral peaks of the record, optionally pooled with
#' other proxies on the same time scale, are assigned to the obliquity (k+s3, k+s6) and
#' precession (k+g5, k+g2, k+g4) terms. Every combination of candidate peaks is tried
#' over a wide range of k (no window is centred on a reference value), and the
#' combinations in which all terms give one consistent k are kept. Because the terms
#' share one k, a correct assignment gives the same k from every term. The result is
#' compared with AstroGeo22 (Farhat et al. 2022) and a Devonian range, and combined with
#' two further checks of the tuned record: the spectral peaks of the eccentricity,
#' obliquity and precession terms, and the amplitude modulation (405 kyr from the
#' 100 kyr envelope, 100 kyr from the precession envelope) against phase-randomised
#' surrogates. This is the updated version of the test (version 3), which works better
#' for Silurian/Devonian data than a search around the AstroGeo22 value.
#'
#' @param URAN_time Tuned record to test (required, no default): data frame or matrix with time in kyr in column 1 and
#' the proxy value in column 2, e.g. the output of a tuning function. It is sorted, linearly interpolated to 1 kyr,
#' detrended and standardised before the analysis.
#' @param others Other tuned proxies on the same time scale (list of data frames/matrices with the same two-column
#' layout as \code{URAN_time}, or \code{NULL}); they are cut to the time range of \code{URAN_time} and their peaks are
#' pooled with those of \code{URAN_time} (a term found in several proxies adds support). \code{Default=NULL}.
#' @param tbw Time-bandwidth product of the MTM spectrum (numeric, > 0; \code{tbw} of \code{astrochron::mtm}). It sets
#' the resolution: the precession terms k+g5 and k+g2 only separate when \code{tbw <= (record length in kyr) / 1600}
#' (about 3 for 5.5 Myr). A larger \code{tbw} merges them into one broad peak and k cannot be measured (a warning is
#' printed). \code{Default=2}.
#' @param pad Zero-padding factor of the MTM spectrum (numeric, >= 1; \code{padfac} of \code{astrochron::mtm}).
#' \code{Default=5}.
#' @param siglevel Confidence level a spectral peak must reach to be used (numeric, between 0 and 1): the AR1
#' confidence level at the peak (and, with \code{ftest = TRUE}, the harmonic F-test level). Only peaks with periods
#' between 12 and 40 kyr are considered. \code{Default=0.9}.
#' @param ftest Also require the harmonic F-test confidence level to reach \code{siglevel} (logical). The default
#' \code{FALSE} uses the AR1 level only, because the F-test tests for a pure sinusoid and rejects broad, strong peaks
#' of a tuned record. \code{Default=FALSE}.
#' @param tol Half-width of the search window around each predicted term period, as a fraction of the period
#' (numeric, e.g. \code{0.04} = +/- 4\%). \code{Default=0.04}.
#' @param scan Remove intervals with little precession imprint (set them to zero, with a 10-kyr cosine taper at the
#' edges) before the spectrum is built (logical); the pre-scan is controlled by \code{scan_band}, \code{scan_win},
#' \code{scan_thr} and \code{scan_min}. \code{Default=TRUE}.
#' @param scan_band Precession band of the pre-scan, in kyr (numeric vector \code{c(shortest, longest)} period); the
#' record is band-passed (\code{astrochron::taner}) in this band and its Hilbert envelope gives the band variance.
#' \code{Default=c(15, 23)}.
#' @param scan_win Length of the running window over which the precession-band variance is averaged, in kyr (numeric,
#' > 0; at most a third of the record is used). \code{Default=200}.
#' @param scan_thr An interval is weak when its averaged precession-band variance is below \code{scan_thr} times its
#' median (numeric, >= 0). \code{Default=0.5}.
#' @param scan_min Shortest piece kept or removed by the pre-scan, in kyr (numeric, >= 0): weak pieces shorter than
#' this are kept and kept pieces shorter than this are removed. \code{Default=100}.
#' @param n_surr Number of phase-randomised surrogates for the amplitude-modulation test (integer, > 0).
#' \code{Default=200}.
#' @param seed Seed of the random number generator (integer), so that the surrogates can be reproduced.
#' \code{Default=42}.
#' @param plot Draw the comparison figure with \code{\link{plot_test_k}} (logical). \code{Default=TRUE}.
#' @param verbose Also print the peak tables and all peak assignments of the precession-constant test (logical, TRUE
#' or FALSE); the compact report is always printed. \code{Default=FALSE}.
#' @param age_ma Age of the record in Ma (numeric); the AstroGeo22 reference k and its envelope are taken at this age
#' (\code{\link{astrogeo22}}) unless \code{k_ref} and \code{k_ref_range} are given. \code{Default=386}.
#' @param k_ref Reference precession constant, as a period in kyr (numeric); used for the resolution check of
#' \code{tbw}, the old search windows (\code{k_search = NULL}) and the plot. \code{NULL} uses AstroGeo22 at
#' \code{age_ma}. \code{Default=NULL}.
#' @param k_ref_range Range of the reference precession constant in kyr (numeric vector \code{c(lower, upper)}); the test
#' passes when the 1-sigma interval of the estimate overlaps it (or the 20.70 kyr composite of Zeeden et al., 2023).
#' \code{NULL} uses the AstroGeo22 envelope at \code{age_ma}. \code{Default=NULL}.
#' @param k_dev Plausible range of k for the test verdict, in kyr (numeric vector \code{c(lower, upper)}); the default is
#' the Devonian range. \code{Default=c(20.4, 21.6)}.
#' @param k_search Range of k in which peaks are searched, in kyr (numeric vector \code{c(lower, upper)}), independent of
#' AstroGeo22. The search window of each term runs from its period at \code{k_search[1]} to its period at
#' \code{k_search[2]} (plus or minus \code{tol}). \code{NULL} gives the old windows around the AstroGeo22 period, which
#' can only return k within about 1 kyr of AstroGeo22 and so cannot test it. \code{Default=c(19.5, 23.5)}.
#' @param constraints Path to a CSV file of published constraints (character, or \code{NULL} for none), with columns
#' \code{set}, \code{observable}, \code{age_ma}, \code{a_RE}, \code{a_sd}, \code{lod_h}, \code{lod_sd}, \code{k_arcsec}
#' and \code{k_sd}, drawn in panels E-G of the figure. \code{Default=NULL}.
#'
#' @return Invisibly a list with the combined k (\code{k}, \code{k_obliquity}, \code{k_precession},
#' \code{k_ranges}), the peak assignment tables (\code{k_terms}, \code{solutions}), the verdicts and
#' score of the three tests, and the spectra used for plotting.
#'
#' @references
#' Zeeden et al. (2023), see \code{\link{zeeden2023_table1}}.
#' Farhat, M., Auclair-Desrotour, P., Boue, G., Laskar, J. (2022). The resonant tidal evolution of the
#' Earth-Moon distance. Astronomy & Astrophysics 665, L1.
#'
#' @examples
#' \donttest{
#' # synthetic tuned record with a known precession constant (k = 21 kyr)
#' test_tuned_check(k_true = 21, tbw = 2, len = 3000)
#' }
#' @importFrom astrochron taner hilbert detrend
#' @importFrom stats approx cor filter median pchisq runif sd arima.sim
#' @importFrom graphics abline arrows axis layout legend lines mtext par plot.new plot.window points polygon rect segments text
#' @export
test_tuned <- function(URAN_time, others = NULL, tbw = 2, pad = 5, siglevel = 0.9, ftest = FALSE, tol = 0.04,
    scan = TRUE, scan_band = c(15, 23), scan_win = 200, scan_thr = 0.5, scan_min = 100,
    n_surr = 200, seed = 42, plot = TRUE, verbose = FALSE,
    age_ma = 386, k_ref = NULL, k_ref_range = NULL, k_dev = c(20.4, 21.6), k_search = c(19.5, 23.5),
    constraints = NULL) {
  # settings
  ag <- astrogeo22(age_ma)                                      # AstroGeo22 at the age of the record (Farhat et al. 2022)
  if (is.null(k_ref)) {
    k_ref <- 1296 / ag$p_arcsec
  }
  if (is.null(k_ref_range)) {
    k_ref_range <- 1296 / c(ag$p_max_arcsec, ag$p_min_arcsec)
  }
  # k_dev: plausible range (Devonian default: AstroGeo22 envelope + composite of Zeeden et al. (2023; 20.70))
  fr <- c("k+s3" = -18.85, "k+s6" = -26.348, "k+g5" = 4.257452, "k+g2" = 7.452, "k+g4" = 17.916)
  amp <- c("k+s3" = 1, "k+s6" = 0.3, "k+g5" = 1, "k+g2" = 0.87, "k+g4" = 0.7) # relative amplitudes (La2004), for the blended k+g5/g2 term
  ob <- names(fr) %in% c("k+s3", "k+s6")
  targets <- c(405.7, 130.7, 123.8, 98.9, 94.9, 30.6, 19.8, 18.9, 16.4)
  names(targets) <- c("405", "130.7", "123.8", "98.9", "94.9", "obl 30.6", "prec 19.8", "prec 18.9", "prec 16.4")
  prep <- function(x) {
    x <- as.data.frame(x)[, 1:2]
    x <- x[order(x[, 1]), ]
    x <- linterp(x, dt = 1, genplot = FALSE, verbose = FALSE)
    x <- astrochron::detrend(x, genplot = FALSE, verbose = FALSE)
    x[, 2] <- as.numeric(scale(x[, 2]))
    x
  }
  tf <- function(x, lo, hi) astrochron::taner(x, flow = 1 / hi, fhigh = 1 / lo, roll = 10^10, demean = TRUE, detrend = TRUE,
      addmean = FALSE, genplot = FALSE, verbose = FALSE)
  term_period <- function(k, t) 1296 / (1296 / k + fr[[t]])
  k_from <- function(P, t) 1296 / (1296 / P - fr[[t]])
  u <- prep(URAN_time)
  dur <- diff(range(u[, 1]))
  set.seed(seed)

  # 1. precession constant
  recs <- c(list(U = URAN_time), others)
  win <- range(u[, 1])
  for (ir in seq_along(recs)) {
    xr <- prep(recs[[ir]])
    recs[[ir]] <- xr[xr[, 1] >= win[1] & xr[, 1] <= win[2], ]
  }
  hb <- tbw / diff(win) # MTM half-bandwidth (cycles/kyr)
  # largest tbw that still resolves k+g5 and k+g2: full bandwidth 2W <= half their separation (synthetic tests with
  # test_tuned_check(): tbw <= 3 recovers k within ~0.05 kyr for a 5.5-Myr record; tbw 4 already gives errors of ~0.5 kyr, tbw 10-20 fails)
  sep52 <- 1 / term_period(k_ref, "k+g2") - 1 / term_period(k_ref, "k+g5")
  tbw_max <- max(1, floor(sep52 * diff(win) / 4))
  # MTM power spectrum with AR1 background and confidence levels (astrochron::mtm, output = 1)
  spec <- function(x) mtm(x, tbw = tbw, padfac = pad, ar1 = TRUE, output = 1, genplot = FALSE, verbose = FALSE)
  # pre-scan: running precession-band variance of the standardised record (taner bandpass, Hilbert envelope; envelope^2 / 2
  # is the band variance), averaged over scan_win; intervals below scan_thr x its median are weak. Weak pieces shorter than
  # scan_min are kept, kept pieces shorter than scan_min are removed. Removed intervals are set to zero (with a 10-kyr cosine
  # taper at the edges), so the spectrum is built from the good parts only, at their true times (phase kept).
  rmean <- function(v, w) {
    r <- as.numeric(stats::filter(v, rep(1 / w, w), sides = 2))
    r[is.na(r)] <- approx(which(!is.na(r)), r[!is.na(r)], which(is.na(r)), rule = 2)$y
    r
  }
  prescan <- function(x) {
    n <- nrow(x)
    w <- min(scan_win, floor(n / 3))
    env <- astrochron::hilbert(tf(x, scan_band[1], scan_band[2]), genplot = FALSE, verbose = FALSE)[, 2]
    share <- rmean(env^2 / 2, w)
    good <- share >= scan_thr * median(share)
    for (pass in 1:2) {
      r <- rle(good)
      short <- r$lengths < scan_min
      if (pass == 1) {
        r$values[short & !r$values] <- TRUE
      } else {
        r$values[short & r$values] <- FALSE
      }
      good <- inverse.rle(r)
    }
    if (!any(good)) {
      good[] <- TRUE
    }
    data.frame(time = x[, 1], share = share, good = good)
  }
  mask <- function(x, g) {
    y <- x[, 2]
    y <- (y - mean(y[g$good])) / sd(y[g$good])
    wt <- as.numeric(g$good)
    d <- rmean(wt, 10)
    wt <- ifelse(g$good, pmin(1, 2 * d), 0)
    wt <- 0.5 - 0.5 * cos(pi * wt) # cosine taper at the edges
    cbind(x[, 1], y * wt)
  }
  scans <- vector("list", length(recs))
  names(scans) <- names(recs)
  for (ir in seq_along(recs)) {
    scans[[ir]] <- if (scan) { prescan(recs[[ir]]) } else { data.frame(time = recs[[ir]][, 1], share = NA, good = TRUE) }
  }
  rec_nm <- names(recs)
  masked <- vector("list", length(rec_nm))
  for (ir in seq_along(rec_nm)) {
    masked[[ir]] <- mask(recs[[rec_nm[ir]]], scans[[rec_nm[ir]]])
  }
  names(masked) <- names(recs)
  spl <- vector("list", length(masked))
  names(spl) <- names(masked)
  for (ir in seq_along(masked)) {
    spl[[ir]] <- spec(masked[[ir]])
  }
  # significant peaks from the astrochron::mtm spectrum (output = 1): AR1 CL >= siglevel and, with ftest = TRUE, harmonic
  # F-test CL >= siglevel within +/- half the bandwidth (Meyers 2012). Peaks are maxima of the spectrum relative to the AR1
  # background (power / AR1 fit, as plotted) that are the highest point within +/- half the bandwidth (the MTM resolution:
  # wiggles closer than that are not separate peaks). The peak period is the centre of the peak top (power / AR1 fit within
  # 5% of the maximum, power-weighted): the apex for a sharp peak, the middle of a flat-topped MTM peak (large tbw), so k is
  # not pushed to the edge of a flat top by the slope of the AR1 background.
  # (astrochron output = 3 returns the F-test maxima, which often sit on the flank of the power peak, not on its top.)
  # strength = power / AR1 fit at the peak
  peaks_of <- function(nm) {
    m <- spl[[nm]]
    f <- m$Frequency
    n <- length(f)
    df <- f[2] - f[1]
    R <- m$Power / m$AR1_fit
    w <- max(1, round(hb / df))
    apex <- logical(n)
    hcl <- numeric(n)
    for (ii in seq_len(n)) {
      jw <- max(1, ii - w):min(n, ii + w)
      apex[ii] <- R[ii] == max(R[jw]) && ii > 1 && ii < n
      hcl[ii] <- max(m$Harmonic_CL[jw])
    }
    i <- which(apex & m$AR1_CL >= 100 * siglevel & (!ftest | hcl >= 100 * siglevel) & f > 1 / 40 & f < 1 / 12)
    if (!length(i)) {
      return(NULL)
    }
    fc <- numeric(length(i))
    for (ii in seq_along(i)) {
      a <- i[ii]
      lo <- a
      hi <- a
      while (lo > max(1, a - w) && R[lo - 1] >= 0.95 * R[a]) {
        lo <- lo - 1
      }
      while (hi < min(n, a + w) && R[hi + 1] >= 0.95 * R[a]) {
        hi <- hi + 1
      }
      fc[ii] <- sum(f[lo:hi] * R[lo:hi]) / sum(R[lo:hi])
    }
    data.frame(record = nm, period = 1 / fc, harmonic_cl = hcl[i], ar1_cl = m$AR1_CL[i], strength = R[i])
  }
  pk_l <- vector("list", length(rec_nm))
  for (ir in seq_along(rec_nm)) {
    pk_l[ir] <- list(peaks_of(rec_nm[ir]))
  }
  pk <- do.call(rbind, pk_l)
  # peaks of several records closer than the half bandwidth are one spectral peak: merged at their power-weighted mean
  # frequency, strength = summed power / AR1 of the records that show it (a term found in several proxies adds support)
  if (!is.null(pk) && length(recs) > 1) {
    pk <- pk[order(1 / pk$period), ]
    cl <- cumsum(c(1, diff(1 / pk$period) > hb))
    pk_all <- pk
    grp <- split(pk, cl)
    mrg <- vector("list", length(grp))
    names(mrg) <- names(grp)
    for (ig in seq_along(grp)) {
      d <- grp[[ig]]
      j <- which.max(d$strength)
      mrg[[ig]] <- data.frame(record = d$record[j], period = 1 / sum(d$strength / d$period) * sum(d$strength), harmonic_cl = max(d$harmonic_cl),
          ar1_cl = max(d$ar1_cl), strength = sum(d$strength), n_records = length(unique(d$record)))
    }
    pk <- do.call(rbind, mrg)
  }
  use <- pk
  if (is.null(use)) {
    use <- data.frame(record = character(0), period = numeric(0), harmonic_cl = numeric(0), ar1_cl = numeric(0), strength = numeric(0))
  }
  # term sets: k+g5 and k+g2 separate, and (when the bandwidth cannot separate them) one blended term at their
  # amplitude-weighted mean frequency; the set that explains more peaks is used
  g52 <- (amp[["k+g5"]] * fr[["k+g5"]] + amp[["k+g2"]] * fr[["k+g2"]]) / (amp[["k+g5"]] + amp[["k+g2"]])
  sets <- list(fr)
  if (2 * hb > 1 / term_period(k_ref, "k+g2") - 1 / term_period(k_ref, "k+g5")) {
    sets[[2]] <- c(fr[c("k+s3", "k+s6")], "k+g5/g2" = g52, fr["k+g4"])
  }
  tp <- function(k, o) 1296 / (1296 / k + o)
  kf <- function(P, o) 1296 / (1296 / P - o)
  ke <- function(P, o) {
    d <- P^2 * hb / 2
    abs(kf(P + d, o) - kf(P - d, o)) / 2
  }
  is_obl <- function(n) n %in% c("k+s3", "k+s6")
  # Peak assignment and ambiguity. The search windows of neighbouring terms overlap (k+g5 / k+g2) and a window often holds
  # several significant peaks, so "the strongest peak per window" is ambiguous. All terms share one k, so a correct
  # assignment gives the same k from every term (equivalently: it reproduces the k-independent frequency spacing of the
  # terms, e.g. f(k+g2) - f(k+g5) = (g2 - g5) / 1296 cycles/kyr).
  # 1. Every combination of candidate peaks (or no peak) per term is tried, keeping the period order of the terms and using
  #    each peak once. A combination is consistent when the k values of its terms agree (chi-square p >= 0.05; error per
  #    term = peak position (quarter bandwidth) + assumed uncertainty of the g/s frequency in the Devonian, sig_g).
  # 2. Consistent combinations are grouped into options with distinct k (> 2 sigma apart); per option the combination with
  #    the most terms, then the strongest peaks, is kept. Precession options need >= 2 terms, obliquity options >= 1 (k+s6
  #    is weak, about 0.3 x k+s3, and often absent).
  # 3. Precession and obliquity are independent: each precession option is combined with the obliquity option that agrees
  #    with it (within 2 sigma; most terms, then strongest). The solutions are ranked by the number of terms that agree,
  #    then by summed peak strength (power / AR1).
  # 4. The best solution is unambiguous only when it has more agreeing terms than any other solution could reach if
  #    obliquity options up to 3 sigma away were also accepted. All solutions are combined into one k +/- sd per group
  #    (obliquity, precession, obliquity + precession), see below.
  sig_g <- c("k+s3" = 0.2, "k+s6" = 0.05, "k+g5" = 0.01, "k+g2" = 0.05, "k+g4" = 0.2, "k+g5/g2" = 0.05) # "/yr, assumed
  ke2 <- function(P, o, t) sqrt(ke(P, o)^2 + (kf(P, o)^2 / 1296 * sig_g[[t]])^2)
  row_of <- function(t, o, b) if (is.null(b)) {
    data.frame(term = t, peak = NA, record = NA, cl = NA, strength = NA, k = NA, k_err = NA)
  } else {
    data.frame(term = t, peak = b$period, record = b$record, cl = b$harmonic_cl, strength = b$strength, k = kf(b$period, o), k_err = ke2(b$period, o, t))
  }
  est_k <- function(a) {
    q <- !is.na(a$k)
    if (!any(q)) {
      return(c(k = NA, sd = NA, se = NA))
    }
    w <- 1 / a$k_err[q]^2
    c(k = sum(w * a$k[q]) / sum(w), sd = max(sqrt(1 / sum(w)), if (sum(q) > 1) { sd(a$k[q]) } else { 0 }), se = sqrt(1 / sum(w)))
  }
  chi_p <- function(a) {
    q <- !is.na(a$k)
    if (sum(q) < 2) {
      return(1)
    }
    w <- 1 / a$k_err[q]^2
    m <- sum(w * a$k[q]) / sum(w)
    pchisq(sum(w * (a$k[q] - m)^2), sum(q) - 1, lower.tail = FALSE)
  }
  solutions <- function(frs) {
    nm <- names(frs)
    cands <- vector("list", length(nm))
    for (it in seq_along(nm)) {
      tm <- nm[it]
      cands[[it]] <- c(NA, if (is.null(k_search)) {
        which(abs(use$period / tp(k_ref, frs[[tm]]) - 1) <= tol)
      } else {
        which(use$period >= tp(k_search[1], frs[[tm]]) * (1 - tol) & use$period <= tp(k_search[2], frs[[tm]]) * (1 + tol))
      })
    }
    g <- as.matrix(expand.grid(cands))
    out <- list()
    for (r in seq_len(nrow(g))) {
      id <- g[r, ]
      q <- !is.na(id)
      po <- id[order(-tp(k_ref, unlist(frs)))]
      po <- po[!is.na(po)] # peaks in the order of the predicted periods
      if (!any(q) || any(duplicated(id[q])) || is.unsorted(-use$period[po], strictly = TRUE)) {
        next
      }
      rows <- vector("list", length(nm))
      for (j in seq_along(nm)) {
        rows[[j]] <- row_of(nm[j], frs[[nm[j]]], if (is.na(id[j])) { NULL } else { use[id[j], ] })
      }
      a <- do.call(rbind, rows)
      a$AstroGeo22 <- tp(k_ref, frs[a$term])
      out[[length(out) + 1]] <- list(table = a, k = est_k(a), n = sum(q), p = chi_p(a), str = sum(a$strength, na.rm = TRUE))
    }
    out
  }
  empty_tab <- function(sub) {
    nm <- names(fr)[sub(names(fr))]
    rows <- vector("list", length(nm))
    for (it in seq_along(nm)) {
      rows[[it]] <- row_of(nm[it], fr[[nm[it]]], NULL)
    }
    a <- do.call(rbind, rows)
    a$AstroGeo22 <- tp(k_ref, fr[nm])
    a
  }
  # numeric field f of every element of the list l (in place of sapply(l, `[[`, f))
  num_field <- function(l, f) {
    v <- numeric(length(l))
    for (q in seq_along(l)) {
      v[q] <- l[[q]][[f]]
    }
    v
  }
  options_of <- function(sub, min_n) {
    sol <- vector("list", length(sets))
    for (iset in seq_along(sets)) {
      s <- sets[[iset]]
      sol[[iset]] <- solutions(s[sub(names(s))])
    }
    sol <- unlist(sol, recursive = FALSE)
    ok <- logical(length(sol))
    for (q in seq_along(sol)) {
      ok[q] <- sol[[q]]$p >= 0.05 && sol[[q]]$n >= min_n
    }
    sol <- sol[which(ok)]
    if (!length(sol)) {
      return(list())
    }
    sol <- sol[order(-num_field(sol, "n"), -num_field(sol, "str"))]
    keep <- list()
    for (x in sol) {
      apart <- logical(length(keep))
      for (q in seq_along(keep)) {
        b <- keep[[q]]
        apart[q] <- abs(x$k[["k"]] - b$k[["k"]]) > 2 * sqrt(b$k[["se"]]^2 + x$k[["se"]]^2)
      }
      if (!length(keep) || all(apart)) {
        keep[[length(keep) + 1]] <- x
      }
    }
    keep
  }
  dist_k <- function(a, b) abs(a$k[["k"]] - b$k[["k"]]) / sqrt(a$k[["se"]]^2 + b$k[["se"]]^2)
  opt_p <- options_of(function(n) !is_obl(n), 2)
  opt_o <- options_of(is_obl, 1)
  obl_found <- length(opt_o) > 0
  best_of <- function(l) if (length(l)) { l[[order(-num_field(l, "n"), -num_field(l, "str"))[1]]] } else { NULL }
  sol_of <- function(P) {
    d <- numeric(length(opt_o))
    for (q in seq_along(opt_o)) {
      d[q] <- dist_k(opt_o[[q]], a = P)
    }
    O <- best_of(opt_o[d <= 2])
    O3 <- best_of(opt_o[d <= 3])
    tab <- rbind(if (is.null(O)) { empty_tab(is_obl) } else { O$table }, P$table)
    list(prec = P, obl = O, dist = if (is.null(O)) {
      if (length(d)) {
        min(d)
      } else {
        NA
      }
    } else {
      dist_k(P, O)
    }, table = tab,
        k = if (is.null(O)) { P$k } else { est_k(tab) }, n = P$n + if (is.null(O)) { 0 } else { O$n }, n3 = P$n + if (is.null(O3)) { 0 } else { O3$n },
        str = P$str + if (is.null(O)) { 0 } else { O$str })
  }
  sols <- vector("list", length(opt_p))
  for (ip in seq_along(opt_p)) {
    sols[[ip]] <- sol_of(opt_p[[ip]])
  }
  if (length(sols)) {
    sols <- sols[order(-num_field(sols, "n"), -num_field(sols, "str"))]
  }
  unique_best <- length(sols) == 1 || (length(sols) > 1 && sols[[1]]$n > max(num_field(sols[-1], "n3")))
  amb <- if (length(sols) > 1 && !unique_best) { "unresolved" } else if (length(sols) > 1) { "resolved" } else { "none" }
  # best solution (or the strongest obliquity option when no precession option exists)
  S <- if (length(sols)) { sols[[1]] } else { NULL }
  o_best <- if (!is.null(S) && !is.null(S$obl)) { S$obl } else { best_of(opt_o) }
  k_obl <- if (is.null(o_best)) { c(k = NA, sd = NA, se = NA) } else { o_best$k }
  k_prec <- if (is.null(S)) { c(k = NA, sd = NA, se = NA) } else { S$prec$k }
  k_both <- if (!is.null(S) && !is.null(S$obl)) { S$k } else { c(k = NA, sd = NA, se = NA) }
  kt <- if (!is.null(S)) { S$table } else { rbind(if (is.null(o_best)) { empty_tab(is_obl) } else { o_best$table }, empty_tab(function(n) !is_obl(n))) }
  k_ag <- k_both[["k"]]
  k_ag_sd <- k_both[["sd"]]
  dev_ok <- function(x) !is.na(x[1]) && x[1] - 2 * x[2] <= k_dev[2] && x[1] + 2 * x[2] >= k_dev[1]
  prec_ok <- !is.null(S) && dev_ok(k_prec)
  obl_ok <- !is.null(S) && !is.null(S$obl) && dev_ok(k_obl)
  ag_found <- prec_ok
  used <- function(a) {
    a <- a[!is.na(a$k), ]
    if (nrow(a)) {
      paste(sprintf("%s %.1f", a$term, a$peak), collapse = ", ")
    } else {
      "no peaks"
    }
  }
  # k of each solution: obliquity + precession when they agree, otherwise precession alone
  sol_k <- vector("list", length(sols))
  for (i in seq_along(sols)) {
    x <- sols[[i]]
    sol_k[[i]] <- data.frame(option = i, k = x$k[["k"]], sd = x$k[["sd"]], what = if (is.null(x$obl)) { "precession alone" } else { "obliquity + precession" },
        k_prec = x$prec$k[["k"]], sd_prec = x$prec$k[["sd"]], k_obl = if (is.null(x$obl)) { NA } else { x$obl$k[["k"]] }, sd_obl = if (is.null(x$obl)) { NA } else { x$obl$k[["sd"]] },
        terms = used(x$table), n = x$n, strength = x$str)
  }
  sol_k <- if (length(sol_k)) { do.call(rbind, sol_k) } else { NULL }
  # One k and uncertainty per group (obliquity, precession, obliquity + precession), combining all solutions of that group
  # (model averaging): each solution i is a normal distribution N(k_i, sd_i) with weight w_i = its summed spectral power
  # (power / AR1 of its peaks, so better-supported assignments count more). The combined k is the weighted mean and the
  # combined sd follows from the law of total variance, sd^2 = sum(w_i * (sd_i^2 + (k_i - k)^2)): it holds both the error
  # of each solution and the spread between the solutions (the peak-assignment ambiguity). Reported as k +/- sd (1 sigma)
  # and the 2-sigma range. With one solution this is that solution's k +/- sd.
  rng <- function(what, l, w) {
    if (!length(l)) {
      return(NULL)
    }
    k <- numeric(length(l))
    sd <- numeric(length(l))
    for (q in seq_along(l)) {
      k[q] <- l[[q]]$k[["k"]]
      sd[q] <- l[[q]]$k[["sd"]]
    }
    w <- w / sum(w)
    m <- sum(w * k)
    s2 <- sqrt(sum(w * (sd^2 + (k - m)^2)))
    data.frame(what = what, k = m, sd = s2, lo = m - 2 * s2, hi = m + 2 * s2, options = length(l), k_best = k[1], sd_best = sd[1], k_min = min(k), k_max = max(k))
  }
  # Only the solutions still in play are combined: the best solution alone when the assignment is resolved (it has more
  # agreeing terms than any other solution could reach), otherwise all solutions that could reach its number of terms.
  # Averaging over solutions that were already rejected would only widen the uncertainty with assignments that do not fit.
  if (amb == "unresolved") {
    ok <- logical(length(sols))
    for (q in seq_along(sols)) {
      ok[q] <- sols[[q]]$n3 >= sols[[1]]$n
    }
    inplay <- sols[which(ok)]
  } else {
    inplay <- sols[seq_len(min(1, length(sols)))]
  }
  if (length(inplay)) {
    obl_l <- list()
    for (x in inplay) {
      if (!is.null(x[["obl"]])) {
        obl_l[[length(obl_l) + 1]] <- x[["obl"]]
      }
    }
  } else if (is.null(o_best)) {
    obl_l <- list()
  } else {
    obl_l <- list(o_best)
  }
  both_l <- list()
  prec_l <- vector("list", length(inplay))
  prec_w <- numeric(length(inplay))
  for (q in seq_along(inplay)) {
    x <- inplay[[q]]
    prec_l[[q]] <- x[["prec"]]
    prec_w[q] <- x$prec$str
    if (!is.null(x$obl)) {
      both_l[[length(both_l) + 1]] <- x
    }
  }
  both_k <- vector("list", length(both_l))
  for (q in seq_along(both_l)) {
    both_k[[q]] <- list(k = both_l[[q]]$k)
  }
  k_ranges <- rbind(rng("obliquity", obl_l, num_field(obl_l, "str")),
      rng("precession", prec_l, prec_w),
      rng("obliquity + precession", both_k, num_field(both_l, "str")))
  k_shown <- if (is.null(k_ranges)) { NA } else if ("obliquity + precession" %in% k_ranges$what) { "obliquity + precession" } else { k_ranges$what[nrow(k_ranges)] }
  # Test 1 verdict on the combined k (obliquity + precession, otherwise precession; peak-assignment ambiguity is part of its
  # sd): PASS when it lies in the plausible range (k_dev) and its 1-sigma interval overlaps the AstroGeo22 envelope or
  # the Devonian composite (20.70). Without a precession solution (>= 2 terms with one consistent k): FAIL.
  kc <- if (is.na(k_shown)) { NULL } else { k_ranges[k_ranges$what == k_shown, ] }
  k_status <- if (is.null(kc) || kc$what == "obliquity") {
    "FAIL"
  } else {
    ok1 <- kc$k >= k_dev[1] && kc$k <= k_dev[2] && ((kc$k - kc$sd <= k_ref_range[2] && kc$k + kc$sd >= k_ref_range[1]) || abs(kc$k - 20.70) <= kc$sd)
    if (!ok1) {
      "FAIL"
    } else if (kc$what == "precession" && !obl_found) {
      "PASS (precession only)"
    } else if (kc$what == "precession") {
      "PASS (obliquity disagrees)"
    } else {
      "PASS"
    }
  }
  k_pass <- grepl("^PASS", k_status)

  # 2. spectral peaks of the whole tuned record
  s <- spec(u)
  sig <- logical(length(targets))
  names(sig) <- names(targets)
  for (it in seq_along(targets)) {
    q <- abs(s[, 1] - 1 / targets[[it]]) <= tbw / dur
    sig[it] <- any(s$Power[q] > s$AR1_95_power[q])
  }
  to <- timeOpt(u, sedmin = 99.9, sedmax = 100.1, numsed = 3, fit = 1, flow = 1 / 22, fhigh = 1 / 15, roll = 10^7,
      targetE = c(405.7, 130.7, 123.8, 98.9, 94.9), targetP = c(19.8, 18.9, 16.4), output = 1, genplot = FALSE, verbose = FALSE)

  # 3. amplitude modulation, against phase-randomised surrogates of the tuned record
  am <- function(x) {
    e405 <- tf(x, 350, 450)
    e100 <- tf(x, 85, 135)
    pr <- tf(x, 15, 22)
    a <- tf(astrochron::hilbert(e100, genplot = FALSE, verbose = FALSE), 350, 450)
    b <- tf(astrochron::hilbert(pr, genplot = FALSE, verbose = FALSE), 85, 135)
    k <- x[, 1] >= min(x[, 1]) + 300 & x[, 1] <= max(x[, 1]) - 300
    c(r405 = cor(e405[k, 2], a[k, 2]), r100 = cor(e100[k, 2], b[k, 2]))
  }
  surr <- function(x) {
    n <- nrow(x)
    X <- fft(x[, 2])
    h <- floor((n - 1) / 2)
    X[2:(h + 1)] <- Mod(X[2:(h + 1)]) * exp(1i * runif(h, 0, 2 * pi))
    X[n:(n - h + 1)] <- Conj(X[2:(h + 1)])
    cbind(x[, 1], Re(fft(X, inverse = TRUE)) / n)
  }
  obs <- am(u)
  nul <- matrix(NA_real_, n_surr, 2, dimnames = list(NULL, c("r405", "r100")))
  for (s in seq_len(n_surr)) {
    nul[s, ] <- am(surr(u))
  }
  p_am <- c(r405 = mean(nul[, 1] >= obs[1]), r100 = mean(nul[, 2] >= obs[2]))

  # report (compact; verbose = TRUE adds the peak tables and all solutions of test 1)
  # Test 2 (spectral peaks): PASS when the 405-kyr peak, >= 2 of the 4 short-eccentricity peaks and >= 2 of the 3 precession
  # peaks rise above AR1 95% (minor terms are often weak in a tuned record, so not all 9 are required).
  # Test 3 (amplitude modulation): PASS p < 0.05; WEAK PASS 0.05 <= p < 0.10 (marginal: the modulation is present but the
  # surrogates cannot fully exclude chance); FAIL p >= 0.10. A weak pass counts as passed but is flagged in the score.
  ecc_ok <- sig[["405"]] && sum(sig[c("130.7", "123.8", "98.9", "94.9")]) >= 2 && sum(sig[c("prec 19.8", "prec 18.9", "prec 16.4")]) >= 2
  am_v <- function(p) if (p < 0.05) { "PASS" } else if (p < 0.10) { "WEAK PASS" } else { "FAIL" }
  yn <- function(x) if (is.na(x)) { "not tested" } else if (x) { "PASS" } else { "FAIL" }
  ks <- function(x) if (is.na(x[1])) { "-" } else { sprintf("%.2f +/- %.2f kyr", x[1], x[2]) }
  kr1 <- function(w) {
    r <- k_ranges[k_ranges$what == w, ]
    if (is.null(r) || !nrow(r)) {
      "-"
    } else {
      sprintf("%.2f +/- %.2f", r$k, r$sd)
    }
  }
  cat("\n=== Age-model test:", round(dur), "kyr ===\n\n")
  cat(sprintf("1. Precession constant (MTM tbw %g; %s)%s -> %s\n", tbw,
      if (scan) { sprintf("%.0f%% kept after pre-scan", 100 * mean(scans[[1]]$good)) } else { "no pre-scan" },
      if (tbw > tbw_max) { sprintf(" WARNING: tbw too large, use <= %d", tbw_max) } else { "" }, k_status))
  if (!is.null(kc)) {
    cat(sprintf("   k = %.2f +/- %.2f kyr (1 sigma; %s%s)\n", kc$k, kc$sd, kc$what,
        if (kc$options > 1) { sprintf(", %d peak assignments combined", kc$options) } else { "" }))
  }
  cat(sprintf("   obliquity %s | precession %s | AstroGeo22 %g Ma %.2f (%.2f-%.2f) | Devonian composite 20.70\n",
      kr1("obliquity"), kr1("precession"), age_ma, k_ref, k_ref_range[1], k_ref_range[2]))
  if (amb == "unresolved") {
    cat("   note: the peak assignment is ambiguous; this is included in the uncertainty (verbose = TRUE lists the solutions)\n")
  }
  if (verbose) {
    cat("   peaks 12-40 kyr:", if (nrow(use)) { paste(sprintf("%.1f", sort(use$period, decreasing = TRUE)), collapse = ", ") } else { "none" }, "\n")
    cat("   best solution: term AstroGeo22 peak F-test CL power/AR1 k from peak\n")
    for (q in seq_len(nrow(kt))) {
      cat(sprintf("   %-7s %5.1f kyr %s %s %s %s\n", kt$term[q], kt$AstroGeo22[q],
          if (is.na(kt$peak[q])) { "    -    " } else { sprintf("%5.1f kyr", kt$peak[q]) }, if (is.na(kt$cl[q])) { "   -  " } else { sprintf("%5.1f%%", kt$cl[q]) },
          if (is.na(kt$strength[q])) { "    -  " } else { sprintf("%7.2f", kt$strength[q]) }, if (is.na(kt$k[q])) { "no peak" } else { sprintf("%.2f +/- %.2f kyr", kt$k[q], kt$k_err[q]) }))
    }
    if (length(sols)) {
      cat("   all solutions: # k precession obliquity terms power/AR1\n")
      for (i in seq_along(sols)) {
        x <- sols[[i]]
        cat(sprintf("   %d %-14s %-35s %-33s %3d %6.1f\n", i, ks(x$k), used(x$prec$table),
            if (is.null(x$obl)) { sprintf("none (nearest %.1f sigma)", x$dist) } else { sprintf("%s (%.1f sigma)", used(x$obl$table), x$dist) }, x$n, x$str))
      }
    }
  }
  cat(sprintf("\n2. Spectral peaks above AR1 95%%: %d of %d (405 %s, short eccentricity %d/4, obliquity %s, precession %d/3) -> %s\n",
      sum(sig), length(sig), if (sig[["405"]]) { "+" } else { "-" }, sum(sig[c("130.7", "123.8", "98.9", "94.9")]), if (sig[["obl 30.6"]]) { "+" } else { "-" },
      sum(sig[c("prec 19.8", "prec 18.9", "prec 16.4")]), yn(ecc_ok)))
  cat(sprintf("   TimeOpt r2 (precession envelope): %.3f\n", to[2, 2]))
  cat(sprintf("\n3. Amplitude modulation (%d surrogates)\n", n_surr))
  cat(sprintf("   405 from 100-kyr envelope:   r = %.2f, p = %.3f -> %s\n", obs[1], p_am[1], am_v(p_am[1])))
  cat(sprintf("   100 from precession envelope: r = %.2f, p = %.3f -> %s\n", obs[2], p_am[2], am_v(p_am[2])))
  verdicts <- c(k_status, yn(ecc_ok), am_v(p_am[1]), am_v(p_am[2]))
  score <- sum(grepl("PASS", verdicts))
  n_weak <- sum(verdicts == "WEAK PASS")
  cat(sprintf("\nScore: %d of %d tests passed%s\n", score, length(verdicts), if (n_weak) { sprintf(" (%d weak)", n_weak) } else { "" }))
  sol_tabs <- vector("list", length(sols))
  for (q in seq_along(sols)) {
    sol_tabs[[q]] <- sols[[q]][["table"]]
  }
  out <- list(duration = dur, k_terms = kt, k = k_both, k_obliquity = k_obl, k_precession = k_prec, k_status = k_status, k_ambiguity = amb,
      solutions = sol_k, solution_tables = sol_tabs, k_ranges = k_ranges, k_shown = k_shown, significant = sig, timeopt_r2env = to[2, 2], am = obs, am_p = p_am,
      score = score, verdicts = verdicts, peaks_test = ecc_ok, mtm = s,
      spec = list(scan = scans, raw = recs, scan_thr = scan_thr, tbw = tbw, tbw_max = tbw_max, pad = pad, siglevel = siglevel, ftest = ftest, peaks = pk, spectra = spl, records = names(recs), window = win, tol = tol, hb = hb,
          age_ma = age_ma, k_ref = k_ref, k_ref_range = k_ref_range, k_search = k_search))
  if (plot) {
    plot_test_k(out, constraints)
  }
  invisible(out)
}


#' @title Comparison figure of the precession-constant test
#'
#' @description Plots the result of \code{\link{test_tuned}}: the MTM spectrum with the AstroGeo22 periods,
#' search windows and the peaks used (A), k, Earth-Moon distance and length of day against AstroGeo22 and
#' the Devonian composite (B-D), AstroGeo22 through time with the estimate and optional published
#' constraints (E-G), and the pre-scan of the record.
#'
#' @param res Result to plot (required, no default): the list returned by \code{\link{test_tuned}}.
#' @param constraints Path to a CSV file of published constraints drawn in panels E-G (character, or \code{NULL} for
#' none); see \code{constraints} in \code{\link{test_tuned}} for the columns. \code{Default=NULL}.
#' @param age_ma Age of the record in Ma (numeric), where the estimate is drawn in panels E-G and at which AstroGeo22 is
#' evaluated. \code{Default=res$spec$age_ma}, the age used in \code{\link{test_tuned}}.
#' @return Called for its plot.
#' @importFrom grDevices adjustcolor hcl.colors
#' @export
plot_test_k <- function(res, constraints = NULL, age_ma = res$spec$age_ma) {
  fr <- c("k+s3" = -18.85, "k+s6" = -26.348, "k+g5" = 4.257452, "k+g2" = 7.452, "k+g4" = 17.916)
  term_period <- function(k, t) 1296 / (1296 / k + fr[[t]])
  k_ag <- res$spec$k_ref
  k_ag_rng <- res$spec$k_ref_range
  k_zee <- 20.70
  col_t <- c("k+s3" = "#7b6fb0", "k+s6" = "#3f3a8c", "k+g5" = "#1b9e77", "k+g2" = "#d95f02", "k+g4" = "#e7298a", "k+g5/g2" = "#7f7f00")
  sp <- res$spec
  tol <- sp$tol
  kt <- res$k_terms
  # one k per group (obliquity, precession, obliquity + precession), combined over all solutions: k +/- 1 sigma (thick) and
  # +/- 2 sigma (thin)
  kr <- res$k_ranges
  cl_est <- c("obliquity" = "grey25", "precession" = "#56b4e9", "obliquity + precession" = "#0072b2")
  em <- function(k) em_state(k)
  conv <- function(k, col) if (col == "k") { k } else { em(k)[[col]] }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
      1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4,
      5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
      8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8), 4, byrow = TRUE), heights = c(0.5, 1, 1, 0.22))
  par(las = 1, tcl = -0.25, mgp = c(1.8, 0.45, 0), cex.axis = 0.8, cex.lab = 0.85)
  # A: spectrum, AstroGeo22 periods (+/- tol), Devonian composite, significant peaks, peaks of the best solution (filled
  # triangles) and of the other solutions (open triangles)
  par(mar = c(3.2, 4, 2.6, 0.8))
  yy <- vector("list", length(sp$spectra))
  names(yy) <- names(sp$spectra)
  for (j in seq_along(sp$spectra)) {
    m <- sp$spectra[[j]]
    f <- m[, 1] > 1 / 40 & m[, 1] < 1 / 12
    yy[[j]] <- log10(m$Power[f] / m$AR1_fit[f])
  }
  yy <- unlist(yy)
  yr <- range(c(yy, 0), na.rm = TRUE)
  yr[2] <- yr[2] + (if (is.null(sp$k_search)) { 0.06 } else { 0.25 }) * diff(yr)
  plot(NA, xlim = c(12, 40), ylim = yr, xlab = "Period (kyr)", ylab = "log10 power / AR1 fit", xaxs = "i")
  ks <- sp$k_search
  dy0 <- diff(yr)
  for (t in names(fr)) {
    pp <- term_period(k_ag, t)
    if (is.null(ks)) {
      rect(pp * (1 - tol), -1e3, pp * (1 + tol), 1e3, col = adjustcolor(col_t[t], 0.18), border = NA)
    } else {
      j <- match(t, names(fr))                                  # search window of each term (k_search), as a bar at the top
      segments(term_period(ks[1], t) * (1 - tol), yr[2] - (0.02 + 0.035 * j) * dy0, term_period(ks[2], t) * (1 + tol),
          yr[2] - (0.02 + 0.035 * j) * dy0, col = col_t[t], lwd = 4, lend = 1)
    }
    abline(v = pp, col = col_t[t], lty = 2)
    abline(v = term_period(k_zee, t), col = "red3", lty = 3)
    mtext(t, 3, 0.1, at = pp, col = col_t[t], cex = 0.65, font = 2)
  }
  cs <- if (length(sp$spectra) > 1) { hcl.colors(length(sp$spectra), "Dark 3") } else { "black" }
  for (j in seq_along(sp$spectra)) {
    m <- sp$spectra[[j]]
    f <- m[, 1] > 1 / 40 & m[, 1] < 1 / 12
    lines(1 / m[f, 1], log10(m$Power[f] / m$AR1_fit[f]), col = cs[j], lwd = 1)
    lines(1 / m[f, 1], log10(m$AR1_90_power[f] / m$AR1_fit[f]), col = "red3", lty = 2, lwd = 0.8)
    lines(1 / m[f, 1], log10(m$AR1_99_power[f] / m$AR1_fit[f]), col = "red3", lty = 3, lwd = 0.8)
  }
  on_curve <- function(P, rec) {
    y <- numeric(length(P))
    for (q in seq_along(P)) {
      m <- sp$spectra[[rec[q]]]
      y[q] <- approx(1 / m[, 1], log10(m$Power / m$AR1_fit), P[q])$y
    }
    y
  }
  if (!is.null(sp$peaks) && nrow(sp$peaks)) {
    points(sp$peaks$period, on_curve(sp$peaks$period, sp$peaks$record), pch = 1, cex = 0.9, col = "grey20")
  }
  dy <- diff(par("usr")[3:4])
  oth <- vector("list", length(res$solution_tables[-1]))
  for (j in seq_along(oth)) {
    a <- res$solution_tables[-1][[j]]
    oth[[j]] <- a[!is.na(a$peak), ]
  }
  oth <- do.call(rbind, oth)
  if (!is.null(oth) && nrow(oth)) {
    oth <- oth[!duplicated(round(oth$peak, 3)) & !(round(oth$peak, 3) %in% round(kt$peak, 3)), ]
    points(oth$peak, on_curve(oth$peak, oth$record) + 0.03 * dy, pch = 25, bg = "white", col = col_t[oth$term], cex = 1.2, lwd = 1.3)
  }
  s0 <- !is.na(kt$peak)
  if (any(s0)) {
    points(kt$peak[s0], on_curve(kt$peak[s0], kt$record[s0]) + 0.03 * dy, pch = 25, bg = col_t[kt$term[s0]], cex = 1.3)
  }
  if (length(sp$spectra) > 1) {
    legend("topright", sp$records, col = cs, lwd = 1, bg = "white", cex = 0.6, ncol = 3)
  }
  too_wide <- sp$tbw > sp$tbw_max
  mtext(sprintf("A  astrochron MTM (tbw %g, %d-%d kyr%s) - test 1: %s%s%s", sp$tbw, round(sp$window[1]), round(sp$window[2]),
      if (is.null(sp$k_search)) { "" } else { sprintf("; bars: search windows for k %g-%g", sp$k_search[1], sp$k_search[2]) }, res$k_status,
      if (!is.null(res$solutions) && nrow(res$solutions) > 1) { sprintf(" (%d solutions%s)", nrow(res$solutions), if (res$k_ambiguity == "unresolved") { ", ambiguous" } else { "" }) } else { "" },
      if (too_wide) { sprintf(" - tbw too large (use tbw <= %d)", sp$tbw_max) } else { "" }), 3, 1.2, adj = 0, font = 2, cex = 0.8, col = if (too_wide) { "red3" } else { "black" })
  # B-D: obliquity, precession and obliquity + precession vs AstroGeo22 (grey band, dashed) and the Devonian composite (red dotted)
  dot <- function(col, ref, ref_rng, zee, xlab, lab) {
    n <- if (is.null(kr)) { 0 } else { nrow(kr) }
    v <- if (n) {
      data.frame(what = kr$what, val = conv(kr$k, col), b1 = conv(kr$k - kr$sd, col), b2 = conv(kr$k + kr$sd, col),
          r1 = conv(kr$lo, col), r2 = conv(kr$hi, col), options = kr$options)
    } else {
      NULL
    }
    xr <- range(c(ref_rng, zee, if (n) { c(v$r1, v$r2) }), na.rm = TRUE)
    xr <- xr + c(-0.12, 0.12) * diff(xr)
    plot(NA, xlim = xr, ylim = c(0.5, 3.5), yaxt = "n", xlab = xlab, ylab = "")
    rect(ref_rng[1], -9, ref_rng[2], 99, col = "grey88", border = NA)
    abline(v = ref, lty = 2, col = "grey30")
    abline(v = zee, lty = 3, col = "red3")
    if (n) {
      y <- c("obliquity" = 3, "precession" = 2, "obliquity + precession" = 1)[v$what]
      cl <- cl_est[v$what]
      segments(v$r1, y, v$r2, y, lwd = 1.2, col = cl)
      segments(v$b1, y, v$b2, y, lwd = 3.5, col = cl)
      points(v$val, y, pch = 16, col = cl, cex = 1.1)
      text(xr[1], y + 0.3, sprintf("%s%s", v$what, ifelse(v$options > 1, sprintf(" (%d solutions combined)", v$options), "")), adj = 0, cex = 0.6)
    }
    mtext(lab, 3, 1.2, adj = 0, font = 2, cex = 0.8)
  }
  par(mar = c(3.2, 0.8, 2.6, 0.8))
  dot("k", k_ag, k_ag_rng, k_zee, "k (kyr)", "B  precession constant")
  ag_age <- astrogeo22(age_ma)
  zee <- em_state(k_zee)
  dot("a_RE", ag_age$a_RE, c(ag_age$a_min_RE, ag_age$a_max_RE), zee$a_RE, "Earth-Moon distance (Earth radii)", "C  Earth-Moon distance")
  dot("lod_solar_h", ag_age$lod_h, c(ag_age$lod_min_h, ag_age$lod_max_h), zee$lod_solar_h, "length of day (h)", "D  length of day")
  # E-G: AstroGeo22 through time, published constraints, this record (obliquity + precession; otherwise the group available)
  ages <- seq(0, 560, by = 2)
  A <- astrogeo22(ages)
  pc <- if (!is.null(constraints) && file.exists(constraints)) { utils::read.csv(constraints) } else { NULL }
  me <- res$k_shown
  b <- if (!is.na(me)) { kr[kr$what == me, ] } else { NULL }
  curve_panel <- function(y, ylo, yhi, py, psd, col, ylab, lab) {
    v <- if (!is.null(b)) { c(val = conv(b$k, col), b1 = conv(b$k - b$sd, col), b2 = conv(b$k + b$sd, col), r1 = conv(b$lo, col), r2 = conv(b$hi, col)) } else { NULL }
    yr <- range(c(y, v), na.rm = TRUE)
    yr <- yr + c(-0.15, 0.15) * diff(yr)
    par(mar = c(3.4, 3.6, 2.6, 0.8))
    plot(NA, xlim = c(0, 580), ylim = yr, xlab = "age (Ma)", ylab = ylab)
    rect(age_ma - 8, -1e3, age_ma + 8, 1e3, col = "grey88", border = NA)
    polygon(c(ages, rev(ages)), c(ylo, rev(yhi)), col = adjustcolor("grey40", 0.25), border = NA)
    lines(ages, y, lwd = 2.2, col = "grey20")
    if (!is.null(pc)) {
      c22 <- pc$set == "cyclo22"
      tid <- grepl("tidal", pc$set) | grepl("tidal", pc$observable)
      cl <- ifelse(c22, "#2a78d6", ifelse(tid, "#e69f00", "#b2182b"))
      ph <- ifelse(c22, 16, ifelse(tid, 15, 17))
      arrows(pc$age_ma, py - psd, pc$age_ma, py + psd, angle = 90, code = 3, length = 0.02, col = cl)
      points(pc$age_ma, py, pch = ph, col = cl)
    }
    if (!is.null(v)) {
      cc <- cl_est[[me]]
      arrows(age_ma, v[["r1"]], age_ma, v[["r2"]], angle = 90, code = 3, length = 0.03, lwd = 1.2, col = cc)
      segments(age_ma, v[["b1"]], age_ma, v[["b2"]], lwd = 4, col = cc)
      points(age_ma, v[["val"]], pch = 23, bg = cc, col = "black", cex = 1.5)
    }
    mtext(lab, 3, 1.2, adj = 0, font = 2, cex = 0.8)
  }
  kk <- if (!is.null(pc)) { 1296 / pc$k_arcsec } else { NULL }
  ksd <- if (!is.null(pc)) { kk * pc$k_sd / pc$k_arcsec } else { NULL }
  curve_panel(A$a_RE, A$a_min_RE, A$a_max_RE, pc$a_RE, pc$a_sd, "a_RE", "Earth-Moon distance (Earth radii)", "E  Earth-Moon distance")
  curve_panel(A$lod_h, A$lod_min_h, A$lod_max_h, pc$lod_h, pc$lod_sd, "lod_solar_h", "length of solar day (h)", "F  length of day")
  curve_panel(1296 / A$p_arcsec, 1296 / A$p_max_arcsec, 1296 / A$p_min_arcsec, kk, ksd, "k", "precession constant k (kyr)", "G  precession constant")
  # keys at fixed fractions of the strip width (not legend(), whose spacing is frozen in user units and clips or spreads
  # out when the plot window is resized)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  key <- function(y, tag, labs, lty, pch, col, bg = NA, lwd = 1) {
    x <- 0.06 + (seq_along(labs) - 1) * 0.94 / length(labs)
    text(0.005, y, tag, adj = 0, cex = 0.75, font = 2)
    for (q in seq_along(labs)) {
      if (!is.na(lty[q])) {
        segments(x[q], y, x[q] + 0.018, y, lty = lty[q], col = col[q], lwd = rep(lwd, length.out = length(labs))[q])
      }
      if (!is.na(pch[q])) {
        points(x[q] + 0.009, y, pch = pch[q], col = col[q], bg = rep(bg, length.out = length(labs))[q], cex = 1)
      }
      text(x[q] + 0.022, y, labs[q], adj = 0, cex = 0.7)
    }
  }
  key(0.78, "A", c("AR1 90%", "AR1 99%", sprintf("significant peak (AR1%s %g%%)", if (sp$ftest) { " + F-test" } else { "" }, 100 * sp$siglevel), "peaks of best solution", "peaks of other solutions",
      sprintf("AstroGeo22 period (k = %.2f)", k_ag), "Devonian composite (k = 20.70)"),
      lty = c(2, 3, NA, NA, NA, 2, 3), pch = c(NA, NA, 1, 25, 25, NA, NA), col = c("red3", "red3", "grey20", "black", "black", "grey30", "red3"), bg = c(NA, NA, NA, "grey60", "white", NA, NA))
  key(0.25, "B-G", c("obliquity", "precession", "obliquity + precession", "combined k: thick 1 sigma, thin 2 sigma", "AstroGeo22 (Farhat et al. 2022)", "Devonian composite (Zeeden et al. 2023)",
      "published: AstroGeo22 set / tidal / new"),
      lty = c(1, 1, 1, NA, 1, 3, NA), pch = c(16, 16, 16, NA, NA, NA, 17), col = c("grey25", "#56b4e9", "#0072b2", NA, "grey20", "red3", "#b2182b"), lwd = c(3, 3, 3, 1, 2.2, 1, 1))
  # pre-scan: record with removed intervals (grey) and running precession-band variance
  g <- sp$scan[[1]]
  x <- sp$raw[[1]]
  par(mar = c(3, 4, 2, 4))
  plot(x[, 1], x[, 2], type = "l", lwd = 0.5, col = "grey30", xlab = "time (kyr)", ylab = "value (std)", xaxs = "i")
  r <- rle(g$good)
  e <- cumsum(r$lengths)
  st <- e - r$lengths + 1
  if (any(!r$values)) {
    rect(g$time[st[!r$values]], -1e3, g$time[e[!r$values]], 1e3, col = adjustcolor("grey50", 0.35), border = NA)
  }
  if (!all(is.na(g$share))) {
    par(new = TRUE)
    plot(g$time, g$share, type = "l", col = "#d95f02", lwd = 1.5, axes = FALSE, xlab = "", ylab = "", xaxs = "i", log = "y")
    abline(h = sp$scan_thr * median(g$share), col = "#d95f02", lty = 2)
    axis(4, col.axis = "#d95f02")
    mtext("precession variance (log)", 4, 2.5, las = 0, cex = 0.6, col = "#d95f02")
  }
  mtext(sprintf("Pre-scan (%s): grey = removed, weak precession (%.0f%% kept); orange = running precession-band variance, dashed = threshold",
      names(sp$scan)[1], 100 * mean(g$good)), 3, 0.5, adj = 0, font = 2, cex = 0.75)
}

#' @title Check the precession-constant test on a synthetic record
#'
#' @description Builds synthetic records (five orbital terms, a 405 kyr cycle and AR(1) noise) with a
#' known precession constant and prints the k recovered by \code{\link{test_tuned}}. A check that the
#' chosen \code{tbw} resolves the terms.
#'
#' @param k_true Known precession constants of the synthetic records, in kyr (numeric vector); one record is built and
#' tested for each value. \code{Default=c(20.7, 21.0, 21.2, 21.4)}.
#' @param tbw MTM time-bandwidth product passed to \code{\link{test_tuned}} (numeric, > 0); see \code{tbw} there.
#' \code{Default=2}.
#' @param noise Standard deviation multiplier of the AR(1) noise (coefficient 0.7) added to the orbital signal
#' (numeric, >= 0); the orbital terms have amplitudes up to 1. \code{Default=1.5}.
#' @param seed Seed of the random number generator (integer), reset before each synthetic record so that the phases
#' and noise can be reproduced. \code{Default=1}.
#' @param len Length of the synthetic records in kyr, sampled every 1 kyr (numeric, > 0). \code{Default=5500}.
#' @param ... Further arguments passed on to \code{\link{test_tuned}} (e.g. \code{age_ma}, \code{k_search});
#' \code{tbw}, \code{n_surr = 2} and \code{plot = FALSE} are already set.
#' @return Prints one line per \code{k_true}.
#' @export
test_tuned_check <- function(k_true = c(20.7, 21.0, 21.2, 21.4), tbw = 2, noise = 1.5, seed = 1, len = 5500, ...) {
  fr <- c(-18.85, -26.348, 4.257452, 7.452, 17.916)
  amp <- c(1, 0.3, 1, 0.87, 0.7)
  t <- 0:len
  for (k in k_true) {
    set.seed(seed)
    P <- 1296 / (1296 / k + fr)
    cyc <- matrix(NA_real_, length(t), 5)
    for (i in 1:5) {
      cyc[, i] <- amp[i] * cos(2 * pi * t / P[i] + runif(1, 0, 2 * pi))
    }
    x <- rowSums(cyc) + 0.8 * cos(2 * pi * t / 405) +
        noise * as.numeric(arima.sim(list(ar = 0.7), length(t)))
    invisible(utils::capture.output(r <- test_tuned(data.frame(t, x), tbw = tbw, n_surr = 2, plot = FALSE, ...)))
    f <- function(x) if (is.na(x[1])) { "-" } else { sprintf("%.2f", x[1]) }
    cat(sprintf("true k %.2f -> obliquity %s, precession %s, obliquity + precession %s\n", k, f(r$k_obliquity), f(r$k_precession), f(r$k)))
  }
}
