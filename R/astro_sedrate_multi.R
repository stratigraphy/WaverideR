#' @title Combined astronomical test for several proxies and sections
#'
#' @description Runs \code{\link{astro_sedrate_test}} on every series (proxy) and fit and
#' combines the evidence: within a section all series share one sedimentation rate, so their
#' same-rate evidence is combined rate by rate; sections have their own rates and are combined
#' at the end.
#'
#' @details
#' \strong{Which tests enter.} For every series and fit the TimeOpt significance (S_opt) and the
#' tuned-band peak score (S_peaks) of \code{astro_sedrate_test} enter as separate components
#' (not their joint score, which loses power on short records); a component enters the combination
#' at a rate only where it is testable there (enough samples per carrier cycle and enough modulating cycles
#' after the edge trim). Tests that cannot work on a series (for example the 405-kyr fit on a record
#' that spans fewer than two 405-kyr cycles at that rate) therefore do not dilute the evidence of
#' the others. Weights can be given per series (for example to down-weight a proxy known to be
#' diagenetically overprinted); they must be chosen before looking at the results.
#'
#' \strong{Combination within a section.} \code{combine = "stouffer"} sums the
#' standard-normal scores of the same-rate p-values, z = sum(w z_i) / sqrt(w' R w), where R is
#' the correlation between the components' z-curves over the testable rates (clamped to 0-1;
#' proxies of one section share noise, and fit 1 and fit 2 of one proxy share data, so the
#' components are not independent). This is the most powerful combination when several proxies
#' carry the signal. \code{combine = "minp"} (default) takes the best component at each rate (most
#' powerful when only some proxies or tests carry the signal, the usual case; it is calibrated with
#' independent noise series per component and is therefore slightly conservative for correlated
#' proxies). Either statistic is calibrated with the noise series of the
#' components (same-rate p, rates within \code{null_smooth} steps pooled, exponential tail), and the
#' family-wise p over the rates of the section uses the minimum-p method.
#'
#' \strong{Combination of sections.} Section-level family-wise p-values are combined with
#' Fisher's method (sections are independent records): p_all.
#'
#' Independent age control (biostratigraphic duration bounds, magnetostratigraphy, radioisotopic
#' ages) is given through \code{sedmin} and \code{sedmax} per section; a narrower,
#' independently justified range reduces the number of rates tested and so raises the attainable
#' significance (Kemp, 2016).
#'
#' @param data The series to test (named list; each element an evenly spaced two-column matrix or data frame
#' with depth in metres in column 1 and proxy value in column 2). Series of the same section must share the
#' depth scale. Unnamed elements are named \code{"series1"}, \code{"series2"}, ... Required, no default
#' (\code{Default=NULL} gives an error).
#' @param section Section label of every series (character or numeric vector, same length as \code{data});
#' series with the same label share one sedimentation rate and are combined rate by rate. \code{NULL} puts all
#' series in one section. \code{Default=NULL}.
#' @param sedmin,sedmax Lowest and highest trial sedimentation rate in cm/kyr (numeric, > 0): one value for all
#' sections, or one value per section (named by section label, or in the order of \code{unique(section)}).
#' \code{Default=1} and \code{Default=20}.
#' @param numsed Number of trial rates per section (integer, > 0), log-spaced unless \code{linLog = 0} is passed
#' through \code{...}. \code{Default=60}.
#' @param fit TimeOpt fits run for every series (numeric vector with \code{1} (precession modulated by
#' eccentricity), \code{2} (short eccentricity modulated by the 405-kyr cycle) or both); every fit of every series
#' is a separate component. \code{Default=c(1, 2)}.
#' @param combine How the components of one section are combined at each rate (character): \code{"minp"} takes the
#' best component, \code{"stouffer"} sums their correlation-corrected standard-normal scores (see Details).
#' \code{Default=c("minp", "stouffer")}, i.e. \code{"minp"}.
#' @param scores Evidence used per series and fit, entered as separate components (character, one or both of
#' \code{"TimeOpt"} (the TimeOpt significance S_opt) and \code{"peaks"} (the tuned-band peak score S_peaks)).
#' \code{Default=c("TimeOpt", "peaks")}.
#' @param weights Weights of the series in the Stouffer combination (named numeric vector, names as in
#' \code{data}); series not named get weight 1. \code{NULL} gives all series equal weight. \code{Default=NULL}.
#' @param age_ma Age of the record in Ma (numeric) for the AstroGeo22 target periods; passed to
#' \code{\link{astro_sedrate_test}}. \code{NULL} uses present-day periods. \code{Default=NULL}.
#' @param targetE,targetP,targetO Your own eccentricity, precession and obliquity target periods in kyr (numeric
#' vectors), passed to \code{\link{astro_sedrate_test}}; \code{NULL} uses the La2004/AstroGeo22 periods.
#' \code{Default=NULL}.
#' @param n_null Number of noise series simulated per series and fit (integer, > 0); passed to
#' \code{\link{astro_sedrate_test}}. \code{Default=200}.
#' @param null_smooth Number of neighbouring rates on each side pooled for the same-rate p and smoothing of the
#' noise statistics (integer, >= 0); passed to \code{\link{astro_sedrate_test}} and used for the combined score.
#' \code{Default=2}.
#' @param run_multicore Simulate the noise series in parallel on several cores (logical, TRUE or FALSE); passed to
#' \code{\link{astro_sedrate_test}}. \code{Default=FALSE}.
#' @param n_cores Number of cores for the parallel computation (integer); \code{NULL} chooses automatically (see
#' \code{\link{astro_sedrate_test}}). \code{Default=NULL}.
#' @param verbose Print a progress message for every series and fit (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param ... Further arguments passed to \code{\link{astro_sedrate_test}}, e.g. \code{linLog}, \code{coco},
#' \code{tune} or \code{null_fit} (not \code{keep_null}, \code{sensitivity} or \code{target_sensitivity},
#' which are set by this function).
#' @return A list of class \code{astro_sedrate_multi} with
#' \describe{
#'   \item{sections}{per section: data frame of rate, combined score S_multi (-log10 same-rate p),
#'   same-rate p, family-wise p and the number of components used}
#'   \item{best}{per section: best rate, S_multi and family-wise p}
#'   \item{components}{per series and fit: best rate and family-wise p of its own score}
#'   \item{p_all}{Fisher combination of the section family-wise p-values}
#'   \item{runs}{the \code{astro_sedrate_test} results}
#' }
#' @references
#' Kemp, D.B. (2016). Optimizing significance testing of astronomical forcing in
#' cyclostratigraphy. Paleoceanography 31, 1516-1531 <doi:10.1002/2016PA002963>
#'
#' Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. and Williams, R.M. (1949). The
#' American Soldier, Vol. 1: Adjustment during Army Life. Princeton University Press.
#'
#' Westfall, P.H. and Young, S.S. (1993). Resampling-Based Multiple Testing: Examples and Methods
#' for p-Value Adjustment. Wiley, New York, ISBN 978-0-471-55761-6.
#' @examples
#' \donttest{
#' set.seed(1)
#' d1 <- mag[1:300, ]
#' d2 <- cbind(d1[, 1], d1[, 2] + rnorm(300, sd = sd(d1[, 2])))
#' res <- astro_sedrate_multi(list(a = d1, b = d2), sedmin = 0.4, sedmax = 1.5, numsed = 12,
#'     fit = 1, n_null = 10, coco = FALSE)
#' res$best
#' }
#' @export
astro_sedrate_multi <- function(data = NULL,
    section = NULL,
    sedmin = 1,
    sedmax = 20,
    numsed = 60,
    fit = c(1, 2),
    combine = c("minp", "stouffer"),
    scores = c("TimeOpt", "peaks"),
    weights = NULL,
    age_ma = NULL,
    targetE = NULL,
    targetP = NULL,
    targetO = NULL,
    n_null = 200,
    null_smooth = 2,
    run_multicore = FALSE,
    n_cores = NULL,
    verbose = FALSE,
    ...) {
  combine <- match.arg(combine)
  scores <- match.arg(scores, several.ok = TRUE)
  if (!is.list(data) || is.data.frame(data) || length(data) < 1) {
    stop("astro_sedrate_multi: data must be a list of series", call. = FALSE)
  }
  if (is.null(names(data)) || any(names(data) == "")) {
    names(data) <- paste0("series", seq_along(data))
  }
  if (is.null(section)) {
    section <- rep("1", length(data))
  }
  section <- as.character(section)
  if (length(section) != length(data)) {
    stop("astro_sedrate_multi: section must have one label per series", call. = FALSE)
  }
  secs <- unique(section)
  pick <- function(v, s) {
    if (length(v) == 1) {
      return(v)
    }
    if (!is.null(names(v))) {
      return(v[[s]])
    }
    v[match(s, secs)]
  }
  w_all <- rep(1, length(data))
  names(w_all) <- names(data)
  if (!is.null(weights)) {
    w_all[names(weights)] <- weights
  }
  runs <- list()
  comp <- NULL
  sec_out <- list()
  best <- NULL
  for (s in secs) {
    idx <- which(section == s)
    lo <- pick(sedmin, s)
    hi <- pick(sedmax, s)
    Zo <- NULL
    Zs <- list()
    wv <- numeric(0)
    rates <- NULL
    for (i in idx) {
      for (f in fit) {
        nm <- paste0(names(data)[i], "_fit", f)
        if (verbose) {
          message("astro_sedrate_multi: ", nm)
        }
        r <- tryCatch(astro_sedrate_test(data[[i]], sedmin = lo, sedmax = hi, numsed = numsed,
            fit = f, age_ma = age_ma, targetE = targetE, targetP = targetP, targetO = targetO,
            n_null = n_null, null_smooth = null_smooth, keep_null = TRUE, sensitivity = FALSE,
            target_sensitivity = FALSE, run_multicore = run_multicore, n_cores = n_cores, ...),
            error = function(e) { conditionMessage(e) })
        if (is.character(r)) {
          comp <- rbind(comp, data.frame(section = s, component = nm, best_sedrate = NA, S = NA,
              p_fw = NA, note = r))
          next
        }
        runs[[nm]] <- r
        sc <- r$scores
        rates <- sc$sedrate
        for (typ in scores) {
          Sv <- if (typ == "TimeOpt") { sc$S_opt } else { sc$S_peaks }
          Pn <- if (typ == "TimeOpt") { r$null_p$S_opt } else { r$null_p$S_peaks }
          pf <- if (typ == "TimeOpt") { sc$p_S } else { sc$p_peaks }
          if (is.null(Pn) || !any(is.finite(Sv))) {
            next
          }
          zo <- stats::qnorm(pmin(pmax(10^-Sv, 1e-15), 1 - 1e-15), lower.tail = FALSE)
          zs <- stats::qnorm(pmin(pmax(Pn, 1e-15), 1 - 1e-15), lower.tail = FALSE)
          Zo <- cbind(Zo, zo)
          Zs[[length(Zs) + 1]] <- zs
          wv <- c(wv, w_all[i])
          jb <- which.max(Sv)
          comp <- rbind(comp, data.frame(section = s, component = paste(nm, typ),
              best_sedrate = rates[jb], S = Sv[jb], p_fw = pf[jb], note = ""))
        }
      }
    }
    if (is.null(Zo)) {
      next
    }
    nc <- ncol(Zo)
    nr <- nrow(Zo)
    nsim <- ncol(Zs[[1]])
    # correlation between the components' z-curves (shared data and noise)
    R <- diag(nc)
    if (nc > 1) {
      for (a in 1:(nc - 1)) {
        for (b in (a + 1):nc) {
          ok <- is.finite(Zo[, a]) & is.finite(Zo[, b])
          rr <- if (sum(ok) > 5) { suppressWarnings(stats::cor(Zo[ok, a], Zo[ok, b])) } else { 0 }
          if (!is.finite(rr)) {
            rr <- 0
          }
          R[a, b] <- min(1, max(0, rr))
          R[b, a] <- R[a, b]
        }
      }
    }
    comb_one <- function(z, null_sim) {
      ok <- is.finite(z)
      if (!any(ok)) {
        return(NA_real_)
      }
      if (combine == "minp") {
        return(max(z[ok]))
      }
      w <- wv[ok]
      den <- if (null_sim) { sqrt(sum(w^2)) } else { sqrt(as.numeric(t(w) %*% R[ok, ok, drop = FALSE] %*% w)) }
      sum(w * z[ok]) / den
    }
    co <- rep(NA_real_, nr)
    cs <- matrix(NA_real_, nr, nsim)
    ncomp <- integer(nr)
    for (q in seq_len(nr)) {
      co[q] <- comb_one(Zo[q, ], FALSE)
      ncomp[q] <- sum(is.finite(Zo[q, ]))
      for (k in seq_len(nsim)) {
        zk <- numeric(nc)
        for (cc in seq_len(nc)) {
          zk[cc] <- Zs[[cc]][q, k]
        }
        zk[!is.finite(Zo[q, ])] <- NA
        cs[q, k] <- comb_one(zk, TRUE)
      }
    }
    okr <- is.finite(co)
    for (q in which(okr)) {
      if (sum(is.finite(cs[q, ])) < nsim / 2) {
        okr[q] <- FALSE
      }
    }
    sg <- wr_ast_sigspec(co, cs, okr, null_smooth)
    Sm <- -log10(sg$p)
    Sm[!okr] <- NA
    tab <- data.frame(sedrate = rates, S_multi = Sm, p_same_rate = sg$p, p_fw = sg$fw,
        n_components = ncomp)
    sec_out[[s]] <- tab
    if (any(is.finite(Sm))) {
      jb <- which.max(Sm)
      best <- rbind(best, data.frame(section = s, best_sedrate = rates[jb], S_multi = Sm[jb],
          p_fw = sg$fw[jb], n_series = length(idx)))
    }
  }
  p_all <- NA_real_
  if (!is.null(best) && any(is.finite(best$p_fw))) {
    pp <- best$p_fw[is.finite(best$p_fw)]
    p_all <- stats::pchisq(-2 * sum(log(pp)), 2 * length(pp), lower.tail = FALSE)
  }
  out <- list(sections = sec_out, best = best, components = comp, p_all = p_all, runs = runs,
      settings = list(combine = combine, fit = fit, n_null = n_null, numsed = numsed))
  class(out) <- "astro_sedrate_multi"
  out
}

#' @export
print.astro_sedrate_multi <- function(x, ...) {
  cat("astro_sedrate_multi:", nrow(x$components), "components,", length(x$sections),
      "section(s), combination", x$settings$combine, "\n\nper section:\n")
  print(x$best, digits = 3, row.names = FALSE)
  cat("\ncomponents (own best rate and family-wise p):\n")
  print(x$components[, c("section", "component", "best_sedrate", "S", "p_fw")], digits = 3,
      row.names = FALSE)
  if (length(x$sections) > 1) {
    cat("\nall sections combined (Fisher): p =", signif(x$p_all, 3), "\n")
  }
  invisible(x)
}

#' @title Plot an astro_sedrate_multi result
#' @description Combined significance spectrum (-log10 same-rate p) per section with its
#' family-wise 5 \% level, and the significance spectra of the single components (thin lines).
#' @param x Result to plot: an object of class \code{astro_sedrate_multi}, as returned by
#' \code{\link{astro_sedrate_multi}}. Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{plot()}.
#' @return A plot.
#' @export
plot.astro_sedrate_multi <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  ns_ <- length(x$sections)
  graphics::par(mfrow = c(ns_, 1), mar = c(4.2, 4.5, 2.2, 1))
  for (s in names(x$sections)) {
    tab <- x$sections[[s]]
    comps <- x$components$component[x$components$section == s]
    allv <- c(3, tab$S_multi[is.finite(tab$S_multi)])
    for (cn in comps) {
      if (!is.null(x$runs[[cn]])) {
        sc <- x$runs[[cn]]$scores
        v <- ifelse(is.finite(sc$S_joint), sc$S_joint, sc$S_opt)
        allv <- c(allv, v[is.finite(v)])
      }
    }
    plot(tab$sedrate, tab$S_multi, type = "n", log = "x", ylim = c(0, max(allv) * 1.08), las = 1,
        xlab = "sedimentation rate (cm/kyr)", ylab = expression(-log[10] ~ p ~ "(same rate)"),
        main = paste("section", s))
    for (cn in comps) {
      if (!is.null(x$runs[[cn]])) {
        sc <- x$runs[[cn]]$scores
        graphics::lines(sc$sedrate, ifelse(is.finite(sc$S_joint), sc$S_joint, sc$S_opt), col = "grey70")
      }
    }
    graphics::lines(tab$sedrate, tab$S_multi, lwd = 2.4, col = "#0b0b0b")
    b <- x$best[x$best$section == s, ]
    if (nrow(b)) {
      graphics::points(b$best_sedrate, b$S_multi, pch = 19)
      graphics::mtext(sprintf("best %.3g cm/kyr, family-wise p %.3g", b$best_sedrate, b$p_fw),
          side = 3, line = 0.3, cex = 0.8)
    }
    graphics::legend("topleft", c("combined", "single components"), col = c("#0b0b0b", "grey70"),
        lwd = c(2.4, 1), bty = "n", cex = 0.8)
  }
  invisible(NULL)
}
