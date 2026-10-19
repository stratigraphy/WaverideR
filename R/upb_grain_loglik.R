#' Eruption-age likelihoods from single-zircon U-Pb dates (for astroMCBchron)
#'
#' @description
#' \code{upb_grain_loglik()} turns single-zircon dates of one or more ash beds into one
#' eruption-age log-likelihood function per bed, for the \code{date_loglik} argument of
#' \code{\link{astroMCBchron}}. Instead of a weighted mean, every grain is modelled as a
#' draw from a mixture of
#' \itemize{
#'   \item a crystallisation interval: the grain crystallised uniformly between the eruption
#'   age \emph{te} and \emph{te} + \emph{D} (antecrystic growth before eruption), with the
#'   duration \emph{D} marginalised over a uniform prior between 0 and \code{D_max};
#'   \item inherited or recycled grains (probability \code{p_inherit}): uniform over
#'   \code{W_inherit} before the eruption;
#'   \item grains with residual Pb loss (probability \code{p_pbloss}): uniform over
#'   \code{W_pbloss} after the eruption,
#' }
#' each convolved with the analytical uncertainty of the grain. The likelihood of \emph{te}
#' is computed on a grid, so no grain has to be rejected beforehand. With \code{D_max} = 0
#' and no mixture components the result is the likelihood of a weighted mean; with a large
#' \code{D_max} the eruption age approaches the youngest concordant grains (as in the
#' Bayesian eruption-age estimator of Keller et al., 2018). Which regime applies depends on
#' the ratio of the crystallisation interval to the analytical uncertainty of the grains,
#' which is why the choice is left to the data rather than fixed beforehand.
#'
#' The functions return a log-likelihood (not a posterior): no prior is placed on \emph{te},
#' so the stratigraphic order and the astrochronology in \code{astroMCBchron} act as the
#' prior. \code{astroMCBchron} normalises each function numerically and combines it with
#' its outlier mixture.
#'
#' @param grains data.frame with columns \code{id} (bed), \code{age} (single-grain date,
#'   in \code{age_unit}) and either \code{sd} (1 sigma) or \code{sd2} (2 sigma, analytical).
#' @param D_max maximum crystallisation interval (age units). Default 0.6 (Ma).
#' @param p_inherit,p_pbloss probabilities of inherited and Pb-loss grains. Defaults 0.1 and 0.05.
#' @param W_inherit,W_pbloss widths of the inherited and Pb-loss components (age units).
#'   Defaults 20 and 3 (Ma).
#' @param step grid step of the eruption age (age units). Default 0.002 (2 kyr for Ma).
#' @param nD number of grid values for the marginalisation over \emph{D}. Default 61.
#' @return A list with \code{loglik} (named list of functions, one per bed, of the eruption
#'   age in age units) and \code{summary} (data.frame per bed: number of grains, maximum-
#'   likelihood eruption age, mean and 2 sd of the normalised likelihood, youngest grain and
#'   weighted mean with MSWD for comparison).
#' @references
#' Keller, C.B., Schoene, B., Samperton, K.M., 2018. A stochastic sampling approach to zircon
#' eruption age interpretation. Geochemical Perspectives Letters 8, 31-35.
#' @examples
#' set.seed(2)
#' g <- data.frame(id = rep(c("ashA", "ashB"), each = 8),
#'                 age = c(430 + runif(8, 0, 0.1) + rnorm(8, 0, 0.04), 431 + runif(8, 0, 0.1) + rnorm(8, 0, 0.04)),
#'                 sd = 0.04)
#' g$age[3] <- 433.2                                     # one inherited grain
#' gl <- upb_grain_loglik(g)
#' gl$summary
#' @export
upb_grain_loglik <- function(grains, D_max = 0.6, p_inherit = 0.1, p_pbloss = 0.05, W_inherit = 20, W_pbloss = 3,
                             step = 0.002, nD = 61) {
  grains <- as.data.frame(grains)
  if (is.null(grains$id) || is.null(grains$age)) stop("grains needs columns 'id' and 'age'")
  s_all <- if (!is.null(grains$sd)) grains$sd else if (!is.null(grains$sd2)) grains$sd2 / 2 else stop("grains needs 'sd' (1 sigma) or 'sd2' (2 sigma)")
  if (p_inherit + p_pbloss >= 1) stop("p_inherit + p_pbloss must be < 1")
  ids <- unique(grains$id); out <- list(); summ <- list()
  for (b in ids) { k <- grains$id == b & is.finite(grains$age) & is.finite(s_all); x <- grains$age[k]; s <- s_all[k]
    if (length(x) < 1) next
    tg <- seq(min(x) - max(1.2, 6 * max(s)), max(x) + max(0.3, 3 * max(s)), by = step); Dg <- seq(0, D_max, length.out = nD)
    comp <- function(te, D) if (D < 1e-9) stats::dnorm(x, te, s) else (stats::pnorm((te + D - x) / s) - stats::pnorm((te - x) / s)) / D
    bg <- function(te) p_inherit * (stats::pnorm((te + W_inherit - x) / s) - stats::pnorm((te - x) / s)) / W_inherit +
      p_pbloss * (stats::pnorm((te - x) / s) - stats::pnorm((te - W_pbloss - x) / s)) / W_pbloss
    ll <- vapply(tg, function(te) { h <- bg(te); v <- vapply(Dg, function(D) sum(log((1 - p_inherit - p_pbloss) * comp(te, D) + h + 1e-300)), 0)
      m <- max(v); m + log(mean(exp(v - m))) }, 0)
    ll <- ll - max(ll); w <- exp(ll); w <- w / sum(w); mu <- sum(w * tg); sdv <- sqrt(sum(w * (tg - mu)^2))
    f <- local({ tg <- tg; ll <- ll; fi <- stats::approxfun(tg, ll, rule = 1)
      function(a) { v <- fi(a); out <- ifelse(is.finite(v), v, NA_real_); bad <- !is.finite(out)
        if (any(bad)) out[bad] <- min(ll) - 50 - 1e3 * pmin(abs(a[bad] - min(tg)), abs(a[bad] - max(tg)))
        out } })
    ww <- 1 / s^2; wm <- sum(ww * x) / sum(ww); msw <- if (length(x) > 1) sum(ww * (x - wm)^2) / (length(x) - 1) else NA
    out[[as.character(b)]] <- f
    summ[[length(summ) + 1]] <- data.frame(id = b, n = length(x), te_ml = tg[which.max(ll)], te_mean = mu, te_2sd = 2 * sdv,
                                           youngest = min(x), youngest_2sd = 2 * s[which.min(x)], wm_all = wm, wm_all_2sd = 2 / sqrt(sum(ww)), mswd_all = msw) }
  list(loglik = out, summary = do.call(rbind, summ))
}
