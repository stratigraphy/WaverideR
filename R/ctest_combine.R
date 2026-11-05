#' @title Combine several CTest results by multiplying their densities
#'
#' @description Reproduces the combination step of mmc8
#' (reconstruct_parameters_Devonian.R), using a base-R weighted quantile
#' instead of Hmisc::wtd.quantile.
#'
#' @param results CTest results to combine (list, ideally named): each element is a \code{ctest} object from
#' \code{\link{ctest}} or a numeric vector of samples (e.g. precession constants in years). The kernel densities
#' of all elements are multiplied. Required, no default.
#' @param what Column of the \code{samples} data frame of the \code{ctest} objects that is combined (character),
#' e.g. \code{"k"} (precession constant, yr), \code{"SR"} (sedimentation rate, m/kyr), \code{"obl_ks3"} or
#' \code{"g2"}; ignored for numeric vectors. \code{Default="k"}.
#' @param n Number of points of the grid on which the densities are evaluated and multiplied (integer, > 1).
#' \code{Default=1000}.
#' @return list: grid `x`, combined density, mean (median of the combined
#'   density) and sd (half the inter-quartile span, as in the paper)
#' @export
ctest_combine <- function(results, what = "k", n = 1000){
  nres <- length(results)
  vals <- vector("list", nres)
  names(vals) <- names(results)
  mins <- numeric(nres)
  maxs <- numeric(nres)
  for (j in seq_len(nres)){
    r <- results[[j]]
    if (is.numeric(r)){
      vals[j] <- list(r)
    } else {
      vals[j] <- list(r$samples[[what]])
    }
    mins[j] <- min(vals[[j]])
    maxs[j] <- max(vals[[j]])
  }
  MM <- c(min(mins), max(maxs))
  x <- seq(MM[1], MM[2], length.out = n)
  dens <- NULL
  for (j in seq_len(nres)){
    dj <- stats::density(vals[[j]], from = MM[1], to = MM[2], n = n)$y
    if (j == 1){
      dens <- dj
    } else {
      dens <- dens * dj
    }
  }
  dens <- dens / sum(dens)
  wq <- function(prob){
    cw <- cumsum(dens) / sum(dens)
    stats::approx(cw, x, xout = prob, ties = "ordered", rule = 2)$y
  }
  q <- wq(c(0.25, 0.5, 0.75))
  list(x = x, density = dens, median = q[2],
      sd = ((q[2] - q[1]) + (q[3] - q[2])) / 2,
      shared_range = c(max(mins), min(maxs)),
      quantiles = stats::setNames(q, c("25%", "50%", "75%")))
}
