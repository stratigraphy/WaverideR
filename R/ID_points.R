#' @title Interactively select points in a plot
#'
#' @description interactively select points
#'Used in the \code{\link{track_period_wavelet}} function to track the ridge/period of a certain cycle
#'
#'
#'@param x x coordinates of the plotted points among which points are selected (numeric vector), or any object
#'accepted by \code{\link[grDevices]{xy.coords}} (e.g. a two-column matrix of x and y). Required, no default.
#'@param y y coordinates of the plotted points (numeric vector of the same length as \code{x});
#'\code{NULL} when \code{x} already holds both coordinates. \code{Default=NULL}.
#'@param n Maximum number of points that can be selected (integer, > 0); selection stops earlier when the user
#'ends the \code{\link[graphics]{identify}} session. \code{Default=length(x)}.
#'@param pch Plotting symbol (integer \code{pch} code, see \code{\link[graphics]{points}}) used to mark each
#'selected point in white. \code{Default=19}.
#'@param ... Not used; accepted for compatibility only.
#'
#' @importFrom grDevices xy.coords
#' @importFrom graphics identify
#' @importFrom graphics points
#' @keywords internal
#' @noRd


ID_points <- function(x,
    y = NULL,
    n = length(x),
    pch = 19,
    ...) {
  oldop <- options(warn = -1)
  on.exit(options(oldop), add = TRUE)
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  sel <- rep(FALSE, length(x))
  res <- integer(0)
  while (sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = F)
    if (!length(ans)) {
      break
    }
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch, col = "white")
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  options(oldop)
  return(res)
}


