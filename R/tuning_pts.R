#' @title Select tuning points in a plot
#'
#' @description select a point in a plot
#'Used in the \code{\link{astro_anchor}} function to select points
#'
#' @param x x coordinates of the points that can be selected (numeric vector), or any object accepted by
#' \code{\link[grDevices]{xy.coords}} (e.g. a two-column matrix) when \code{y} is \code{NULL}. Required, no default.
#' @param y y coordinates of the points (numeric vector of the same length as \code{x}); \code{NULL} takes them from
#' \code{x}. \code{Default=NULL}.
#' @param n Number of points to select (integer); currently not used, one point is selected per call.
#' \code{Default=1}.
#' @param pch Plotting symbol (integer code of \code{\link[graphics]{points}}) used to mark the selected point in red.
#' \code{Default=19}.
#' @param ... Further arguments passed to \code{\link[graphics]{identify}} (e.g. \code{tolerance}).
#'
#'@return
#'Returns a matrix with the select points their xy coordinates
#'
#' @importFrom grDevices xy.coords
#' @importFrom graphics identify
#' @importFrom graphics points
#' @keywords internal
#' @noRd

tuning_pts <- function(x,
    y = NULL,
    n = 1,
    pch = 19,
    ...) {
  oldop <- options(warn = -1)
  on.exit(options(oldop), add = TRUE)

  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  sel <- rep(FALSE, length(x))
  res <- integer(0)
  ans <- identify(x[!sel], y[!sel], n = 1, plot = F,
      ...)
  ans <- which(!sel)[ans]
  points(x[ans], y[ans], pch = pch, col = "red")
  sel[ans] <- TRUE
  res <- c(res, ans)
  options(oldop)
  return(res)
}
