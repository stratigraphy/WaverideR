#' @title Select tuning points in a plot
#'
#' @description select a point in a plot
#'Used in the \code{\link{astro_anchor}} function to select points
#'
#'@param x x axis
#'@param y y axis
#'@param n number of points
#'@param pch size of point
#'
#'@return
#'Returns a matrix with the select points their xy coordinates
#'
#' @importFrom grDevices xy.coords
#' @importFrom graphics identify
#' @importFrom graphics points
#' @keywords internal
#' @noRd

tuning_pts <- function(x, y = NULL, n = 1, pch = 19,
                       ...) {
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
  res
}
