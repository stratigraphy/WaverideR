#' @title Interactively select points in a plot
#'
#' @description interactively select points
#'Used in the \code{\link{track_period_wavelet}} function to track the ridge/period of a certain cycle
#'
#'
#'@param x x axis
#'@param y y axis
#'@param n number of points
#'@param pch size of point
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
  defaultW <- getOption("warn")
  options(warn = -1)
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  sel <- rep(FALSE, length(x))
  res <- integer(0)
  while (sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n = 1, plot = F)
    if (!length(ans))
      break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], pch = pch, col = "white")
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  options(warn=defaultW)
  return(res)
}



