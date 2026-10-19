#' @title Plot the result of bayes_changepoint
#'
#' @description Plots the data, the posterior-mean regression model with its
#' pointwise credible interval, and the posterior probability of a change
#' point (right axis), as in the Acycle change point figure (bayes_cpt.m),
#' together with the posterior distribution of the number of change points.
#'
#' @param bcp Output of \code{\link{bayes_changepoint}}.
#' @param x_lab Label of the x-axis \code{Default="depth"}
#' @param y_lab Label of the proxy axis \code{Default="proxy"}
#' @param plot_k Add a panel with the posterior distribution of the number
#' of change points \code{Default=TRUE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Eric Ruggieri,
#' Mingsong Li.
#'
#' @references
#' Ruggieri, E. (2013). A Bayesian approach to detecting change points in
#' climatic records. International Journal of Climatology, 33, 520-528.
#' <doi:10.1002/joc.3447>
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' y <- c(rep(0, 60), rep(1.5, 60)) + stats::rnorm(120, sd = 0.5)
#' bcp <- bayes_changepoint(data = cbind(1:120, y), k_max = 4, d_min = 5,
#'                          v_0 = 3, n_sim = 200, seed = 1)
#' plot_bayes_changepoint(bcp)
#' }
#'
#' @return No return value, called for its side effect (a base graphics plot).
#'
#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics polygon
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics barplot
#' @importFrom graphics layout
#' @importFrom graphics legend
#' @importFrom grDevices adjustcolor
plot_bayes_changepoint <- function(bcp = NULL,
                                   x_lab = "depth",
                                   y_lab = "proxy",
                                   plot_k = TRUE){
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  if (plot_k == TRUE){
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  }
  graphics::par(mar = c(4, 4, 2, 4))
  x <- bcp$data[, 1]
  y <- bcp$data[, 2]
  ci <- bcp$model_ci
  yr <- range(c(y, ci[, 2], ci[, 4]), finite = TRUE)
  graphics::plot(x, y, type = "l", lwd = 2, col = "steelblue", ylim = yr,
                 xlab = x_lab, ylab = y_lab, main = "Change points")
  ok <- is.finite(ci[, 2]) & is.finite(ci[, 4])
  if (any(ok)){
    graphics::polygon(c(x[ok], rev(x[ok])), c(ci[ok, 2], rev(ci[ok, 4])),
                      col = grDevices::adjustcolor("green3", alpha.f = 0.3),
                      border = NA)
  }
  graphics::lines(x, bcp$model[, 2], col = "darkgreen", lty = 2, lwd = 2)
  graphics::par(new = TRUE)
  graphics::plot(x, bcp$cpt[, 2], type = "l", col = "red", lwd = 2,
                 ylim = c(0, 2), axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, at = seq(0, 1, by = 0.2), col = "red", col.axis = "red")
  graphics::mtext("posterior probability", side = 4, line = 2.5, col = "red")
  if (plot_k == TRUE){
    graphics::par(mar = c(4, 5.5, 2, 1))
    kp <- bcp$k_posterior
    graphics::barplot(kp[, 2], names.arg = kp[, 1], col = "grey60",
                      xlab = "number of change points",
                      ylab = "posterior probability")
  }
  invisible(NULL)
}
