#' @title Discriticizes lithologs
#'
#' @description Discriticizes lithologs to allow further time-series analysis first the
#' Greatest common divisor/highest common factor is calculated which is then used to discriticize
#' the litholog to an evenly sampled data series. The function is designed to place the boundary
#' at the original depth level of the bed boundaries. The Greatest common divisor/highest common factor can
#' be a very small number as such the discriticized data set can be large which impacts computational
#' performance later on therefore a linear interpolation option is added to downscale the data to allow
#' for computational efficiency later on. This is made to discriticize lithologs created using the
#' 'StratigrapheR' package. as such the same data format for input is used.
#' eg. column 1 is bottom of the bed, column 2 is top of bed, column is depth rank/proxy value
#'
#'@references
#' Wouters, S., Da Silva, A.-C., Boulvain, F., and Devleeschouwer, X.. 2021.
#' StratigrapheR: Concepts for Litholog Generation in R.
#' The R Journal. \doi{<doi:10.32614/RJ-2021-039>}
#'
#'
#' @param litholog litholog input matrix with 3 columns column 1 is bottom of the bed,
#'  column 2 is top of bed, column is depth rank/proxy value
#' @param subset_fact subset factor which is x times the greatest common divider \code{Default=100}.
#' @param lin_interp Linear interpolation of the data set \code{Default=FALSE}
#' @param dt step size  \code{Default=NULL}.
#' @param genplot generate plot  \code{Default=FALSE}
#' @param x_lab label for the y-axis \code{Default="rank"}
#'@param y_lab label for the y-axis \code{Default="depth (m)"}
#'@param keep_editable Keep option to add extra features after plotting  \code{Default=FALSE}
#'
#' @examples
#'# Convert depth rank record to a discrete proxy record to allow for further
#'# analysis in which discrete time series are needed
#'depth_rank_example_disc <- lithlog_disc(litholog = depth_rank_example,
#'            subset_fact = 100,
#'            genplot = FALSE,
#'            x_lab = "rank",
#'            y_lab = "depth (m)",
#'            keep_editable=FALSE)
#'
#'
#'@return
#'Returns a matrix with 2 columns, the first column is depth the second columns is the depth/rank proxy
#'If genplot is \code{Default=TRUE} then a plot of the discriticizes time series is plotted
#' @export


lithlog_disc <- function(litholog = NULL,
                         subset_fact = 100,
                         lin_interp = FALSE,
                         dt = NULL,
                         genplot = FALSE,
                         x_lab = "rank",
                         y_lab = "depth (m)",
                         keep_editable = FALSE) {
  a <- -1
  m <- 0.999

  while (round(m, 0) != round(m, 3)) {
    a = a + 1
    x <- c(litholog[, 2] - litholog[, 1]) * (10 ^ a)
    m = min(x)
    while (any(round(x %% m, 4) > 0)) {
      m = m - 1
    }
  }

  m <- m / (10 ^ a)

  subsetsample_nr <- m / subset_fact
  bed.example <- litholog
  col_interp <- matrix(data = NA,
                       nrow = 0,
                       ncol = 2)
  colnames(col_interp) <- c("new_nrs", "val")
  col_interp <- na.omit(col_interp)

  for (i in 1:nrow(bed.example)) {
    if (i == 1) {
      bottom <- bed.example[i, 1]
      top <- bed.example[i, 2]
      new_nrs <-
        seq(bottom, top, subsetsample_nr)
      val <-
        c(rep(bed.example[i, 3], times = length(new_nrs) - 1), ((bed.example[i, 3] +
                                                                   bed.example[i + 1, 3]) / 2))
      col_interp <- rbind(col_interp, cbind(new_nrs, val))
    }

    if (i > 1 & i != nrow(bed.example)) {
      bottom <- bed.example[i, 1]
      top <- bed.example[i, 2]

      new_nrs <-
        seq(bottom + subsetsample_nr,
            top,
            subsetsample_nr)
      val <-
        c(rep(bed.example[i, 3], times = length(new_nrs) - 1), ((bed.example[i, 3] +
                                                                   bed.example[i + 1, 3]) / 2))
      col_interp <- rbind(col_interp, cbind(new_nrs, val))
    }

    if (i == nrow(bed.example)) {
      bottom <- bed.example[i, 1]
      top <- bed.example[i, 2]
      new_nrs <- seq(bottom + subsetsample_nr,
                     top, subsetsample_nr)
      val <- c(rep(bed.example[i, 3], times = length(new_nrs)))
      col_interp <- rbind(col_interp, cbind(new_nrs, val))
    }
  }


  if (lin_interp == TRUE) {
    col_interp <- na.omit(col_interp)
    yleft_out <- col_interp[1, 2]
    yright_out <- col_interp[nrow(col_interp), 2]
    seq <-
      seq(
        from = min(col_interp[, 1]),
        to = max(col_interp[, 1]),
        by = dt
      )
    app <- approx(
      x = col_interp[, 1],
      y = col_interp[, 2],
      xout = seq,
      method = "linear",
      yleft = yleft_out,
      yright = yright_out
    )
    col_interp <- as.data.frame(cbind(app$x, app$y))
  }



  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    plot(
      col_interp[, 2],
      col_interp[, 1],
      type = "l",
      col = "black",
      xlab = x_lab,
      ylab = y_lab
    )
  }
return(col_interp)
}
