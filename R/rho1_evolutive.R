#' @title Evolutionary lag-1 autocorrelation coefficient (rho1)
#'
#' @description The \code{\link{rho1_evolutive}} function calculates the
#' lag-1 autocorrelation coefficient of an AR(1) model in moving windows, as
#' used in the sedimentary noise model of Li et al. (2018a,b).
#' With \code{n_sim = 0} a single deterministic run is made
#' (\code{acycleEvolutionaryRho1.m}): the series is linearly interpolated on
#' a grid anchored at the first coordinate with interval
#' \code{sample_min}, windows of \code{round(window_min/sample_min)} points
#' are linearly detrended and rho1 = sum(z[t] z[t-1]) / sum(z[t-1]^2) of the
#' demeaned residuals z is reported at the midpoint of each window.
#' With \code{n_sim > 0} the Acycle Monte Carlo is run (AC.m menu
#' "Lag-1 autocorrelation coefficient (rho1)"): for every simulation a
#' sampling interval and a window are drawn uniformly from
#' [\code{sample_min}, \code{sample_max}] and [\code{window_min},
#' \code{window_max}], the data are linearly interpolated
#' (\code{interpolate.m}), evolutionary rho1 is computed with
#' \code{erhoAR1.m} (window of round(window/dt) points, linear detrending,
#' \code{rhoAR1ML.m}), window centers are optionally shifted by a random
#' offset uniform in [-window/2, window/2], and the curves are linearly
#' interpolated on a regular grid and summarized by the percentiles 2.5, 5,
#' 10, 15.865, 25, 50, 75, 84.135, 90, 95 and 97.5.
#' Translated from Acycle v3.1 acycleEvolutionaryRho1.m,
#' acycleInterpolateSeries.m, erhoAR1.m, rhoAR1ML.m, interpolate.m and AC.m
#' (menu_rho_Callback) (Mingsong Li; rhoAR1 by Dorothee Husson).
#'
#' @param data Matrix with two columns: depth/time and proxy value (sorted,
#' non-finite rows removed). \code{Default=NULL}
#' @param window_min Window length (single run) or minimum window length
#' (Monte Carlo). If NULL: single run min(400, 0.3 * data length); Monte
#' Carlo 400 if 0.3 * data length > 400, otherwise 0.3 * data length.
#' \code{Default=NULL}
#' @param window_max Maximum window length (Monte Carlo). If NULL: 500 if
#' 0.3 * data length > 400, otherwise 0.4 * data length. \code{Default=NULL}
#' @param sample_min Sampling interval (single run) or minimum sampling
#' interval (Monte Carlo). If NULL the 95th percentile of the data sampling
#' intervals. \code{Default=NULL}
#' @param sample_max Maximum sampling interval (Monte Carlo). If NULL twice
#' the 95th percentile of the data sampling intervals. \code{Default=NULL}
#' @param n_sim Number of Monte Carlo simulations (at least 50, as in
#' Acycle); 0 gives a single deterministic run. \code{Default=1000}
#' @param n_out Number of grid points of the Monte Carlo output. If NULL
#' min(1000, number of data points). \code{Default=NULL}
#' @param shift_grid Randomly shift the window centers of every simulation
#' (Acycle "shift grids" > 1). \code{Default=TRUE}
#' @param seed Seed for the random number generator (NULL = not set).
#' \code{Default=NULL}
#' @param verbose Show a progress bar. \code{Default=FALSE}
#' @param genplot Plot the result. \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle), rhoAR1 by Dorothee Husson.
#'
#' @references
#' Li, M., Hinnov, L.A., Huang, C., Ogg, J.G., 2018a. Sedimentary noise and
#' sea levels linked to land-ocean water exchange and obliquity forcing.
#' Nature Communications 9, 1004. \doi{10.1038/s41467-018-03454-y}
#'
#' Li, M., Huang, C., Hinnov, L., Chen, W., Ogg, J., Tian, W., 2018b.
#' Astrochronology of the Anisian stage (Middle Triassic) at the Guandao
#' reference section, South China. Earth and Planetary Science Letters 482,
#' 591-606. \doi{10.1016/j.epsl.2017.11.042}
#'
#' Mudelsee, M., 2002. TAUEST: a computer program for estimating persistence
#' in unevenly spaced weather/climate time series. Computers & Geosciences
#' 28, 69-72.
#'
#' @examples
#' \donttest{
#' set.seed(3)
#' t <- cumsum(runif(800, 0.5, 1.5))
#' x <- as.numeric(arima.sim(list(ar = 0.8), n = 800))
#' single <- rho1_evolutive(data = cbind(t, x), n_sim = 0)
#' mc <- rho1_evolutive(data = cbind(t, x), n_sim = 50, seed = 1,
#'   genplot = TRUE)
#' }
#' @return For \code{n_sim = 0} a list with \code{rho1} (matrix [window
#' center, rho1]), \code{window_points}, \code{sample_interval} and
#' \code{window}. For the Monte Carlo a list with
#' \code{percentile} (matrix [x, 11 percentiles] on the grid points where
#' they are defined), \code{median} (matrix [x, median]),
#' \code{percentiles}, \code{grid}, \code{rho1} (n_out x n_sim matrix of all
#' realisations) and \code{draws} (sampling intervals, windows, shifts).
#'
#' @export
#' @importFrom stats runif
#' @importFrom stats quantile
#' @importFrom stats approx
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics polygon
#' @importFrom graphics lines
#' @importFrom grDevices rgb

rho1_evolutive <- function(data = NULL,
                           window_min = NULL,
                           window_max = NULL,
                           sample_min = NULL,
                           sample_max = NULL,
                           n_sim = 1000,
                           n_out = NULL,
                           shift_grid = TRUE,
                           seed = NULL,
                           verbose = FALSE,
                           genplot = FALSE){
  # ---- internal helper functions (local to this function) ----

  acy_detrend_linear <- function(x, t = NULL){
    x <- as.numeric(x)
    n <- length(x)
    if (n < 2){
      return(x - mean(x))
    }
    if (is.null(t)){
      t <- seq_len(n)
    }
    X <- cbind(1, t)
    beta <- qr.solve(X, x)
    x - as.numeric(X %*% beta)
  }

  acy_check_two_column <- function(data, name = "data"){
    data <- as.matrix(data)
    if (ncol(data) < 2){
      stop(paste0(name, " must have at least two columns: [depth/time, value]"))
    }
    data <- data[, 1:2, drop = FALSE]
    keep <- is.finite(data[, 1]) & is.finite(data[, 2])
    data <- data[keep, , drop = FALSE]
    data <- data[order(data[, 1]), , drop = FALSE]
    data
  }

  # @keywords internal
  # @noRd
  dynot_prctile <- function(x, p){
    # MATLAB prctile (midpoint rule, NaN ignored) == quantile type 5
    x <- x[!is.na(x)]
    if (length(x) == 0){
      return(rep(NA_real_, length(p)))
    }
    as.numeric(stats::quantile(x, probs = p / 100, type = 5, names = FALSE))
  }

  # @keywords internal
  # @noRd
  dynot_row_prctile <- function(mat, p){
    out <- matrix(NA_real_, nrow(mat), length(p))
    for (r in seq_len(nrow(mat))){
      out[r, ] <- dynot_prctile(mat[r, ], p)
    }
    out
  }

  # @keywords internal
  # @noRd
  dynot_interp1 <- function(x, y, xi){
    # MATLAB interp1 linear, NaN outside the range
    if (length(x) < 2){
      stop("At least two points are needed for interpolation")
    }
    stats::approx(x, y, xout = xi, method = "linear", rule = 1,
                  ties = "ordered")$y
  }

  # @keywords internal
  # @noRd
  dynot_erho <- function(data, window){
    # erhoAR1.m with rhoAR1ML.m
    dt <- data[2, 1] - data[1, 1]
    n <- nrow(data)
    npts <- pda_round(window / dt)
    m <- n - npts + 1
    if (m < 1 || npts < 2){
      stop("Window is too long (or too short) for the interpolated series")
    }
    vals <- data[, 2]
    rho <- numeric(m)
    for (m1 in seq_len(m)){
      x <- acy_detrend_linear(vals[m1:(m1 + npts - 1)])
      xm <- x - sum(x) / npts
      rho[m1] <- sum(xm[-1] * xm[-npts]) / sum(xm[-npts]^2)
    }
    xgrid <- pda_linspace(data[1, 1] + window / 2, data[n, 1] - window / 2, m)
    cbind(xgrid, rho)
  }

  # @keywords internal
  # @noRd
  dynot_interp_strict <- function(data, dt){
    # acycleInterpolateSeries.m (linear, integer-index grid)
    x <- data[, 1]
    n <- length(x)
    span <- x[n] - x[1]
    if (dt > span){
      stop("The sampling interval must not exceed the coordinate span")
    }
    raw <- span / dt
    nearest <- pda_round(raw)
    tol <- 8 * pda_eps(max(1, abs(raw)))
    near <- abs(raw - nearest) <= tol
    if (near == TRUE){
      count <- nearest
    } else {
      count <- floor(raw)
    }
    grid <- x[1] + (0:count) * dt
    if (near == TRUE){
      grid[length(grid)] <- x[n]
    }
    if (length(grid) == n && all(grid == x)){
      return(cbind(grid, data[, 2]))
    }
    cbind(grid, dynot_interp1(x, data[, 2], grid))
  }

  # @keywords internal
  # @noRd
  dynot_rho1_strict <- function(data, dt, window){
    # acycleEvolutionaryRho1.m
    if (!(is.finite(dt) && dt > 0)){
      stop("The sampling interval must be finite and positive")
    }
    if (!(is.finite(window) && window > 0)){
      stop("The window must be finite and positive")
    }
    interp <- dynot_interp_strict(data, dt)
    n <- nrow(interp)
    k <- pda_round(window / dt)
    if (k < 3){
      stop("round(window/sample interval) must be at least three points")
    }
    if (k > n){
      stop("The window exceeds the interpolated series")
    }
    n_win <- n - k + 1
    basis <- (1:k) - (k + 1) / 2
    basis_ss <- sum(basis^2)
    rho <- numeric(n_win)
    for (w in seq_len(n_win)){
      v <- interp[w:(w + k - 1), 2]
      cv <- v - sum(v) / k
      slope <- sum(basis * cv) / basis_ss
      dtr <- cv - slope * basis
      z <- dtr - sum(dtr) / k
      den <- sum(z[-k]^2)
      if (!(den > 0)){
        stop(paste0("Window ", w, " has a degenerate rho1 denominator"))
      }
      rho[w] <- sum(z[-1] * z[-k]) / den
    }
    first <- interp[1:n_win, 1]
    last <- interp[k:n, 1]
    center <- first + (last - first) / 2
    out <- cbind(center, rho)
    colnames(out) <- c("center", "rho1")
    out
  }

  # @keywords internal
  # @noRd
  dynot_band_plot <- function(percentile, n_interval, median_index, x_range,
                              x_lab = "depth/time", y_lab = "", reverse = FALSE){
    colorcode <- grDevices::rgb(c(221, 201, 176, 126, 67),
                                c(234, 227, 219, 201, 180),
                                c(224, 209, 188, 146, 100), maxColorValue = 255)
    x <- percentile[, 1]
    pm <- percentile[, -1, drop = FALSE]
    np <- ncol(pm)
    y_lim <- range(pm, na.rm = TRUE)
    if (reverse == TRUE){
      y_lim <- rev(y_lim)
    }
    graphics::plot(x, pm[, median_index], type = "n", xlim = x_range,
                   ylim = y_lim, xlab = x_lab, ylab = y_lab)
    for (i in seq_len(min(n_interval, 5))){
      graphics::polygon(c(x, rev(x)), c(pm[, np + 1 - i], rev(pm[, i])),
                        col = colorcode[i], border = NA)
    }
    graphics::lines(x, pm[, median_index], col = grDevices::rgb(0, 120, 0,
                    maxColorValue = 255), lwd = 1.5, lty = 2)
    invisible(NULL)
  }

  # ---------------------------------------------------------------------------
  # internal helpers
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  pda_eps <- function(x){
    x <- abs(x)
    if (x == 0){
      return(2^-1074)
    }
    2^(floor(log2(x)) - 52)
  }

  # @keywords internal
  # @noRd
  pda_round <- function(x){
    sign(x) * floor(abs(x) + 0.5)
  }

  # @keywords internal
  # @noRd
  pda_linspace <- function(a, b, n){
    if (n < 1){
      return(numeric(0))
    }
    if (n == 1){
      return(b)
    }
    out <- a + (0:(n - 1)) * (b - a) / (n - 1)
    out[1] <- a
    out[n] <- b
    out
  }

  # @keywords internal
  # @noRd
  pda_colon <- function(a, b, d){
    # MATLAB/Octave colon operator a:d:b for d > 0
    q <- (b - a) / d
    n <- floor(q + 3 * .Machine$double.eps * max(1, abs(q))) + 1
    if (n < 1){
      return(numeric(0))
    }
    out <- a + (0:(n - 1)) * d
    if (out[n] > b){
      out[n] <- b
    }
    out
  }

  # @keywords internal
  # @noRd
  pda_interpolate <- function(data, rate){
    # interpolate.m: linear interpolation on data(1,1):rate:data(end,1)
    n <- nrow(data)
    xe <- pda_colon(data[1, 1], data[n, 1], rate)
    ye <- stats::approx(data[, 1], data[, 2], xout = xe, method = "linear",
                        rule = 1, ties = "ordered")$y
    cbind(xe, ye)
  }

  # ---- main code ----

  data <- acy_check_two_column(data)
  n_data <- nrow(data)
  if (n_data < 2){
    stop("data must contain at least two rows")
  }
  if (any(diff(data[, 1]) <= 0)){
    stop("data coordinates must be unique")
  }
  sample_rate <- diff(data[, 1])
  data_length <- data[n_data, 1] - data[1, 1]
  samp95 <- dynot_prctile(sample_rate, 95)

  if (n_sim == 0){
    if (is.null(window_min)){
      if (0.3 * data_length > 400){
        window_min <- 400
      } else {
        window_min <- 0.3 * data_length
      }
    }
    if (is.null(sample_min)){
      sample_min <- samp95
    }
    rho <- dynot_rho1_strict(data, sample_min, window_min)
    if (genplot == TRUE){
      graphics::plot(rho[, 1], rho[, 2], type = "l", xlab = "depth/time",
                     ylab = "RHO in AR(1)")
    }
    return(list(rho1 = rho,
                window_points = pda_round(window_min / sample_min),
                sample_interval = sample_min,
                window = window_min))
  }

  if (!(is.finite(n_sim) && n_sim == trunc(n_sim))){
    stop("n_sim must be an integer")
  }
  if (n_sim < 50){
    stop("Number of simulations is too few (Acycle requires at least 50), try 1000")
  }
  if (is.null(n_out)){
    if (n_data > 1000){
      n_out <- 1000
    } else {
      n_out <- n_data
    }
  }
  if (is.null(window_min) || is.null(window_max)){
    if (0.3 * data_length > 400){
      w1 <- 400
      w2 <- 500
    } else {
      w1 <- 0.3 * data_length
      w2 <- 0.4 * data_length
    }
    if (is.null(window_min)){
      window_min <- w1
    }
    if (is.null(window_max)){
      window_max <- w2
    }
  }
  if (is.null(sample_min)){
    sample_min <- samp95
  }
  if (is.null(sample_max)){
    sample_max <- 2 * samp95
  }

  if (!is.null(seed)){
    set.seed(seed)
  }
  u_sample <- stats::runif(n_sim)
  u_window <- stats::runif(n_sim)
  samplez <- sample_min + (sample_max - sample_min) * u_sample
  window_sim <- window_min + (window_max - window_min) * u_window
  y_grid <- pda_linspace(data[1, 1], data[n_data, 1], n_out)
  powy <- matrix(0, n_out, n_sim)
  u_shift <- rep(NA_real_, n_sim)
  shift <- rep(0, n_sim)

  if (verbose == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (i in seq_len(n_sim)){
    window <- window_sim[i]
    data_even <- pda_interpolate(data, samplez[i])
    rhox <- dynot_erho(data_even, window)
    if (shift_grid == TRUE){
      u_shift[i] <- stats::runif(1)
      shift[i] <- -1 * window / 2 + window * u_shift[i]
    }
    powy[, i] <- dynot_interp1(rhox[, 1] + shift[i], rhox[, 2], y_grid)
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, i)
    }
  }

  percent <- c(2.5, 5, 10, 15.865, 25, 50, 75, 84.135, 90, 95, 97.5)
  powyp <- dynot_row_prctile(powy, percent)
  mask <- !is.na(powyp[, 1])
  percentile <- cbind(y_grid[mask], powyp[mask, , drop = FALSE])
  colnames(percentile) <- c("x", paste0("p", percent))
  med <- cbind(y_grid[mask], powyp[mask, 6])
  colnames(med) <- c("x", "median")

  if (genplot == TRUE){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    dynot_band_plot(percentile, 5, 6, range(data[, 1]), x_lab = "depth/time",
                    y_lab = "RHO in AR(1)", reverse = FALSE)
  }

  list(percentile = percentile,
       median = med,
       percentiles = percent,
       grid = y_grid,
       rho1 = powy,
       draws = list(samplez = samplez, windowz = window_sim, shift = shift,
                    u_sample = u_sample, u_window = u_window,
                    u_shift = u_shift))
}
