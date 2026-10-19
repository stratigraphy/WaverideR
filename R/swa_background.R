#' @title Smoothed Window Average (SWA) spectral background
#'
#' @description Estimates the background of a power spectrum with the
#' Smoothed Window Average (SWA) method of Kemp (2016), as refined by
#' Weedon (2020, 2022). The log10 power spectrum is averaged in
#' non-overlapping windows of odd width (11 to 99 ordinates, always keeping
#' at least three complete windows). The window averages are linearly
#' interpolated, with two alternative extrapolations towards zero frequency.
#' Among all candidate widths/extrapolations the fit with the fewest
#' non-decreasing fitted steps (and, in case of a tie, the smallest RMSE)
#' is selected. A constant (white-noise) and a quadratic high-frequency tail
#' are then tested, and finally the fit is repeatedly Hanning smoothed as
#' long as the RMSE degrades by less than 1\%.
#' Translated from Acycle v3.1 specswa.m (Mingsong Li, after the FORTRAN
#' code of Graham Weedon).
#'
#' @param freq Strictly increasing frequency vector (at least 33 values).
#' @param log_power log10 power values belonging to \code{freq}.
#' @param n_data Number of data points of the analysed series. Retained for
#' compatibility with Acycle; it does not influence the fitted background.
#' @param verbose Print fit diagnostics. \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li,
#' after the SWA FORTRAN program of Graham Weedon.
#'
#' @references
#' Kemp, D.B., 2016. Optimizing significance testing of astronomical forcing
#' in cyclostratigraphy. Paleoceanography 31, 1516-1531.
#'
#' Weedon, G.P., 2020. Spectral background estimation for cyclostratigraphy
#' (SWA). In: Weedon, G.P. 2022, Problems with the current practice of
#' spectral analysis in cyclostratigraphy: Avoiding false detection of
#' regular cyclicity. Earth-Science Reviews 235, 104261.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 512
#' x <- stats::filter(stats::rnorm(n), 0.7, method = "recursive")
#' x <- as.numeric(x) + sin(2 * pi * (1:n) / 20)
#' sp <- stats::spec.pgram(x, taper = 0, detrend = FALSE, plot = FALSE)
#' bg <- swa_background(freq = sp$freq, log_power = log10(sp$spec),
#'   n_data = n)
#' plot(sp$freq, sp$spec, type = "l", log = "y")
#' lines(sp$freq, bg$swa, col = "red", lwd = 2)
#' }
#'
#' @return A list with elements
#' \code{swa} (SWA background in linear power units, same length as
#' \code{freq}), \code{n_windows} (selected window width in ordinates;
#' second output of specswa.m) and \code{diagnostics} (list with
#' \code{initial_variant}, \code{initial_nondecreasing_count},
#' \code{initial_rmse}, \code{selected_window}, \code{tail_model}
#' ("extrapolated", "constant" or "quadratic"), \code{n_smoothing},
#' \code{final_rmse}, \code{n_data}, \code{n_frequency}).
#'
#' @export
#' @importFrom stats sd
swa_background <- function(freq = NULL, log_power = NULL, n_data = NULL,
                           verbose = FALSE){
  # ---- internal helper functions (local to this function) ----

  # @title Quadratic least-squares fit, MATLAB polyfit(x, y, 2) with [p,S,mu]
  # @keywords internal
  # @noRd
  swa_polyfit2 <- function(x, y){
    mu1 <- mean(x)
    mu2 <- stats::sd(x)
    xs <- (x - mu1) / mu2
    V <- cbind(xs^2, xs, 1)
    q <- qr(V)
    p <- qr.coef(q, y)
    as.numeric(V %*% p)
  }

  # ---- main code ----
  freq <- as.numeric(freq)
  xlogp <- as.numeric(log_power)
  if (length(freq) < 1 || any(!is.finite(freq))){
    stop("freq must be a non-empty finite numeric vector")
  }
  if (length(xlogp) < 1 || any(!is.finite(xlogp))){
    stop("log_power must be a non-empty finite numeric vector")
  }
  if (is.null(n_data) || length(n_data) != 1 || !is.finite(n_data) ||
      n_data <= 0 || n_data != round(n_data)){
    stop("n_data must be a positive integer")
  }
  if (any(diff(freq) <= 0)){
    stop("freq must be strictly increasing")
  }
  nout <- length(xlogp)
  if (length(freq) != nout || nout < 33){
    stop(paste0("SWA requires matching frequency/log-power vectors with at ",
                "least 33 ordinates (three initial 11-bin windows)"))
  }

  winwidth <- 9
  rmsemin <- 1.0e+20
  revmin <- Inf
  selected_finish <- NULL
  selected_winwidth <- NULL
  selected_last_mean <- NULL
  selected_variant <- NULL
  selected_rev <- NULL
  selected_rmse <- NULL
  winfinal <- NULL
  fitsp <- numeric(nout)
  fitsp1 <- numeric(nout)
  fitsp2 <- numeric(nout)

  while (winwidth < 99 && floor(nout / (winwidth + 2)) >= 3){
    winwidth <- winwidth + 2
    wincentre <- floor(winwidth / 2) + 1
    ncycle <- floor(nout / winwidth)
    # window log-power averages
    midw <- colSums(matrix(xlogp[seq_len(ncycle * winwidth)],
                           nrow = winwidth)) / winwidth

    # low frequencies: project slope(s) back towards zero frequency
    diff1 <- (midw[2] - midw[1]) / winwidth
    diff2 <- (midw[3] - midw[2]) / winwidth
    delta <- diff2 - diff1
    diff0 <- diff1 - delta
    steps <- wincentre - (1:(wincentre - 1))
    fitsp1[1:(wincentre - 1)] <- midw[1] - steps * diff1
    fitsp2[1:(wincentre - 1)] <- midw[1] - steps * diff0
    fitsp1[wincentre] <- midw[1]
    fitsp2[wincentre] <- midw[1]

    # interpolation between window averages
    ii <- wincentre
    for (k in seq_len(ncycle - 1)){
      d0 <- (midw[k + 1] - midw[k]) / winwidth
      idx <- (ii + 1):(ii + winwidth)
      val <- midw[k] + (1:winwidth) * d0
      fitsp1[idx] <- val
      fitsp2[idx] <- val
      ii <- ii + winwidth
    }

    # high frequencies: slope(s) of the last windows
    diff1 <- (midw[ncycle] - midw[ncycle - 1]) / winwidth
    diff2 <- (midw[ncycle - 1] - midw[ncycle - 2]) / winwidth
    delta <- diff2 - diff1
    diff0 <- diff1 - delta
    finish <- ii
    if (finish < nout){
      steps <- 1:(nout - finish)
      fitsp1[(finish + 1):nout] <- midw[ncycle] + steps * diff1
      fitsp2[(finish + 1):nout] <- midw[ncycle] + steps * diff0
    }

    rmsef1 <- sqrt(sum((xlogp - fitsp1)^2) / (nout - 1))
    rmsef2 <- sqrt(sum((xlogp - fitsp2)^2) / (nout - 1))
    rev1 <- sum(diff(fitsp1) >= 0)
    rev2 <- sum(diff(fitsp2) >= 0)

    if (rev1 < revmin || (rev1 == revmin && rmsef1 < rmsemin)){
      revmin <- rev1
      rmsemin <- rmsef1
      winfinal <- winwidth
      selected_finish <- finish
      selected_winwidth <- winwidth
      selected_last_mean <- midw[ncycle]
      selected_variant <- 1
      selected_rev <- rev1
      selected_rmse <- rmsef1
      fitsp <- fitsp1
    }
    if (rev2 < revmin || (rev2 == revmin && rmsef2 < rmsemin)){
      revmin <- rev2
      rmsemin <- rmsef2
      winfinal <- winwidth
      selected_finish <- finish
      selected_winwidth <- winwidth
      selected_last_mean <- midw[ncycle]
      selected_variant <- 2
      selected_rev <- rev2
      selected_rmse <- rmsef2
      fitsp <- fitsp2
    }
  }

  if (is.null(selected_finish)){
    stop("SWA could not select a finite background fit")
  }
  finish <- selected_finish
  winwidth <- selected_winwidth
  if (verbose == TRUE){
    cat(sprintf("\nwinfinal = %d\nrmsemin = %f\n", winfinal, rmsemin))
  }
  rmsel1 <- rmsemin
  tail_model <- "extrapolated"

  # constant (white noise) end of spectrum
  fitsp3 <- fitsp
  fitsp3[finish:nout] <- selected_last_mean
  rmse <- sqrt(sum((xlogp - fitsp3)^2) / (nout - 1))
  if (rmse < rmsel1){
    fitsp <- fitsp3
    rmsel1 <- rmse
    tail_model <- "constant"
    if (verbose == TRUE){
      cat(sprintf("  Using fixed final power in end window, rmse = %f\n",
                  rmsel1))
    }
  } else {
    if (verbose == TRUE){
      cat(sprintf("  White noise end of log spectrum: RMSE = %f (not better)\n",
                  rmse))
    }
  }

  # quadratic end fit (MATLAB polyfit with centring and scaling)
  tail_start <- max(1, finish - floor(winwidth / 2))
  tail_index <- tail_start:nout
  x_tail <- freq[tail_index]
  y_tail <- xlogp[tail_index]
  fitsp3 <- fitsp
  if (length(tail_index) >= 3 && length(unique(x_tail)) >= 3){
    fitsp3[tail_index] <- swa_polyfit2(x_tail, y_tail)
    rmse <- sqrt(sum((xlogp - fitsp3)^2) / (nout - 1))
  } else {
    rmse <- Inf
  }
  if (rmse < rmsel1){
    fitsp <- fitsp3
    rmsel1 <- rmse
    tail_model <- "quadratic"
    if (verbose == TRUE){
      cat(sprintf("  Switching to quadratic fit for end window, rmse = %f\n",
                  rmsel1))
    }
  } else {
    if (verbose == TRUE){
      cat(sprintf("  Quadratic fit for end window: RMSE = %f (not better)\n",
                  rmse))
    }
  }

  # Hanning smoothing while RMSE degrades by less than 1 %
  nsmth <- 0
  rmset <- rmsel1
  while (nsmth < 500){
    candidate <- numeric(nout)
    candidate[1] <- 0.5 * fitsp[1] + 0.5 * fitsp[2]
    candidate[nout] <- 0.5 * fitsp[nout - 1] + 0.5 * fitsp[nout]
    mid <- 2:(nout - 1)
    candidate[mid] <- 0.25 * fitsp[mid - 1] + 0.5 * fitsp[mid] +
      0.25 * fitsp[mid + 1]
    candidate_rmse <- sqrt(sum((xlogp - candidate)^2) / (nout - 1))
    if (candidate_rmse < rmsel1 * 1.01){
      fitsp <- candidate
      rmset <- candidate_rmse
      nsmth <- nsmth + 1
    } else {
      break
    }
  }
  if (verbose == TRUE){
    cat(sprintf("\nAfter smoothing RMSE = %.6f\n", rmset))
  }

  diagnostics <- list(initial_variant = selected_variant,
                      initial_nondecreasing_count = selected_rev,
                      initial_rmse = selected_rmse,
                      selected_window = winfinal,
                      tail_model = tail_model,
                      n_smoothing = nsmth,
                      final_rmse = rmset,
                      n_data = n_data,
                      n_frequency = nout)
  list(swa = 10^fitsp, n_windows = winfinal, diagnostics = diagnostics)
}
