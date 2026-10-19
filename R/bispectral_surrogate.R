#' @title FT and IAAFT surrogates for bispectral significance testing
#'
#' @description The \code{\link{bispectral_surrogate}} function generates
#' null-hypothesis surrogate series. Fourier-transform (FT) phase-randomised
#' surrogates preserve the periodogram exactly (all Fourier amplitudes,
#' including DC and the even-length Nyquist bin) while randomising the phases.
#' Iterative amplitude-adjusted Fourier-transform (IAAFT) surrogates
#' (Schreiber and Schmitz, 1996) preserve the observed marginal distribution
#' exactly and iteratively approach the observed Fourier amplitudes; the
#' iteration stops early when the relative spectral error changes by less than
#' 1e-10. The reported relative spectral error is evaluated only on the
#' independent positive frequencies (DC and the even-length Nyquist bin
#' excluded).
#' Translated from Acycle v3.1 bispectralSurrogate.m (Mingsong Li).
#'
#' @param data Input: a numeric vector, or a 2-column matrix [depth/time,
#' value] of which the second column is used.
#' @param type "iaaft" or "ft" (phase randomisation). \code{Default="iaaft"}
#' @param iaaft_iter Maximum number of IAAFT iterations. \code{Default=200}
#' @param n_sim Number of surrogates to generate. \code{Default=1}
#' @param seed Seed for the random number generator (\code{NULL} leaves the
#' generator untouched). \code{Default=NULL}
#' @param verbose Show a progress bar. \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' (Acycle bispectral module).
#'
#' @references
#' Theiler, J., Eubank, S., Longtin, A., Galdrikian, B., Farmer, J.D., 1992.
#' Testing for nonlinearity in time series: the method of surrogate data.
#' Physica D 58, 77-94.
#'
#' Schreiber, T., Schmitz, A., 1996. Improved surrogate data for nonlinearity
#' tests. Physical Review Letters 77, 635-638.
#'
#' @examples
#' \donttest{
#' set.seed(2)
#' x <- exp(0.35 * as.numeric(arima.sim(list(ar = 0.75), n = 512)))
#' sur <- bispectral_surrogate(data = x, type = "iaaft", n_sim = 5, seed = 3)
#' # marginal distribution preserved exactly
#' max(abs(sort(sur$surrogates[, 1]) - sort(x)))
#' sur$relative_spectral_error
#' }
#' @return
#' A list with \code{surrogates} (matrix, one surrogate per column),
#' \code{type}, \code{iterations} (completed IAAFT iterations per surrogate,
#' 1 for FT) and \code{relative_spectral_error} (per surrogate).
#'
#' @export
#' @importFrom stats fft
#' @importFrom stats runif
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar

bispectral_surrogate <- function(data = NULL,
                                 type = "iaaft",
                                 iaaft_iter = 200,
                                 n_sim = 1,
                                 seed = NULL,
                                 verbose = FALSE){
  # ---- internal helper functions (local to this function) ----

  # @keywords internal
  # @noRd
  bisp_canonical_surrogate <- function(value){
    value <- tolower(trimws(value))
    if (value %in% c("phase", "ft", "fourier")){
      return("ft")
    }
    if (value == "iaaft"){
      return("iaaft")
    }
    stop("surrogate type must be ft (phase) or iaaft.")
  }

  # ---------------------------------------------------------------------------
  # surrogates (bispectralSurrogate.m)
  # ---------------------------------------------------------------------------

  # @keywords internal
  # @noRd
  bisp_draw <- function(n, type){
    # random numbers consumed by one surrogate (same order as MATLAB code)
    if (type == "ft"){
      stats::runif(max(0, floor((n + 1) / 2) - 1))
    } else {
      sample.int(n)
    }
  }

  # @keywords internal
  # @noRd
  bisp_spectral_error <- function(y, target){
    n <- length(y)
    if (n < 3){
      return(0)
    }
    pos <- 2:floor((n + 1) / 2)
    cand <- Mod(stats::fft(y))
    d <- cand[pos] - target[pos]
    sqrt(sum(d^2)) / max(sqrt(sum(target[pos]^2)), .Machine$double.xmin)
  }

  # @keywords internal
  # @noRd
  bisp_surrogate_core <- function(y, type, iterations, draw){
    n <- length(y)
    if (type == "ft"){
      spec <- stats::fft(y)
      target <- Mod(spec)
      npos <- max(0, floor((n + 1) / 2) - 1)
      pos <- 1 + seq_len(npos)
      if (length(pos) > 0){
        spec[pos] <- Mod(spec[pos]) * exp(1i * 2 * pi * draw)
        spec[n - pos + 2] <- Conj(spec[pos])
      }
      sur <- Re(stats::fft(spec, inverse = TRUE)) / n
      return(list(surrogate = sur, iterations = 1,
                  relative_spectral_error = bisp_spectral_error(sur, target)))
    }
    target <- Mod(stats::fft(y))
    sorted_target <- sort(y)
    sur <- y[draw]
    prev <- Inf
    completed <- 0
    err <- NA_real_
    for (i in seq_len(iterations)){
      ph <- Arg(stats::fft(sur))
      proj <- Re(stats::fft(target * exp(1i * ph), inverse = TRUE)) / n
      rk <- order(proj, method = "radix")
      nxt <- numeric(n)
      nxt[rk] <- sorted_target
      sur <- nxt
      err <- bisp_spectral_error(sur, target)
      completed <- i
      if (i > 1 && abs(prev - err) <= 1e-10 * max(1, prev)){
        break
      }
      prev <- err
    }
    list(surrogate = sur, iterations = completed, relative_spectral_error = err)
  }

  # ---- main code ----
  if (is.null(data)){
    stop("data must be a numeric vector or a 2-column matrix")
  }
  if (is.matrix(data) || is.data.frame(data)){
    data <- as.matrix(data)
    if (ncol(data) >= 2){
      y <- as.numeric(data[, 2])
    } else {
      y <- as.numeric(data[, 1])
    }
  } else {
    y <- as.numeric(data)
  }
  if (length(y) == 0 || any(!is.finite(y))){
    stop("Bispectral surrogates require a nonempty finite input series.")
  }
  type <- bisp_canonical_surrogate(type)
  if (type == "iaaft"){
    if (!(length(iaaft_iter) == 1 && is.finite(iaaft_iter) &&
          iaaft_iter >= 1 && iaaft_iter == round(iaaft_iter))){
      stop("iaaft_iter must be a positive integer.")
    }
  }
  if (!is.null(seed)){
    set.seed(seed)
  }
  n <- length(y)
  surrogates <- matrix(NA_real_, n, n_sim)
  iterations <- numeric(n_sim)
  errors <- numeric(n_sim)
  if (verbose == TRUE){
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  }
  for (i in seq_len(n_sim)){
    draw <- bisp_draw(n, type)
    s <- bisp_surrogate_core(y, type, iaaft_iter, draw)
    surrogates[, i] <- s$surrogate
    iterations[i] <- s$iterations
    errors[i] <- s$relative_spectral_error
    if (verbose == TRUE){
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (verbose == TRUE){
    close(pb)
  }
  list(surrogates = surrogates,
       type = type,
       iterations = iterations,
       relative_spectral_error = errors)
}
