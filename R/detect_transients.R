#' Transient (impulsive) features of a record by transient extraction
#'
#' Finds the positions in a record where energy over many periods arrives at the same instant: the signature of an
#' impulsive feature such as a hiatus (a phase jump of all cycles), an erosion surface, a turbidite or a spike. The
#' wavelet coefficients are reassigned along the record to their group delay (the time-reassigned synchrosqueezing /
#' transient-extracting transform, He et al. 2019; Yu et al. 2023), only the coefficients that lie on their own group
#' delay are kept, and their energy is summed over the periods. The sum is compared with a running median and peaks
#' above \code{threshold} times the median are reported, which can be passed as \code{jump_depths} to
#' \code{\link{track_astro_hmm}} or \code{\link{track_astro_auto}}.
#'
#' The transform cannot see a jump of the sedimentation rate itself (the cycles stay continuous); it sees the phase
#' jumps that a gap in time causes.
#'
#' @param data Record: first column depth, second column proxy.
#' @param lowerPeriod,upperPeriod,dj Period range and spacing of the transform.
#' @param omega_nr Morlet omega_0 \code{Default=6} (short wavelets localise transients best).
#' @param threshold Peaks above this multiple of the running median of the transient energy are reported
#' \code{Default=4}.
#' @param min_distance Minimum distance between reported transients (depth units) \code{Default=upperPeriod / 4}.
#' @param window Length (depth units) of the running median \code{Default=2 * upperPeriod}.
#'
#' @return A list with \code{profile} (data frame: depth, energy, ratio to the running median) and \code{transients}
#' (data frame: depth, ratio of the detected peaks, strongest first).
#'
#' @references He, D., Cao, H., Wang, S. and Chen, X. (2019). Time-reassigned synchrosqueezing transform: the
#' algorithm and its applications in mechanical signal processing. Mech. Syst. Signal Process. 117, 255-279.
#'
#' Yu, G. et al. (2023). Synchro-transient-extracting transform for the analysis of signals with both harmonic and
#' impulsive components. arXiv:2306.01320.
#'
#' @examples
#' set.seed(2)
#' t <- seq(0, 2000, by = 0.5)
#' x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 41) + 0.4 * cos(2 * pi * t / 405)
#' keep <- !(t > 800 & t < 900) & !(t > 1500 & t < 1550)        # two gaps in time (hiatuses)
#' z <- seq(0, by = 2, length.out = sum(keep))
#' x <- x_t[keep] + rnorm(length(z), sd = 0.5)
#' tr <- detect_transients(data.frame(z, x), lowerPeriod = 20, upperPeriod = 400)
#' tr$transients
#' @export
detect_transients <- function(data, lowerPeriod, upperPeriod, dj = 1 / 20, omega_nr = 6, threshold = 4,
                              min_distance = NULL, window = NULL) {
  dat <- as.matrix(data)
  dat <- dat[stats::complete.cases(dat[, 1:2]), 1:2, drop = FALSE]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dt <- stats::median(diff(dat[, 1]))
  xout <- seq(dat[1, 1], dat[nrow(dat), 1], by = dt)
  xi <- stats::approx(dat[, 1], dat[, 2], xout)$y
  n <- length(xi)
  xs <- (xi - mean(xi)) / stats::sd(xi)
  J <- as.integer(log2(upperPeriod / lowerPeriod) / dj)
  periods <- lowerPeriod * 2^((0:J) * dj)
  wdef <- .wr_wavelet_equiv("morlet", omega_nr, NULL, ngrid = 8001)
  scales <- periods / wdef$fourier_factor
  pad <- 2^(trunc(log2(n) + 0.5) + 1) - n
  N <- n + pad
  omega.k <- 1:floor(N / 2) * (2 * pi) / (N * dt)
  omega.k <- c(0, omega.k, -omega.k[floor((N - 1) / 2):1])
  fx <- fft(c(xs, rep(0, pad)))
  inv <- function(v) (fft(fx * v, inverse = TRUE) / N)[1:n]
  energy <- numeric(n)                      # energy reassigned to its group delay
  orig <- numeric(n)                        # energy at the original position
  for (j in seq_along(scales)) {
    dau <- wdef$daughter(scales[j], omega.k, dt)
    w <- inv(dau)
    h <- 1e-4 / scales[j]
    ddau <- (wdef$daughter(scales[j], omega.k + h, dt) - wdef$daughter(scales[j], omega.k - h, dt)) / (2 * h)
    wt <- inv(1i * ddau)
    tau <- Re(wt / w)                                # b - t_hat
    e <- Mod(w)^2 / scales[j]
    tb <- round((seq_len(n) - 1) - tau / dt) + 1
    ok <- is.finite(tb) & tb >= 1 & tb <= n & is.finite(e)
    acc <- numeric(n); agg <- rowsum(e[ok], tb[ok]); acc[as.integer(rownames(agg))] <- agg[, 1]
    energy <- energy + acc
    orig <- orig + e
  }
  # a harmonic component is reassigned onto itself (t_hat = b): the ratio of the reassigned to the original energy is
  # about one; an impulsive feature draws the energy of its whole neighbourhood onto one position
  if (is.null(window)) window <- 2 * upperPeriod
  k <- max(3, 2 * floor(window / dt / 2) + 1)
  k <- min(k, 2 * floor((n - 1) / 2) + 1)
  base <- stats::runmed(orig, k, endrule = "median")
  ratio <- energy / pmax(base, 1e-12 * max(orig))
  profile <- data.frame(depth = xout, energy = energy, ratio = ratio)
  if (is.null(min_distance)) min_distance <- upperPeriod / 4
  cand <- which(ratio > threshold)
  trans <- data.frame(depth = numeric(0), ratio = numeric(0))
  if (length(cand)) {
    o <- cand[order(-ratio[cand])]
    for (i in o) {
      if (!nrow(trans) || all(abs(xout[i] - trans$depth) >= min_distance)) trans <- rbind(trans, data.frame(depth = xout[i], ratio = ratio[i]))
    }
  }
  list(profile = profile, transients = trans)
}
