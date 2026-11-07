#' @title Bayesian change point analysis (Ruggieri, 2013)
#'
#' @description The \code{\link{bayes_changepoint}} function applies the
#' Bayesian change point algorithm of Ruggieri (2013) to a proxy record. The
#' record is modelled as a sequence of segments, each described by its own
#' linear regression (intercept + trend) with conjugate normal /
#' scaled-inverse-chi-square priors. The algorithm
#' \enumerate{
#' \item calculates the marginal likelihood (log form) of every possible
#' sub-interval of the data (Eq. 2 of Ruggieri, 2013),
#' \item builds the partition function by forward recursion (dynamic
#' programming, Eqs. 3-4),
#' \item calculates the exact posterior distribution of the number of change
#' points (Eqs. 5-7) and
#' \item draws \code{n_sim} solutions (number of change points, their
#' locations and the segment regressions) by stochastic backtrace.
#' }
#' The posterior probability of a change point at each data point, the
#' posterior-mean regression model and its \eqn{R^2} are returned as in
#' Acycle. In addition the per-sample models are kept, from which pointwise
#' credible intervals of the model are derived (these reflect the
#' uncertainty in the number and location of the change points; each sampled
#' segment uses its posterior-mean regression coefficients, as in Acycle).
#'
#' Translated from Acycle v3.1 bayes_cpt.m, partition_fn.m, pick_k1.m and
#' nCk.m (Mingsong Li, modified from the code of Eric Ruggieri); default
#' hyperparameters as in the Acycle "Changepoint" menu (AC.m).
#'
#' @param data Record to analyse (matrix or data frame with two columns: depth or time in column 1, proxy value
#' in column 2; finite values only, at least two rows). Data points are treated as equally spaced (the
#' regression uses the sample index as predictor, as in Acycle). Required, no default (\code{Default=NULL} gives
#' an error).
#' @param k_max Maximum number of change points considered (integer, >= 1). \code{Default=10}.
#' @param d_min Minimum distance in data points between adjacent change points, i.e. the minimum segment length
#' (integer, >= 1). \code{Default=1}.
#' @param k_0 Hyperparameter (variance scaling) of the normal prior on the regression coefficients (numeric,
#' > 0); smaller values give a vaguer prior. \code{Default=0.01}.
#' @param v_0 Hyperparameter (number of pseudo data points) of the scaled inverse chi-square prior on the
#' segment variance (numeric, > 0); \code{NULL} uses \code{ceiling(abs(x[n] - x[1]) / 44)}, with x the first
#' column of \code{data}, as in Acycle. \code{Default=NULL}.
#' @param sig_0 Hyperparameter (pseudo data variance) of the scaled inverse chi-square prior on the segment
#' variance (numeric, > 0); \code{NULL} uses the variance of the proxy. \code{Default=NULL}.
#' @param n_sim Number of solutions (number and positions of the change points with their segment regressions)
#' drawn by stochastic backtrace (integer, >= 1). \code{Default=500}.
#' @param seed Seed of the random number generator (integer), so that the sampled solutions can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param ci Lower and upper probabilities of the pointwise credible interval of the model (numeric vector of
#' length 2, between 0 and 1). \code{Default=c(0.025,0.975)}.
#' @param stable Use numerically safe log-space arithmetic (shifted log-sum-exp, \code{lgamma}, log binomial
#' coefficients) (logical, TRUE or FALSE). With \code{FALSE} the arithmetic of Acycle is reproduced exactly
#' (\code{log(sum(exp(.)))} in linear space), which under- or overflows for long records or records with a very
#' small or large variance; a warning is given when that happens. \code{Default=FALSE}.
#' @param genplot Plot the result with \code{\link{plot_bayes_changepoint}} (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param verbose Show a text progress bar while the solutions are sampled (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#'
#' @author Michiel Arts (R translation); original MATLAB code: Eric Ruggieri,
#' modified for Acycle by Mingsong Li.
#'
#' @references
#' Ruggieri, E. (2013). A Bayesian approach to detecting change points in
#' climatic records. International Journal of Climatology, 33, 520-528.
#' <doi:10.1002/joc.3447>
#'
#' Li, M., Hinnov, L., Kump, L. (2019). Acycle: Time-series analysis software
#' for paleoclimate research and education. Computers & Geosciences, 127,
#' 12-22. <doi:10.1016/j.cageo.2019.02.011>
#'
#' @examples
#' \donttest{
#' # a record with two changes in mean and one change in trend
#' set.seed(2)
#' x <- 1:150
#' y <- c(rep(0, 50), rep(2, 50), 2 + 0.05 * (1:50)) + stats::rnorm(150, sd = 0.4)
#' bcp <- bayes_changepoint(data = cbind(x, y), k_max = 6, d_min = 5,
#'                          k_0 = 0.01, v_0 = 2, n_sim = 500, seed = 1)
#' bcp$k_posterior
#' plot_bayes_changepoint(bcp)
#' }
#'
#' @return A list with elements
#' \itemize{
#' \item \code{model}: matrix [x, posterior-mean model] (Acycle "BayesRegModel").
#' \item \code{cpt}: matrix [x, posterior probability of a change point]
#' (Acycle "BayesChangepoint"). A change point at index i means that a new
#' segment starts at index i+1.
#' \item \code{r2}: \eqn{R^2 = 1 - \sum(y-model)^2 / \sum y^2} as in Acycle.
#' \item \code{k_posterior}: matrix [number of change points (0..k_max),
#' posterior probability].
#' \item \code{model_ci}: matrix [x, lower, median, upper] pointwise credible
#' interval of the sampled models.
#' \item \code{beta}: matrix [intercept, slope] of the posterior-mean
#' regression coefficients at each data point (predictor = sample index).
#' \item \code{samples}: \code{n_sim} x \code{k_max} matrix with the index of
#' the sampled change points (0 = none), as Acycle's samp_holder.
#' \item \code{samples_x}: same as \code{samples} expressed in the unit of
#' the first column (NA = none).
#' \item \code{n_cpt}: sampled number of change points for every solution.
#' \item \code{sample_models}: N x \code{n_sim} matrix of the sampled models.
#' \item \code{log_partition}: the k_max x N partition function (log form).
#' \item \code{param}: list with the used hyperparameters.
#' }
#'
#' @export
#' @importFrom stats var
#' @importFrom stats runif
#' @importFrom stats quantile
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
bayes_changepoint <- function(data = NULL,
    k_max = 10,
    d_min = 1,
    k_0 = 0.01,
    v_0 = NULL,
    sig_0 = NULL,
    n_sim = 500,
    seed = NULL,
    ci = c(0.025, 0.975),
    stable = FALSE,
    genplot = FALSE,
    verbose = FALSE){
  # Internal helper functions (local to this function)

  # internal helpers

  # @keywords internal
  # @noRd
  bcp_prepare <- function(data, k_max, d_min, k_0, v_0, sig_0, n_sim){
    data <- as.matrix(data)
    if (ncol(data) < 2){
      stop("data must be a matrix with two columns: [depth/time, value]")
    }
    x <- as.numeric(data[, 1])
    y <- as.numeric(data[, 2])
    if (any(!is.finite(x)) || any(!is.finite(y))){
      stop("data must be finite")
    }
    nn <- length(y)
    if (nn < 2){
      stop("data must contain at least two rows")
    }
    if (is.null(v_0)){
      v_0 <- ceiling(abs(x[nn] - x[1]) / 44)
    }
    if (is.null(sig_0)){
      sig_0 <- stats::var(y)
    }
    if (k_max < 1 || k_max != round(k_max)){
      stop("k_max must be a positive integer")
    }
    if (d_min < 1 || d_min != round(d_min)){
      stop("d_min must be a positive integer")
    }
    if (n_sim < 1 || n_sim != round(n_sim)){
      stop("n_sim must be a positive integer")
    }
    if (k_0 <= 0 || v_0 <= 0 || sig_0 <= 0){
      stop("k_0, v_0 and sig_0 must be positive")
    }
    list(x = x, y = y, k_max = k_max, d_min = d_min, k_0 = k_0, v_0 = v_0,
        sig_0 = sig_0, n_sim = n_sim)
  }

  # @keywords internal
  # @noRd
  bcp_lse <- function(v, stable = FALSE){
    # log(sum(exp(v))) as in Acycle, or the shifted (safe) version
    if (length(v) == 0){
      return(-Inf)
    }
    if (stable == FALSE){
      return(log(sum(exp(v))))
    }
    mx <- max(v)
    if (mx == -Inf){
      return(-Inf)
    }
    mx + log(sum(exp(v - mx)))
  }

  # @keywords internal
  # @noRd
  bcp_solve2 <- function(j11, j12, j22, b1, b2){
    # J \ b for the symmetric positive definite 2 x 2 matrix J (MATLAB uses a
    # Cholesky factorisation for this case); vectorised
    r11 <- sqrt(j11)
    r12 <- j12 / r11
    r22 <- sqrt(j22 - r12 * r12)
    z1 <- b1 / r11
    z2 <- (b2 - r12 * z1) / r22
    x2 <- z2 / r22
    x1 <- (z1 - x2 * r12) / r11
    list(b1 = x1, b2 = x2)
  }

  # @keywords internal
  # @noRd
  bcp_det2 <- function(j11, j12, j22){
    # det(J) by LU with partial pivoting (LAPACK dgetrf), vectorised
    swap <- abs(j12) > abs(j11)
    d <- numeric(length(j11))
    if (any(!swap)){
      l <- j12[!swap] * (1 / j11[!swap])
      d[!swap] <- j11[!swap] * (j22[!swap] - l * j12[!swap])
    }
    if (any(swap)){
      l <- j11[swap] * (1 / j12[swap])
      d[swap] <- -(j12[swap] * (j12[swap] - l * j22[swap]))
    }
    d
  }

  # @keywords internal
  # @noRd
  bcp_log_py <- function(y, d_min, k_0, v_0, sig_0, stable = FALSE){
    # log probability density of every sub-interval y[i:j] (Eq. 2)
    nn <- length(y)
    m <- 2
    py <- matrix(-Inf, nn, nn)
    if (stable == TRUE){
      lg0 <- lgamma(v_0 / 2)
    } else {
      lg0 <- log(gamma(v_0 / 2))
    }
    const <- v_0 * log(sig_0 * v_0 / 2) / 2 + m * log(k_0) / 2 - lg0
    for (i in seq_len(nn)){
      j0 <- i + d_min - 1
      if (j0 > nn){
        next
      }
      tt <- as.numeric(i:nn)  # double: tt * tt and its cumulative sum overflow integers for long records
      yy <- y[i:nn]
      zz <- yy - y[i]
      uu <- tt - i
      s_y <- cumsum(yy)
      s_ty <- cumsum(tt * yy)
      s_t <- cumsum(tt)
      s_tt <- cumsum(tt * tt)
      s_z <- cumsum(zz)
      s_zz <- cumsum(zz * zz)
      s_uz <- cumsum(uu * zz)
      s_u <- cumsum(uu)
      s_uu <- cumsum(uu * uu)
      sel <- (j0 - i + 1):(nn - i + 1)
      a <- sel
      j11 <- k_0 + a
      j12 <- s_t[sel]
      j22 <- k_0 + s_tt[sel]
      bh <- bcp_solve2(j11, j12, j22, s_y[sel], s_ty[sel])
      b1 <- bh$b1
      b2 <- bh$b2
      g <- b1 + b2 * i - y[i]
      ssr <- s_zz[sel] - 2 * g * s_z[sel] - 2 * b2 * s_uz[sel] + a * g * g +
          2 * g * b2 * s_u[sel] + b2 * b2 * s_uu[sel]
      ssr[ssr < 0] <- 0
      s_n <- v_0 * sig_0 + k_0 * (b1 * b1 + b2 * b2) + ssr
      v_n <- v_0 + a
      dj <- bcp_det2(j11, j12, j22)
      py[i, (j0):nn] <- const + lgamma(v_n / 2) - v_n * log(s_n / 2) / 2 -
          a * log(2 * pi) / 2 - log(dj) / 2
    }
    py
  }

  # @keywords internal
  # @noRd
  bcp_partition <- function(py, k_max, stable = FALSE){
    # forward recursion, partition_fn.m (Eqs. 3 and 4)
    nn <- ncol(py)
    pp <- matrix(-Inf, k_max, nn)
    if (nn >= 2){
      for (j in 2:nn){
        tmp <- py[1, 1:(j - 1)] + py[2:j, j]
        pp[1, j] <- bcp_lse(tmp, stable)
      }
    }
    if (k_max >= 2){
      for (k in 2:k_max){
        if (k + 1 > nn){
          next
        }
        for (j in (k + 1):nn){
          tmp <- pp[k - 1, 1:(j - 1)] + py[2:j, j]
          pp[k, j] <- bcp_lse(tmp, stable)
        }
      }
    }
    pp
  }

  # @keywords internal
  # @noRd
  bcp_nck <- function(n, k, stable = FALSE){
    # nCk.m: prod(n-k+1:n)/factorial(k)
    if (stable == TRUE){
      return(exp(sum(log((n - k + 1):n)) - lfactorial(k)))
    }
    p <- 1
    for (v in (n - k + 1):n){
      p <- p * v
    }
    p / factorial(k)
  }

  # @keywords internal
  # @noRd
  bcp_k_posterior <- function(pp, py, k_max, d_min, stable = FALSE){
    # posterior distribution of the number of change points (Eqs. 5-7)
    nn <- ncol(py)
    k <- pp[, nn]
    for (i in seq_len(k_max)){
      nk <- nn - (i + 1) * d_min + i
      if (nk >= i){
        if (stable == TRUE){
          lnck <- sum(log((nk - i + 1):nk)) - lfactorial(i)
        } else {
          lnck <- log(bcp_nck(nk, i))
        }
        k[i] <- k[i] + log(0.5) - log(k_max) - lnck
      }
    }
    k <- c(py[1, nn] + log(0.5), k)
    total_k <- bcp_lse(k, stable)
    exp(k - total_k)
  }

  # @keywords internal
  # @noRd
  bcp_pick <- function(p){
    # pick_k1.m: sample an index from probability vector p with one uniform
    r <- stats::runif(1)
    cs <- cumsum(p)
    stop_here <- !(r >= cs)
    stop_here[is.na(stop_here)] <- TRUE
    idx <- which(stop_here)
    if (length(idx) == 0){
      return(length(p))
    }
    idx[1]
  }

  # @keywords internal
  # @noRd
  bcp_segment_beta <- function(y, i1, i2, k_0){
    tt <- as.numeric(i1:i2)
    yy <- y[i1:i2]
    a <- i2 - i1 + 1
    bcp_solve2(k_0 + a, sum(tt), k_0 + sum(tt * tt), sum(yy), sum(tt * yy))
  }

  # @keywords internal
  # @noRd
  bcp_sample <- function(y, py, pp, kprob, n_sim, k_0, stable = FALSE,
      keep_models = TRUE, verbose = FALSE){
    # stochastic backtrace (bayes_cpt.m step 5)
    nn <- length(y)
    k_max <- nrow(pp)
    samp_holder <- matrix(0, n_sim, k_max)
    chgpt_loc <- numeric(nn)
    beta1 <- numeric(nn)
    beta2 <- numeric(nn)
    n_cpt <- numeric(n_sim)
    models <- NULL
    if (keep_models == TRUE){
      models <- matrix(0, nn, n_sim)
    }
    tt_all <- 1:nn
    if (verbose == TRUE){
      pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    }
    for (s in seq_len(n_sim)){
      num <- bcp_pick(kprob) - 1
      n_cpt[s] <- num
      if (num > 0){
        back <- nn
        if (num >= 2){
          for (kk in num:2){
            tmp <- pp[kk - 1, 1:(back - 1)] + py[2:back, back]
            tot <- bcp_lse(tmp, stable)
            tmp <- exp(tmp - tot)
            cp <- bcp_pick(tmp)
            chgpt_loc[cp] <- chgpt_loc[cp] + 1
            samp_holder[s, kk] <- cp
            bh <- bcp_segment_beta(y, cp + 1, back, k_0)
            idx <- (cp + 1):back
            beta1[idx] <- beta1[idx] + bh$b1
            beta2[idx] <- beta2[idx] + bh$b2
            if (keep_models == TRUE){
              models[idx, s] <- bh$b1 + tt_all[idx] * bh$b2
            }
            back <- cp
          }
        }
        tmp <- py[1, 1:(back - 1)] + py[2:back, back]
        tot <- bcp_lse(tmp, stable)
        tmp <- exp(tmp - tot)
        cp <- bcp_pick(tmp)
        chgpt_loc[cp] <- chgpt_loc[cp] + 1
        samp_holder[s, 1] <- cp
        bh <- bcp_segment_beta(y, cp + 1, back, k_0)
        idx <- (cp + 1):back
        beta1[idx] <- beta1[idx] + bh$b1
        beta2[idx] <- beta2[idx] + bh$b2
        if (keep_models == TRUE){
          models[idx, s] <- bh$b1 + tt_all[idx] * bh$b2
        }
        bh <- bcp_segment_beta(y, 1, cp, k_0)
        idx <- 1:cp
        beta1[idx] <- beta1[idx] + bh$b1
        beta2[idx] <- beta2[idx] + bh$b2
        if (keep_models == TRUE){
          models[idx, s] <- bh$b1 + tt_all[idx] * bh$b2
        }
      } else {
        bh <- bcp_segment_beta(y, 1, nn, k_0)
        beta1 <- beta1 + bh$b1
        beta2 <- beta2 + bh$b2
        if (keep_models == TRUE){
          models[, s] <- bh$b1 + tt_all * bh$b2
        }
      }
      if (verbose == TRUE){
        utils::setTxtProgressBar(pb, s)
      }
    }
    if (verbose == TRUE){
      close(pb)
    }
    beta1 <- beta1 / n_sim
    beta2 <- beta2 / n_sim
    chgpt_loc <- chgpt_loc / n_sim
    model <- beta1 + tt_all * beta2
    r2 <- 1 - sum((y - model)^2) / sum(y^2)
    list(samp_holder = samp_holder, chgpt_loc = chgpt_loc, beta1 = beta1,
        beta2 = beta2, model = model, r2 = r2, n_cpt = n_cpt, models = models)
  }

  # Main code
  prep <- bcp_prepare(data, k_max, d_min, k_0, v_0, sig_0, n_sim)
  if (!is.null(seed)){
    set.seed(seed)
  }
  py <- bcp_log_py(prep$y, prep$d_min, prep$k_0, prep$v_0, prep$sig_0,
      stable = stable)
  pp <- bcp_partition(py, prep$k_max, stable = stable)
  kprob <- bcp_k_posterior(pp, py, prep$k_max, prep$d_min, stable = stable)
  if (any(!is.finite(kprob))){
    warning(paste0("Non-finite posterior probabilities: the Acycle ",
        "(linear-space) arithmetic under/overflowed. Use ",
        "stable = TRUE."))
  }
  smp <- bcp_sample(prep$y, py, pp, kprob, prep$n_sim, prep$k_0,
      stable = stable, keep_models = TRUE, verbose = verbose)

  x <- prep$x
  nn <- length(x)
  model_ci <- matrix(NA, nn, 4)
  model_ci[, 1] <- x
  probs <- c(ci[1], 0.5, ci[2])
  for (i in seq_len(nn)){
    row_i <- smp$models[i, ]
    if (all(is.finite(row_i))){
      model_ci[i, 2:4] <- stats::quantile(row_i, probs = probs, names = FALSE)
    }
  }
  colnames(model_ci) <- c("x", "lower", "median", "upper")

  samples_x <- matrix(NA, prep$n_sim, prep$k_max)
  for (s in seq_len(prep$n_sim)){
    for (kk in seq_len(prep$k_max)){
      if (smp$samp_holder[s, kk] > 0){
        samples_x[s, kk] <- x[smp$samp_holder[s, kk]]
      }
    }
  }
  model <- cbind(x = x, model = smp$model)
  cpt <- cbind(x = x, probability = smp$chgpt_loc)
  k_post <- cbind(n_changepoints = 0:prep$k_max, probability = kprob)
  beta <- cbind(intercept = smp$beta1, slope = smp$beta2)

  out <- list(model = model,
      cpt = cpt,
      r2 = smp$r2,
      k_posterior = k_post,
      model_ci = model_ci,
      beta = beta,
      samples = smp$samp_holder,
      samples_x = samples_x,
      n_cpt = smp$n_cpt,
      sample_models = smp$models,
      log_partition = pp,
      data = cbind(x = x, y = prep$y),
      param = list(k_max = prep$k_max, d_min = prep$d_min,
          k_0 = prep$k_0, v_0 = prep$v_0,
          sig_0 = prep$sig_0, n_sim = prep$n_sim,
          ci = ci, stable = stable))
  if (genplot == TRUE){
    plot_bayes_changepoint(out)
  }
  out
}
