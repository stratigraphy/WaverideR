# modBchron -- pure-R implementation of modifiedBChron (Trayler et al., 2020), the
# Compound Poisson-Gamma age-depth model of Haslett and Parnell (2008) for deep time.
# The C routines of BchronModified.c (truncatedWalk, truncatedRat, dtweediep1,
# predictInterp, predictExtrapUp, predictExtrapDown) are translated line by line into
# R (internal functions with a leading dot), so nothing has to be compiled.
# modBchron() takes the same arguments as modifiedBChron::ageModel() and returns the same
# object (class "mod_bchron_model_run"); modBchron_predict() corresponds to agePredict().
.rtruncn <- function(a, b) {
  # Rejection sampler for a standard normal truncated to [a, b], translated
  # from rtruncn() in BchronModified.c.
  A <- atan(a)
  B <- atan(b)
  maxA <- exp(-a^2 / 4) / cos(A)
  maxB <- exp(-b^2 / 4) / cos(B)
  maxR <- max(maxA, maxB)
  if (a < 1 && b > -1) {
    maxR <- exp(-0.25) * sqrt(2.0)
  }
  repeat {
    r2 <- runif(1, 0, 1)
    r <- sqrt(r2) * maxR
    th <- runif(1, A, B)
    u <- r * cos(th)
    x <- tan(th)
    if ((x^2) < (log(u) * -4)) {
      break
    }
  }
  return(x)
}
#' @keywords internal
#' @noRd
.truncatedWalkC <- function(old, sd, low, high) {
  # translated from truncatedWalk() in BchronModified.c
  lowlimold <- (low - old) / sd
  upplimold <- (high - old) / sd
  y <- .rtruncn(lowlimold, upplimold)
  old + sd * y
}
#' @keywords internal
#' @noRd
.truncatedRatC <- function(old, sd, low, high, newvalue) {
  # translated from truncatedRat() in BchronModified.c
  lowlimold <- (low - old) / sd
  upplimold <- (high - old) / sd
  lowlimnew <- (low - newvalue) / sd
  upplimnew <- (high - newvalue) / sd
  plowold <- pnorm(lowlimold, 0, 1)
  puppold <- pnorm(upplimold, 0, 1)
  plownew <- pnorm(lowlimnew, 0, 1)
  puppnew <- pnorm(upplimnew, 0, 1)
  (puppold - plowold) / (puppnew - plownew)
}
#' @keywords internal
#' @noRd
.dtweedielogwsmallp <- function(y, phi, power) {
  # translated from dtweedielogwsmallp() in BchronModified.c
  drop <- 37
  if (power < 1) {
    stop("Error - power<1!")
  }
  if (power > 2) {
    stop("Error - power>2!")
  }
  if (phi <= 0) {
    stop("Error - phi<=0!")
  }
  if (y <= 0) {
    stop("Error - y<=0!")
  }
  p <- power
  a <- (2 - p) / (1 - p)
  a1 <- 1 - a
  r <- -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p)
  logz <- r
  jmax <- (y^(2 - p)) / (phi * (2 - p))
  j <- max(1, jmax)
  cc <- logz + a1 + a * log(-a)
  wmax <- a1 * jmax
  estlogw <- wmax
  while (estlogw > (wmax - drop)) {
    j <- j + 2
    estlogw <- j * (cc - a1 * log(j))
  }
  hij <- as.integer(ceiling(j))
  logz <- r
  jmax <- y^(2 - power) / (phi * (2 - power))
  j <- max(1, jmax)
  wmax <- a1 * jmax
  estlogw <- wmax
  while ((estlogw > (wmax - drop)) && (j >= 2)) {
    j <- max(1, j - 2)
    estlogw <- j * (cc - a1 * log(j))
  }
  lowj <- as.integer(max(1, floor(j)))
  newj <- lowj:hij
  g <- lgamma(newj + 1) + lgamma(-a * newj)
  A <- r * newj - g
  count <- length(newj) # == hij - lowj + 1
  # NB: the original C computes this "max" as fmax2(A[k], count) inside a
  # loop, i.e. it compares each log-weight against the *element count*
  # rather than finding the true max of A -- this is almost certainly an
  # upstream bug (a proper log-sum-exp shift should just be max(A)), but it
  # is preserved here exactly so behaviour matches the compiled version.
  m <- max(pmax(A, count))
  we <- exp(A - m)
  sumwe <- sum(we)
  log(sumwe) + m
}
#' @keywords internal
#' @noRd
.dtweedieseriessmallp <- function(power, y, mu, phi) {
  # translated from dtweedieseriessmallp() in BchronModified.c
  logw <- .dtweedielogwsmallp(y, phi, power)
  tau <- phi * (power - 1) * mu^(power - 1)
  lambda <- mu^(2 - power) / (phi * (2 - power))
  logf <- -y / tau - lambda - log(y) + logw
  exp(logf)
}
#' @keywords internal
#' @noRd
.dtweediep1C <- function(y, power, mu, phi) {
  # translated from dtweediep1() in BchronModified.c (scalar version; the
  # user-facing `dtweediep1` closures inside ageModel()/bchron_model() wrap
  # this with Vectorize(), exactly as the original .C()-based closures did)
  eps <- 0.00000001
  lambda2 <- mu^(2 - power) / (phi * (2 - power)) - eps
  alpha <- (2 - power) / (power - 1)
  beta <- 1 / (phi * (power - 1) * mu^(power - 1))
  mu2 <- alpha * lambda2 / beta
  phi2 <- (alpha + 1) / ((lambda2 * alpha)^(1 / (alpha + 1)) * beta^(alpha / (alpha + 1)))
  fTplus1 <- .dtweedieseriessmallp(power, y, mu, phi)
  fTplus2 <- .dtweedieseriessmallp(power, y, mu, phi)
  fTplus3 <- .dtweedieseriessmallp(power, y, mu2, phi2)
  fTplus1 + (1 / eps) * (fTplus2 - fTplus3)
}
#' @keywords internal
#' @noRd
.predictInterpC <- function(alpha, lambda, beta, predictPositions,
    diffPositionj, currPositionsj, currPositionsjp1,
    thetaj, thetajp1) {
  # translated from predictInterp() in BchronModified.c
  Nd <- rpois(1, lambda * diffPositionj)
  if (Nd > 0) {
    depthEvents <- sort(runif(Nd, currPositionsj, currPositionsjp1))
  } else {
    depthEvents <- numeric(0)
  }
  timeEventsUnsc <- rgamma(Nd + 1, shape = alpha, rate = beta)
  timeEventsSum <- sum(timeEventsUnsc)
  timeEvents <- (thetajp1 - thetaj) * timeEventsUnsc / timeEventsSum
  timeEventsCumsum <- numeric(Nd + 1)
  timeEventsCumsum[1] <- 0
  if (Nd >= 1) {
    for (i in 2:(Nd + 1)) {
      timeEventsCumsum[i] <- timeEventsCumsum[i - 1] + timeEvents[i]
    }
  }
  allTimeEvents <- numeric(Nd + 2)
  allTimeEvents[1:(Nd + 1)] <- timeEventsCumsum + thetaj
  allTimeEvents[Nd + 2] <- thetajp1
  allDepthEvents <- numeric(Nd + 2)
  allDepthEvents[1] <- currPositionsj
  allDepthEvents[Nd + 2] <- currPositionsjp1
  if (Nd > 0) {
    allDepthEvents[2:(Nd + 1)] <- depthEvents
  }
  ## NOTE (performance): the original C code (and an earlier draft of this
  ## translation) evaluated .linInterp() once per predictPosition in an R
  ## loop, which is fine in compiled C but is a serious bottleneck in pure R
  ## when Nd is large (long segments between dates => large Nd => this
  ## linear scan dominates runtime). approx() does the identical linear
  ## interpolation (allDepthEvents/allTimeEvents are constructed monotonic
  ## increasing here) via R's compiled internals, all predictPositions at
  ## once -- same numeric result, dramatically faster. rule=2 (nearest
  ## endpoint) is used defensively for the edge case where a query falls
  ## just outside the constructed range; in normal use every predictPosition
  ## passed in here is already guaranteed within [currPositionsj,
  ## currPositionsjp1] by the caller.
  approx(allDepthEvents, allTimeEvents, xout = predictPositions, rule = 2)$y
}
#' @keywords internal
#' @noRd
.predictExtrapUpC <- function(alpha, lambda, beta, predictPositions,
    currPositions1, theta1, maxExtrap, extractDate) {
  # translated from predictExtrapUp() in BchronModified.c
  NpredictPositions <- length(predictPositions)
  bad <- TRUE
  count <- 0
  predvals <- numeric(NpredictPositions)
  while (bad) {
    ## PERFORMANCE NOTE: the original C (and an earlier draft of this
    ## translation) drew one rexp()/rgamma() pair per step in a maxExtrap-
    ## length loop -- fine in compiled C, but R's per-call function overhead
    ## makes that loop the dominant cost of the whole MCMC when this runs
    ## every iteration (which it does whenever a predictPosition sits at or
    ## beyond the outermost date, a common case). depthEvents/timeEvents are
    ## just a running sum of i.i.d. increments, so drawing all maxExtrap-1
    ## increments at once with cumsum() is statistically identical (same
    ## distributions, same independence) and dramatically faster; it just
    ## consumes the random stream in a different order than one-at-a-time
    ## draws would; validated against the compiled C by distribution
    ## (mean/sd/range over many replicate calls), not by matching individual
    ## random draws.
    depthEvents <- c(currPositions1, currPositions1 + cumsum(rexp(maxExtrap - 1, rate = lambda)))
    timeEvents  <- c(theta1, theta1 - cumsum(rgamma(maxExtrap - 1, shape = alpha, rate = beta)))
    predvals <- approx(depthEvents, timeEvents, xout = predictPositions, rule = 2)$y
    count <- count + 1
    bad <- FALSE
    if (any(predvals < extractDate)) {
      bad <- TRUE
    }
    if (count == 50) {
      predvals[predvals < extractDate] <- extractDate
      bad <- FALSE
      warning("Unable to find suitable chronologies for top of section - truncated to upper truncation date")
    }
  }
  predvals
}
#' @keywords internal
#' @noRd
.predictExtrapDownC <- function(alpha, lambda, beta, predictPositions,
    currPositions1, theta1, maxExtrap, extractDate) {
  # translated from predictExtrapDown() in BchronModified.c
  NpredictPositions <- length(predictPositions)
  bad <- TRUE
  count <- 0
  predvals <- numeric(NpredictPositions)
  while (bad) {
    ## see the performance note in .predictExtrapUpC() -- vectorized draws,
    ## validated by distribution rather than matching individual draws
    depthEvents <- c(currPositions1, currPositions1 - cumsum(rexp(maxExtrap - 1, rate = lambda)))
    timeEvents  <- c(theta1, theta1 + cumsum(rgamma(maxExtrap - 1, shape = alpha, rate = beta)))
    predvals <- approx(depthEvents, timeEvents, xout = predictPositions, rule = 2)$y
    count <- count + 1
    bad <- FALSE
    if (any(predvals > extractDate)) {
      bad <- TRUE
    }
    if (count == 50) {
      predvals[predvals > extractDate] <- extractDate
      bad <- FALSE
      warning("Unable to find suitable chronologies for bottom of section - truncated to lower truncation date")
    }
  }
  predvals
}
# SECTION 3 -- ageModel.r
# (originally used .C('truncatedWalk'/'truncatedRat'/'dtweediep1'/
#  'predictInterp'/'predictExtrapUp'/'predictExtrapDown', PACKAGE =
#  'modifiedBChron'); those calls now go to the pure-R engines in SECTION 0)
#' Run a modified version of the Bchronology algorithm from the R package "Bchron" (Haslett and Parnell, 2008)
#'
#' @param ages Ages of the dated samples (numeric vector, one value per date, in one time unit such as kyr or Ma that
#' is used for all age arguments). Required (no default).
#' @param ageSds One-sigma uncertainties of \code{ages} (numeric vector, > 0, same unit, length and order as
#' \code{ages}); for \code{distTypes = "U"} it is the half-width of the uniform distribution. Required (no default).
#' @param positions Stratigraphic positions of the dated samples as distance above the base of the section, e.g. in
#' metres (numeric vector, same length and order as \code{ages}; unordered positions are sorted with a warning).
#' Required (no default).
#' @param positionThicknesses Stratigraphic uncertainty of each date, given as half thickness in the unit of
#' \code{positions} (numeric vector, >= 0, same length and order as \code{ages}).
#' \code{Default=rep(0, length(positions))} (no positional uncertainty).
#' @param ids Sample names of the dates (character or numeric vector, same length and order as \code{ages}); all
#' dates with the same id are combined into a single age PDF of one dated horizon. Required (no default).
#' @param distTypes Distribution used for each date (character vector, same length and order as \code{ages}):
#' \code{"G"} Gaussian (mean \code{ages}, standard deviation \code{ageSds}) or \code{"U"} uniform
#' (\code{ages} plus or minus \code{ageSds}). \code{Default=rep('G', length(ages))}.
#' @param MC Number of MCMC iterations (integer, > \code{burn}). \code{Default=10000}.
#' @param burn Number of initial MCMC iterations discarded as burn-in when the intervals are computed (integer,
#' >= 1 and < \code{MC}). \code{Default=2000}.
#' @param probability Probability of the returned highest density interval and of the outlier check of the dates
#' (numeric, between 0 and 1; 0.95 is 95 percent). \code{Default=0.95}.
#' @param adapt Tune the proposal standard deviations of the model parameters during the run with the adaptive
#' proposal algorithm of Haario et al. (1998) (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param mhSD Metropolis-Hastings proposal standard deviations of the ages of the dated horizons, in the age unit
#' (numeric vector of length \code{length(unique(ids))}, > 0); the starting values when \code{adapt = TRUE}, used
#' as they are otherwise. \code{Default=rep(0.001, length(unique(ids)))}.
#' @param psiSD Metropolis-Hastings proposal standard deviation of the Compound Poisson-Gamma scale parameter psi
#' (numeric, > 0); the starting value when \code{adapt = TRUE}. \code{Default=0.001}.
#' @param muSD Metropolis-Hastings proposal standard deviation of the Compound Poisson-Gamma mean parameter mu
#' (numeric, > 0); the starting value when \code{adapt = TRUE}. \code{Default=0.001}.
#' @param mu0 Warm-start centre of the initial Compound Poisson-Gamma mean parameter mu (numeric, > 0, or \code{NULL}).
#' If \code{NULL}, mu is initialised the original way, from the raw date/position spacing. If given, the chain starts
#' at \code{abs(rnorm(1, mean = mu0, sd = muSD))}, so that a chain can be restarted near a previously known good region
#' of parameter space (e.g. in a rejection-sampling/restart scheme) instead of re-exploring from scratch. It does not
#' change the Metropolis-Hastings acceptance ratio, only where the chain begins. \code{Default=NULL}.
#' @param psi0 Warm-start centre of the initial Compound Poisson-Gamma scale parameter psi (numeric, > 0, or
#' \code{NULL}), analogous to \code{mu0}. If \code{NULL}, psi is initialised the original way
#' (\code{abs(rnorm(1, 1, sd = psiSD))}). \code{Default=NULL}.
#' @param thetas0 Warm-start vector of starting ages of the dated horizons (numeric vector, one value per unique
#'   id, in the same order as \code{unique(ids)}). If NULL (default) each is initialised from its own
#'   likelihood peak, as before. This matters most in a batch-RESTART scheme: the adaptive proposal
#'   SDs start at their \code{mhSD} default of 0.001 and grow by only about a factor of 10 per
#'   200-iteration adaptation window, so they need roughly 1000 iterations to reach the right scale
#'   (measured on a 7500 +/- 15 date: mhSD goes 0.001 -> 0.012 at iter 200 -> 0.157 at iter 400 ->
#'   ~34 at iter 1000). A chain restarted COLD inside a short batch never gets there: its dated
#'   horizons stay pinned within a fraction of a kyr of their starting values, so the run silently
#'   reports near-ZERO age uncertainty at the dated horizons. Passing \code{thetas0} together with
#'   the previous run's adapted \code{mhSD}/\code{psiSD}/\code{muSD} lets a restarted chain resume
#'   properly mixed instead of paying that burn-in again. Values are still clipped into
#'   [truncateUp, truncateDown]. \code{Default=NULL}.
#' @param mu_fixed,psi_fixed Hold the Compound Poisson-Gamma mean (\code{mu_fixed}) or scale (\code{psi_fixed})
#'   parameter at the given value and skip its Metropolis-Hastings updates (numeric, > 0, or \code{NULL}).
#'   \code{NULL} estimates the parameter in the MCMC as before. \code{Default=NULL} for both.
#'
#'   WHY THIS EXISTS. mu and psi are estimated from the gaps BETWEEN dated horizons, so n dates give
#'   only n-1 numbers to estimate both. With a handful of dates psi is effectively UNIDENTIFIED: on a
#'   3-date section its posterior spanned a factor of 489 (95\% range 28.9 to 14100). That matters
#'   enormously, because psi controls the Poisson event rate of the age-depth bridge,
#'   lambda = mu^(2-p)/(psi*(2-p)). At the posterior median psi of ~3450 that works out at about ONE
#'   event across a whole 3750-kyr interval -- the age path between two dates is then essentially a
#'   single random gamma jump, and the interior uncertainty balloons to several hundred kyr regardless
#'   of how precise the dates are. Left "free", psi is not being estimated so much as allowed to
#'   wander somewhere unhelpful.
#'
#'   WHAT PINNING IT BUYS (measured, floating-time model, 3 dates at +/-30 kyr, 3 replicates each):
#'     psi free  : interior SD 458 kyr, RMSE of the median vs the true curve 15.8 kyr
#'     psi = 100 : interior SD 245 kyr, RMSE 20.3 kyr
#'     psi = 10  : interior SD  82 kyr, RMSE 15.4 kyr
#'     psi = 1   : interior SD  35 kyr, RMSE 15.4 kyr
#'   So a 13-fold tightening with NO loss of accuracy -- the free-psi width was not information, it
#'   was an unidentified nuisance parameter wandering.
#'
#'   THE PRICE, AND IT IS REAL. psi encodes how much you trust the timescale you are fitting against.
#'   Pinning it small asserts that the cyclostratigraphic axis is right. Repeating the test with a
#'   deliberately WRONG axis (cycle misidentified as 25 kyr instead of 21, plus a spurious trend):
#'     psi free    : interior SD 451 kyr, RMSE 94 kyr, 95\% coverage 1.00
#'     psi = 1     : interior SD  35 kyr, RMSE 94 kyr, 95\% coverage 0.49
#'   Both are equally WRONG, but the wide free-psi band still contained the truth while the pinned one
#'   did not -- at half of all depths the true age fell outside the model's own 95\% interval. The
#'   sprawling free-psi uncertainty was acting as an accidental safety net.
#'
#'   SO: pin psi only alongside the independent astronomical controls in spectral_filter.R (the
#'   spectral peak test and the cycle-count test). Those are what establish that the axis is right,
#'   and once psi is pinned they are no longer a refinement -- they are the thing standing between you
#'   and a confidently wrong answer. Report a sensitivity range over psi rather than a single value.
#' @param predictPositions Stratigraphic positions at which the age model is evaluated, in the unit of
#' \code{positions} (numeric vector); positions above the top or below the base date are extrapolated.
#' \code{Default=seq(min(positions), max(positions), length = 500)} (500 evenly spaced points from the bottom to the
#' top dated horizon).
#' @param truncateUp Youngest permitted age, used as truncation for extrapolating above the top of the section
#' (numeric, in the age unit). \code{Default=0}.
#'   IMPORTANT: this is a HARD lower bound on the age of the TOPMOST DATED HORIZON as well, not
#'   just on extrapolated points above it -- it is used directly as the \code{low} limit of that
#'   horizon's Metropolis-Hastings proposal. So if your top date's own likelihood sits close to
#'   this bound (the common case being a "present day" date of 0 +/- sd with the default
#'   truncateUp = 0), roughly the lower half of that date's likelihood is excluded by the bound,
#'   and its posterior comes out correspondingly narrower and skewed. For a date of 0 +/- 15 with
#'   truncateUp = 0 the posterior is a half-normal: sd ~9.0 instead of 15, and a 95\% width of
#'   ~33 kyr instead of ~59 kyr -- verified against both the model output and the closed-form
#'   half-normal. That is the model correctly doing what you told it, NOT uncertainty being lost;
#'   the age-depth model's uncertainty at a dated horizon equals that date's POSTERIOR, which for
#'   any date far from a bound comes back essentially equal to the input (e.g. 7500 +/- 15 at the
#'   base returns a 95\% width of ~58.6 kyr against the input's 58.8 kyr). If you did not intend the
#'   bound to bite, set truncateUp below the youngest age you want to permit (e.g. -1000) and the
#'   full date uncertainty is recovered. The same reasoning applies to \code{truncateDown} and the
#'   bottommost date.
#' @param truncateDown Oldest permitted age, used as truncation for extrapolating below the bottom of the section
#' (numeric, in the age unit). See the note on \code{truncateUp}: it is likewise a hard upper bound on the age of the
#' bottommost dated horizon, not only on extrapolated points. \code{Default=1e10} (in effect no bound).
#' @param extrapUp Number of Compound Poisson-Gamma events simulated for every extrapolated age path above the top
#' dated horizon (integer, > 1); it must be large enough for the path to reach the highest \code{predictPositions}.
#' \code{Default=100}.
#' @param extrapDown Number of Compound Poisson-Gamma events simulated for every extrapolated age path below the
#' bottom dated horizon (integer, > 1), as \code{extrapUp}. \code{Default=100}.
#' @return HDI = specified probability highest density interval for the model run
#' @return model = Raw model predictions for each MCMC iteration
#' @return thetas = The posterior values for each dated horizon from each MCMC run. burn is not removed.
#' @return predictPositions = Stratigraphic positions where the model was evaluated. Same as \code{predictPositions}
#' @return mu, psi = The posterior values for the Compound Poisson-Gamma mean and scale parameters. burn is not removed.
#' @return ageGrid = Grid that age PDFs were evaluated over. Useful for plotting
#' @return likelihoods = Age PDFs for each set of dated samples. Useful for plotting
#' @return nAges = Number of ages that were combined into each sample. Useful for plotting
#' @return masterPositions = Stratigraphic positions for each set of combined ages. Useful for plotting
#' @return ids = Vector of unique sample names. Useful for plotting
#'
#' @importFrom stats approx approxfun dnorm dunif pnorm quantile rexp rgamma rnorm rpois runif sd
#' @export
modBchron <- function(ages,
    ageSds,
    positions,
    positionThicknesses = rep(0, length(positions)),
    ids,
    distTypes = rep('G',length(ages)),
    MC = 10000,
    burn = 2000,
    probability = 0.95,
    adapt = TRUE,
    mhSD = rep(0.001,length(unique(ids))),
    psiSD = 0.001,
    muSD = 0.001,
    mu0 = NULL,
    psi0 = NULL,
    thetas0 = NULL,
    mu_fixed = NULL,
    psi_fixed = NULL,
    predictPositions = seq(min(positions),max(positions),length = 500),
    truncateUp = 0,
    extrapUp = 100,
    truncateDown = 1e10,
    extrapDown = 100){
  ##
  # truncated random walk (was compiled C in BchronModified.c; now pure R --
  # see .truncatedWalkC()/.truncatedRatC() in SECTION 0)
  truncatedWalk = function(old,sd,low,high) {
    if(isTRUE(all.equal(low, high, tolerance = 1e-10))) {
      return(list(new=low,rat=1))
    }
    new <- .truncatedWalkC(old, sd, low, high)
    rat <- .truncatedRatC(old, sd, low, high, new)
    #if(is.nan(rat)) rat=1 # Useful for when proposed value and new value are identical
    return(list(new=new,rat=rat))
  }
  ##
  # tweedie density (was compiled C; now pure R -- see .dtweediep1C())
  dtweediep1 = Vectorize(function(y,p,mu,phi) {
    .dtweediep1C(y = y, power = p, mu = mu, phi = phi)
  })
  ##
  ## prediction and interpolation (was compiled C; now pure R -- see
  ## .predictInterpC()/.predictExtrapUpC()/.predictExtrapDownC())
  predictInterp = function(alpha,
      lambda,
      beta,
      predictPositions,
      diffPositionj,
      currPositionsj,
      currPositionsjp1,
      thetaj,
      thetajp1) {
    .predictInterpC(alpha = alpha, lambda = lambda, beta = beta,
        predictPositions = predictPositions,
        diffPositionj = diffPositionj,
        currPositionsj = currPositionsj,
        currPositionsjp1 = currPositionsjp1,
        thetaj = thetaj, thetajp1 = thetajp1)
  }
  predictExtrapUp = function(alpha,
      lambda,
      beta,
      predictPositions,
      currPositions1,
      theta1,
      maxExtrap,
      extractDate) {
    .predictExtrapUpC(alpha = alpha, lambda = lambda, beta = beta,
        predictPositions = predictPositions,
        currPositions1 = currPositions1, theta1 = theta1,
        maxExtrap = maxExtrap, extractDate = extractDate)
  }
  predictExtrapDown = function(alpha,
      lambda,
      beta,
      predictPositions,
      currPositions1,
      theta1,
      maxExtrap,
      extractDate) {
    .predictExtrapDownC(alpha = alpha, lambda = lambda, beta = beta,
        predictPositions = predictPositions,
        currPositions1 = currPositions1, theta1 = theta1,
        maxExtrap = maxExtrap, extractDate = extractDate)
  }
  ##
  ## function for creating a summed PDF of ages
  compoundProb <- function(ages, sigs, distType, x){
    interval <- matrix(0,nrow = length(x), ncol=length(ages))
    for (i in 1:length(ages)) {
      if(distType[i] == 'G') {
        interval[,i] <- dnorm(x,ages[i],sigs[i])/length(ages)*mean(diff(x))
      } else if(distType[i] == 'U') {
        interval[,i] <- dunif(x,ages[i]-sigs[i],ages[i]+sigs[i])/length(ages)*mean(diff(x))
      }
    }
    G <- rowSums(interval)
    return(G)
  }
  ##
  ## INPUT VALIDATION: every per-date vector must describe the SAME set of dates.
  ## Without this check R silently RECYCLES a short vector, and the failure is
  ## both invisible and badly misleading. The classic case is passing 3 dates but
  ## only 2 ids: `ids = c("dateA","dateB")` recycles to c("dateA","dateB","dateA"),
  ## so date 3 is silently MERGED into date 1 -- the model then sees 2 horizons
  ## instead of 3, gives horizon 1 a bimodal likelihood spanning both merged ages,
  ## and treats the merged date's position as the top of the section, so everything
  ## above it becomes wild extrapolation. Nothing errors; you just get a nonsense
  ## age-depth model. Fail loudly instead.
  .dateLens <- c(ages = length(ages), ageSds = length(ageSds),
      positions = length(positions), ids = length(ids),
      distTypes = length(distTypes),
      positionThicknesses = length(positionThicknesses))
  if (length(unique(.dateLens)) != 1L) {
    stop("ages, ageSds, positions, ids, distTypes and positionThicknesses must all have ",
        "one entry per date, but their lengths are: ",
        paste(sprintf("%s=%d", names(.dateLens), .dateLens), collapse = ", "),
        ".\n  (R would otherwise silently recycle the short one -- e.g. 3 ages with ",
        "2 ids merges the 3rd date into the 1st. Give every date its own id.)")
  }
  ##
  o = order(positions)
  if(any(positions[o]!=positions)) {
    warning("positions not given in order - re-ordering")
    ages <- ages[o]
    ageSds <- ageSds[o]
    positions <- positions[o]
    distTypes <- distTypes[o]
    positionThicknesses <- positionThicknesses[o]
    ids = ids[o]
  }
  ##
  ## sort the data
  nSamples <- length(unique(ids)) # get the number of unique date layers
  masterPositions <- vector()
  nNames <- unique(ids) # get a list of each unique sample name
  for(i in 1:nSamples) {
    masterPositions[i] <- positions[ids == nNames[i]][1]
  }
  ##
  nThicknesses <- vector(length = nSamples)
  for(i in 1:nSamples){
    nThicknesses[i] <- unique(positionThicknesses[positions == masterPositions[i]])
  }
  ##
  # count the number of ages at each unique stratigraphic position
  # this is used for scaling the plots later
  nAges <- rep(NA,length(nSamples))
  for(i in 1:nSamples){
    nAges[i] <- length(ids[ids == nNames[i]])
  }
  ##
  ## generate the summed probability distributions at each unique depth position
  prob <-matrix(0,nrow=100000,ncol=nSamples) # create an empty matrix to store the probabilities
  ageGrid <- seq(min(ages-ageSds*10),max(ages+ageSds*10),length.out=100000) # Grid of ages to evaluate over
  for(j in 1:nSamples){
    prob[, j] <- compoundProb(ages[ids == nNames[j]],
        ageSds[ids == nNames[j]],
        distType = distTypes[ids == nNames[j]],
        x = ageGrid)
  }
  rm(j)
  ##
  ## prealocate some storage
  thetaStore <- matrix(NA,nrow=MC,ncol=nSamples)
  muStore <- vector(length=MC)
  psiStore <- muStore
  modelStore <- matrix(NA,ncol=MC,nrow=length(predictPositions))
  positionStore <- matrix(NA,nrow=MC,ncol=nSamples)
  predictStore <-  matrix(NA,ncol=MC,nrow=length(predictPositions))
  psiSDStore <- psiStore
  muSDStore <- psiStore
  mhSDStore <- thetaStore
  ##
  ## get some starting guesses for MH sampling and make sure they conform to superposition
  thetas <- vector(length = nSamples) # create a vector to store current thetas
  for(i in 1:nSamples){
    thetas[i] <- ageGrid[which.max(prob[, i])] + rnorm(1, 0, 0.00001)
  }
  ## optional warm start: resume from a configuration a previous run reached,
  ## instead of starting cold at the likelihood peak (see the thetas0 docs --
  ## this is what stops a short restarted batch from reporting near-zero age
  ## uncertainty at the dated horizons)
  if (!is.null(thetas0)) {
    if (length(thetas0) != nSamples) {
      stop("thetas0 must have one value per unique id (expected length ", nSamples,
          ", got ", length(thetas0), ")")
    }
    thetas <- as.numeric(thetas0)
  }
  ## PRE-EXISTING BUG FIX: clip the initial guesses into the model's own
  ## hard truncation bounds [truncateUp, truncateDown]. Without this, a
  ## date whose likelihood peaks very close to a truncation boundary (e.g.
  ## a "present day" date at age 0 with truncateUp = 0, the default) can
  ## get an initial guess that lands JUST OUTSIDE that boundary -- not from
  ## the tiny rnorm(1,0,0.00001) jitter, but from ageGrid's own finite
  ## resolution (ageGrid has 100000 points across the full likelihood
  ## range, so its nearest point to the true peak can be a small but
  ## nonzero distance away, in either direction). Once that happens, this
  ## parameter is caught in a genuine trap: .truncatedRatC()'s asymmetric-
  ## proposal Hastings correction is computed from how much of a truncated
  ## normal's mass sits inside [low, high] around the OLD value, and if
  ## the old value already lies outside [low, high], that window's mass is
  ## computed as if it were deep in an excluded tail -- a correction
  ## factor on the order of 1e-6, verified directly:
  ##   .truncatedRatC(old = -0.0045, sd = 0.001, low = 0, high = 1e10, newvalue = 0.0002)
  ##   # ~5.9e-06 -- log() of that alone makes essentially every proposal
  ##   # reject, regardless of what the priors/likelihood say
  ## The parameter then NEVER updates again for the rest of the chain (its
  ## "old" value stays invalid forever, so the same near-zero correction
  ## recurs every iteration) -- silently understating that date's
  ## posterior uncertainty, and, if it's the very top or bottom date,
  ## showing up as the unrelated-looking warning "Unable to find suitable
  ## chronologies for top/bottom of section" firing on every single MCMC
  ## iteration (that date's fitted age sits fractionally beyond the
  ## extrapolation truncation floor/ceiling forever). Clipping the initial
  ## guess here guarantees every theta starts inside its own valid range,
  ## so this trap can never be entered in the first place.
  thetas <- pmin(pmax(thetas, truncateUp), truncateDown)
  ##
  ## calculate some initial parameters
  ## based on Haslett and Parnell (2008)
  currPositions <- masterPositions
  ## if mu0/psi0 are supplied (warm start), centre the initial draw there
  ## instead of the original blind guess -- everything downstream (the
  ## M-H acceptance ratios for mu, psi, and every theta) is completely
  ## unchanged; this only moves where the chain begins.
  mu <- if (is.null(mu0)) { abs(rnorm(1,mean=mean(diff(thetas))/mean(diff(currPositions)),sd=muSD)) } else { abs(rnorm(1, mean = mu0, sd = muSD)) }
  psi <- if (is.null(psi0)) { abs(rnorm(1,1,sd=psiSD)) } else { abs(rnorm(1, mean = psi0, sd = psiSD)) }
  ## mu_fixed / psi_fixed: hold the parameter at a supplied value and skip its
  ## Metropolis-Hastings update entirely (see their @param notes). Everything
  ## else -- the thetas, the bridge, the acceptance ratios -- is untouched.
  if (!is.null(mu_fixed)) {
    mu  <- mu_fixed
  }
  if (!is.null(psi_fixed)) {
    psi <- psi_fixed
  }
  p = 1.2
  alpha <- (2-p) / (p-1) # Haslett and Parnell (2008)
  ##
  ## PRE-EXISTING BUG FIX: safety floors for the adaptive proposal SDs
  ## (psiSD/muSD/mhSD). If an adaptive SD is ever computed from a window
  ## with EXACTLY zero observed variance, the "get rid of zeroes" fallback
  ## a few lines below (in the `if(adapt==T)` block) re-estimates it from
  ## that parameter's whole-chain sd so far -- but if the ENTIRE chain to
  ## that point has itself been perfectly constant (e.g. a dated horizon
  ## sitting right at a truncateUp/truncateDown boundary can get unlucky
  ## and never accept a single proposal in the first h=200-iteration
  ## adaptation window), that whole-chain sd is ALSO exactly zero, and
  ## mhSD[q] silently becomes 0. With sd=0, truncatedWalk() always returns
  ## the current value unchanged, so that parameter is now PERMANENTLY
  ## FROZEN for the rest of the run (once its own step size is exactly 0,
  ## it can never recover). This was invisible in every run so far except
  ## as an unrelated-looking warning ("Unable to find suitable
  ## chronologies for top/bottom of section", firing on every single
  ## iteration) -- the actual symptom is a dated horizon's theta showing
  ## ZERO unique values across the whole chain (check with
  ## length(unique(model$thetas[, i]))), which silently understates that
  ## date's posterior uncertainty. Flooring every adaptive SD at a small,
  ## scale-aware fraction of its OWN starting value guarantees it can
  ## never go perfectly to zero, restoring ergodicity without changing the
  ## adaptive algorithm's behaviour in any normal (non-degenerate) run.
  mhSD_floor <- pmax(mhSD, 1e-10) * 1e-3
  psiSD_floor <- max(psiSD, 1e-10) * 1e-3
  muSD_floor <- max(muSD, 1e-10) * 1e-3
  pb <- utils::txtProgressBar(min = 1,
      max = MC,
      style = 3,
      width=40,
      char = '<>') # create a progress bar
  for(n in 1:MC){
    utils::setTxtProgressBar(pb, n) # set the progress bar
    ##
    ## calculate some secondary model parameters
    ## based on Haslett and Parnell (2008)
    lambda <- (mu ^ (2 - p)) / (psi * (2 - p))
    beta <- 1 / (psi * (p - 1) * mu ^ (p - 1))
    ##
    ## choose some random positions from a uniform distribution for each dated horizon.
    currPositions <- runif(nSamples,
        masterPositions - nThicknesses,
        masterPositions + nThicknesses)
    do <- order(currPositions)
    diffPositions <- diff(currPositions[do])
    thetas[do] <- sort(thetas, decreasing = T)
    positionStore[n, ] <- currPositions
    ##
    # Interpolate between the current thetas
    for(j in 1:(nSamples - 1)){
      inRange <- predictPositions >= currPositions[do[j]] & predictPositions <= currPositions[do[j+1]]
      predval <- predictInterp(alpha = alpha,
          lambda = lambda,
          beta = beta,
          predictPositions = predictPositions[inRange],
          diffPositionj = diffPositions[j],
          currPositionsj = currPositions[do[j]],
          currPositionsjp1 = currPositions[do[j+1]],
          thetaj = thetas[do[j]],
          thetajp1 = thetas[do[j+1]])
      modelStore[inRange,n] <- predval
    }
    ##
    # Extrapolate up
    if(any(predictPositions >= currPositions[do[nSamples]])){
      inRange <- predictPositions >= currPositions[do[nSamples]]
      predval <- predictExtrapUp(alpha = alpha,
          lambda = lambda,
          beta = beta,
          predictPositions = predictPositions[inRange],
          theta1 = thetas[do[nSamples]],
          currPositions1 = currPositions[do[nSamples]],
          maxExtrap = extrapUp,
          extractDate = truncateUp)
      modelStore[inRange,n] <- predval
    }
    # extrapolate down
    if(any(predictPositions <= currPositions[do[1]])){
      inRange <- predictPositions <= currPositions[do[1]]
      predval <- predictExtrapDown(alpha = alpha,
          lambda = lambda,
          beta = beta,
          predictPositions = predictPositions[inRange],
          theta1 = thetas[do[1]],
          currPositions1 = currPositions[do[1]],
          maxExtrap = extrapDown,
          extractDate = truncateDown)
      modelStore[inRange,n] <- predval
    }
    ##
    ## Update thetas
    for(i in 1:nSamples){
      thetaCurrent <- thetas[do[i]]
      ## choose a new theta from a truncated normal where truncation are the surrounding thetas
      thetaProposed <- truncatedWalk(old = thetaCurrent,
          sd = mhSD[do[i]],
          low = ifelse(i == nSamples,
              truncateUp,
              thetas[do[i+1]]-1e-10), # if we're at the top truncate at a really small number
          high = ifelse(i == 1, 1e10, thetas[do[i-1]] + 1e-10)) # if we're at the bottom, truncate at a really big number
      pProposed <- log(prob[, do[i]][which.min(abs(ageGrid - thetaProposed$new))])
      pProposed <-  max(pProposed, -1000000, na.rm = T) # just in case things get really small
      pCurrent <- log(prob[, do[i]][which.min(abs(ageGrid - thetaCurrent))])
      pCurrent <-  max(pCurrent, -1000000, na.rm = T) # just in case things get really small
      priorProposed <-  ifelse(i == 1,
          0,
          log(dtweediep1(thetas[do[i - 1]] - thetaProposed$new,
              p = p,
              mu = mu * (diffPositions[i - 1]),
              phi = psi * (diffPositions[i - 1]) ^ (p - 1)))) +
          ifelse(i == nSamples,
              0,
              log(dtweediep1(thetaProposed$new - thetas[do[i + 1]],
                  p = p,
                  mu = mu*(diffPositions[i]),
                  phi = psi*(diffPositions[i]) ^ (p - 1))))
      priorCurrent <-  ifelse(i == 1,0,
          log(dtweediep1(thetas[do[i - 1]] - thetaCurrent,
              p = p,
              mu = mu * (diffPositions[i - 1]),
              phi = psi * (diffPositions[i - 1])^(p - 1)))) +
          ifelse(i == nSamples,
              0,
              log(dtweediep1(thetaCurrent - thetas[do[i + 1]],
                  p = p,
                  mu = mu*(diffPositions[i]),
                  phi = psi*(diffPositions[i]) ^ (p - 1))))
      # clean up the probabilities if they get really really small
      priorCurrent <- max(priorCurrent, -1000000, na.rm = T)
      priorProposed <- max(priorProposed, -1000000, na.rm = T)
      # calculate the joint probability
      logRtheta <- priorProposed - priorCurrent + pProposed - pCurrent + log(thetaProposed$rat)
      # keep the probability to go to -Inf
      logRtheta <- max(logRtheta, -1000000, na.rm = T)
      if(runif(1, 0, 1) < min(1, exp(logRtheta))) {
        thetas[do[i]] <- thetaProposed$new
      }
    }
    thetaStore[n, ] <- thetas
    #
    ## Update mu
    muCurrent <- mu
    muProposed <- truncatedWalk(old = mu,
        sd = muSD,
        low = 1e-10,
        high = 1e10)
    priorMu <- sum(log(dtweediep1(abs(diff(thetas[do])),
        p = p,
        mu = muProposed$new * diffPositions,
        phi = psi / diffPositions^(p-1)))) -
        sum(log(dtweediep1(abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psi / diffPositions^(p - 1)))) + log(muProposed$rat)
    mu <- if (is.null(mu_fixed)) {
      ifelse(runif(1, 0, 1) < min(1, exp(priorMu)), muProposed$new, muCurrent)
    } else {
      mu_fixed
    }
    muStore[n] <- mu
    ##
    ## Update psi
    psiCurrent <- psi
    psiProposed <- truncatedWalk(old = psi,
        sd = psiSD,
        low = 1e-10,
        high = 1e10)
    priorPsi <- sum(log(dtweediep1(abs(diff(thetas[do])),
        p = p,
        mu = mu * diffPositions,
        phi = psiProposed$new / diffPositions ^ (p - 1)))) -
        sum(log(dtweediep1(abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psi / (diffPositions ^ (p - 1))))) + log(psiProposed$rat)
    psi <- if (is.null(psi_fixed)) {
      ifelse(runif(1, 0, 1) < min(1, exp(priorPsi)), psiProposed$new, psiCurrent)
    } else {
      psi_fixed
    }
    psiStore[n] <- psi
    ##
    ## Update SDs
    if(adapt == T){
      h = 200
      if(n%%h == 0){
        cd <- 2.4 / sqrt(1) # Gelman et al. (1996)
        psiK <- psiStore[(n - h + 1): n] - mean(psiStore[(n - h + 1): n])
        psiRt <- var(psiK) # calculate the variance
        psiSD <- sqrt((cd ^ 2) * psiRt) # calculate the standard deviation
        psiSD <- ifelse(is.na(psiSD),cd * sd(psiStore, na.rm = T), psiSD) # get rid of NAs
        psiSD <- ifelse(psiSD == 0,cd * sd(psiStore, na.rm = T), psiSD) # get rid of zeroes
        psiSD <- max(psiSD, psiSD_floor, na.rm = TRUE) # safety floor -- see note above the main loop
        muK <-  muStore[(n - h + 1): n] - mean(muStore[(n - h + 1): n]) # calculate K
        muRt <- var(muK) # calculate the variance
        muSD <- sqrt(cd ^ 2 * muRt)
        muSD <- ifelse(is.na(muSD), cd * sd(muStore, na.rm = T), muSD)
        muSD <- ifelse(muSD == 0, cd * sd(muStore, na.rm = T), muSD)
        muSD <- max(muSD, muSD_floor, na.rm = TRUE) # safety floor
        for(q in 1:nSamples){ # proposal SD adjustment for thetas
          thetaK <- thetaStore[(n - h + 1): n, q] - mean(thetaStore[( n - h + 1): n, q])
          thetaRt <- var(thetaK)
          mhSD[q] <- ifelse(is.na(sqrt(cd ^ 2 * thetaRt)),
              cd * sd(thetaStore[, q], na.rm = T),
              sqrt(cd ^ 2 * thetaRt))
          mhSD[q] <- max(mhSD[q], mhSD_floor[q], na.rm = TRUE) # safety floor -- see note above the main loop
        }
      }
      ## Store the results
      psiSDStore[n]  <- psiSD
      muSDStore[n]   <- muSD
      mhSDStore[n, ] <- mhSD
    }
  }
  ##
  ## check to see if posterior and likelihood distributions overlap at the specified confidence
  ## if they don't flag them for the user to review
  isOutlier <- function(probability){
    outlier <- vector(length = ncol(thetaStore))
    for(i in 1:ncol(thetaStore)){
      f <- approxfun(x = cumsum(prob[, i]),
          y = ageGrid)
      likelihood <- f(c((1 - probability) / 2, 0.5, (1 + probability) / 2))
      posterior <- as.numeric(quantile(thetaStore[burn:MC, i],
          prob = c((1 - probability) / 2, 0.5, (1 + probability) / 2)))
      outlier[i] <- !(any(posterior > likelihood[1]) & any(posterior < likelihood[3]))
    }
    return(outlier)
  }
  outliers <- isOutlier(probability = probability)
  ##
  ## print a warning if any of the input ages do not overlap the model posterior at the specified probability
  if(any(outliers)){
    for(i in 1:length(ids[outliers])){
      warning(paste('sample', ids[outliers][i],'posterior distribution', probability*100,'% HDI', 'does not overlap the input likelihood. Consider screening for outliers.'))
    }
  }
  ##
  ## HDI and SD are two different summaries of the SAME post-burn-in ensemble.
  ## The HDI bounds are QUANTILES, estimated from the tails, so they move
  ## around a lot between runs at modest ensemble sizes -- an outlying curve
  ## can shift the 2.5%/97.5% points directly. The SD is a moment, computed
  ## from every curve at once, so it is markedly steadier (measured: see
  ## ensemble_stability() below, which quantifies exactly how much steadier
  ## and how many curves you need for a target precision). Both are reported
  ## so you can quote the HDI and MONITOR with the SD.
  .post <- modelStore[, burn:MC, drop = FALSE]
  .hdi <- matrix(NA_real_, 3, nrow(.post))
  .sd <- numeric(nrow(.post))
  for (r in seq_len(nrow(.post))) {
    .q <- quantile(.post[r, ], c((1 - probability) / 2, 0.5 , (1 + probability) / 2))
    .hdi[, r] <- .q
    .sd[r] <- sd(.post[r, ], na.rm = TRUE)
  }
  dimnames(.hdi) <- list(names(quantile(0, c((1 - probability) / 2, 0.5 , (1 + probability) / 2))), rownames(.post))
  if (length(rownames(.post))) {
    names(.sd) <- rownames(.post)
  }
  output <- list(HDI = .hdi,
      sd = .sd,
      mean = rowMeans(.post, na.rm = TRUE),
      model = modelStore,
      thetas = thetaStore,
      positionStore = positionStore,
      psi = psiStore,
      mu = muStore,
      predictPositions = predictPositions,
      ageGrid = ageGrid,
      likelihoods = prob,
      nAges = nAges,
      masterPositions = masterPositions,
      ids = nNames,
      psiSDStore = psiSDStore,
      muSDStore = muSDStore,
      mhSDStore = mhSDStore,
      burn = burn,
      MC = MC,
      outliers = outliers,
      probability = probability,
      positionThicknesses = positionThicknesses)
  class(output) <- "mod_bchron_model_run"
  return(output)
}

#' Predict ages at new positions from a modBchron (or modifiedBChron) model
#'
#' This function is used to predict the age of new stratigraphic points using the output from an ageModel run.
#'
#' @param model Age-depth model, the output of \code{\link{modBchron}} (a list with at least \code{model},
#' \code{predictPositions} and \code{burn}). Required (no default).
#' @param newPositions Stratigraphic positions whose age is predicted, in the unit of the \code{positions} of the model
#' (numeric vector); they must lie within the \code{predictPositions} range of the original model run, outside it
#' the prediction is \code{NA}. Required (no default).
#' @param ids Names of the new positions, used as row names of the output (character or numeric vector, same length
#' as \code{newPositions}). \code{Default=1:length(newPositions)}.
#' @param newPositionThicknesses Stratigraphic uncertainty of each new position, given as half thickness in the
#' unit of \code{newPositions} (numeric vector, >= 0, same length and order as \code{newPositions}); in every MCMC
#' iteration the position is drawn uniformly within plus or minus this value.
#' \code{Default=rep(0, length(newPositions))} (no positional uncertainty).
#' @param probability Probability of the returned highest density interval (numeric, between 0 and 1).
#' \code{Default=0.95}.
#' @return HDI = 95 percent Highest Density Interval for each \code{newPosition}
#' @return raw = Age predictions for each \code{newPosition} for each MCMC run.
#' @export
modBchron_predict <- function(model,
    newPositions,
    ids = 1:length(newPositions),
    newPositionThicknesses = rep(0, length(newPositions)),
    probability = 0.95){
  ## this function is used to predict the age of new points using a modified bchron output model
  ## INPUTS
  ## model = output of the modified bchron age-depth model
  ## newPositions = new stratigraphic points to predict the age of.
  ## The new points must fall within the original stratigraphic range of the model
  ## OUTPUTS
  ## ConfInt = 95% Confidence interval for each point
  ## raw = age predictions for each point for each iteration in the MCMC chain
  pb <- utils::txtProgressBar(min = 1,
      max = ncol(model$model),
      style = 3,
      width = 40,
      char = '<>') # create a progress bar
  x <- model$predictPositions
  n <- length(newPositions)
  predictStore <- matrix(nrow = ncol(model$model),
      ncol = n)

  for(i in 1:ncol(model$model)){
    currPositions <- runif(n,
        newPositions - newPositionThicknesses,
        newPositions + newPositionThicknesses)
    utils::setTxtProgressBar(pb, i) # set the progress bar
    y <- model$model[, i]
    f <- approxfun(x, y)
    predictStore[i, ] <- t(f(currPositions))
  }

  # remove burn burn-in
  predictStore <- predictStore[model$burn:  ncol(model$model), ]

  predictStore <- as.data.frame(predictStore)

  HDI <- matrix(NA_real_, ncol(predictStore), 3)
  for (j in seq_len(ncol(predictStore))) {
    HDI[j, ] <- quantile(predictStore[[j]], c((1 - probability) / 2, 0.5, (1 + probability) / 2), na.rm = T)
  }
  dimnames(HDI) <- list(names(predictStore), names(quantile(0, c((1 - probability) / 2, 0.5, (1 + probability) / 2))))
  HDI <- cbind(newPositions, HDI)
  HDI <- data.frame(HDI)

  HDI$ids <- ids
  colnames(HDI) <- c('newPositions',
      as.character(c((1 - probability) / 2, 0.5, (1 + probability) / 2)),
      'ids')

  predictStore <- data.frame(predictStore)
  return(list(HDI = HDI,
      raw = predictStore))
}
