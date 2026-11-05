#' Joint Bayesian age-depth model of radioisotopic dates and a floating astrochronology (astroMCBchron)
#'
#' @description
#' \code{astroMCBchron()} combines radioisotopic dates (e.g. CA-ID-TIMS U-Pb) with a
#' floating astrochronology (the tracked period of an astronomical cycle, typically the
#' 405-kyr eccentricity cycle, as obtained with \code{\link{track_period_wavelet}},
#' \code{\link{completed_series}} or \code{\link{retrack_wt_MC}}) in one joint posterior
#' distribution. The floating astrochronology is not treated as exact: between each pair of
#' dates its duration is multiplied by a stretch factor whose prior width is set by the
#' tracking uncertainty (the spread between independently tracked proxies), and the duration
#' of every hiatus is estimated together with the age model. The dates enter the likelihood
#' only, so the astrochronology and the dates remain independent sources of information.
#'
#' The model (Arts et al., in prep.):
#' \itemize{
#'   \item the record is divided into intervals bounded by the dated horizons and the top and
#'   base of the record. In interval \emph{j} the floating duration (number of tracked cycles
#'   x \code{target_period}) is multiplied by a stretch factor \emph{s_j}, with
#'   log \emph{s_j} ~ N(0, tau_j^2 + tau0^2). tau_j is the standard error of the floating
#'   duration of the interval (between-proxy spread / sqrt(n proxies)) and tau0 a common
#'   additional slack with a half-normal prior (sd = \code{tau0_sd}). The two outer
#'   intervals, which are not bracketed by dates, are tied to their neighbours
#'   (sd = \code{outer_sd} in log units);
#'   \item each hiatus has a duration H with a normal prior (\code{hiatuses$mean},
#'   \code{hiatuses$sd}), truncated at zero;
#'   \item the stratigraphic position of each date may vary uniformly within
#'   +/- \code{thickness} (give the real half-thickness of the sampled interval, e.g. 0.005 m
#'   for a bed located to the nearest cm); dates that share a \code{depth_group} can in
#'   addition be shifted together by a common depth offset (prior sd \code{depth_offset_sd}),
#'   for dates whose depths come from a different depth scale than the tracks (e.g. outcrop
#'   samples correlated into a core, or box-label depths);
#'   \item the floating astrochronology is drawn from the realisations of the tracked period
#'   (one realisation per segment and MCMC state);
#'   \item the age of every position follows from the age of an anchor date plus the
#'   cumulative stretched floating durations and the hiatuses;
#'   \item date model (\code{date_model}): by default (\code{"mixture"}) every date has a mixture
#'   likelihood, (1 - p) x N(date, sd) + p x outlier component, with the outlier probability
#'   p estimated (beta prior \code{outlier_prior}) and the posterior outlier probability of
#'   each date reported. The outlier component is a Gaussian with \code{outlier_k} x sd
#'   (\code{"scaled"}), or one-sided: \code{"old"} (the date may be too old, inherited or recycled
#'   zircon) or \code{"young"} (residual Pb loss). \code{"two_stage"} reproduces version 1 (a robust
#'   Student-t stage flags dates that move more than \code{outlier_z} sd; flagged dates enter
#'   the final stage with \code{outlier_k} x sd); \code{"gaussian"} uses no outlier treatment.
#'   Instead of a Gaussian, any date can have its own likelihood (\code{date_loglik}, e.g.
#'   from single-zircon dates with \code{\link{upb_grain_loglik}}); it is normalised and
#'   combined with the same outlier mixture;
#'   \item systematic uncertainties: dates from different laboratories or tracers
#'   (\code{group}) get a systematic offset relative to the reference group (the group of the
#'   anchor date), with a prior from the tracer calibration uncertainty of both groups, so
#'   that their relative offset is estimated by the astrochronology. The tracer uncertainty
#'   of the reference group and the 238U decay-constant uncertainty are common to all ages
#'   and are reported as additional envelopes (\code{sd2_tracer_kyr}, \code{sd2_total_kyr});
#'   they do not affect durations;
#'   \item the within-interval uncertainty is added as a monotone gamma (Dirichlet) bridge
#'   between the dates (free gamma walks beyond the outermost dates), with a dispersion per
#'   cycle equal to the relative between-proxy spread of the cycle durations; the reported 2
#'   sd is computed analytically (Rao-Blackwellised) from the posterior draws;
#'   \item optionally, every date is predicted from all other dates and the astrochronology
#'   (leave-one-out), which tests the consistency between dates and cycles. The date model
#'   (outlier probability, or the two-stage screening) is refitted in every leave-one-out run,
#'   so that no information from the full fit is reused.
#' }
#' With \code{bchron = TRUE}, a U-Pb-only Bayesian age-depth model with the same dates
#' (\code{\link{modBchron}}, a pure-R implementation of modifiedBChron; Trayler et al.,
#' 2020; Haslett and Parnell, 2008) is run alongside, on the same prediction grid, so that
#' both models and the reduction in uncertainty can be compared directly.
#'
#' @param dates Radioisotopic dates (data frame with one row per dated horizon, at least two rows). Required columns:
#'   \code{age} (in \code{age_unit}), \code{sd} (1 sigma, same unit as \code{age}, > 0) and \code{position} (depth or
#'   height, same unit as the tracks; two dates may not share a position). Optional columns:
#'   \itemize{
#'   \item \code{id} (character; default \code{"date_1"}, \code{"date_2"}, ...);
#'   \item \code{thickness} (half thickness of the uniform stratigraphic position uncertainty, position units;
#'   default 0);
#'   \item \code{mswd} (MSWD of a weighted-mean date; its sd is multiplied by sqrt(MSWD) when MSWD > 1 and
#'   \code{inflate_mswd = TRUE});
#'   \item \code{group} (laboratory or tracer group, character or numeric; dates of other groups than that of the
#'   anchor get an estimated systematic offset; \code{NA} = the first group);
#'   \item \code{tracer_rel} (1 sd relative tracer-calibration uncertainty of the date; overrides the argument
#'   \code{tracer_rel});
#'   \item \code{depth_group} (label of dates whose depths share one unknown offset, e.g. samples correlated from
#'   another section; \code{NA} = none; needs \code{depth_offset_sd}).
#'   }
#'   Give the internal (analytical, X) uncertainty in \code{sd}: the tracer and decay-constant terms are handled by
#'   \code{group}, \code{tracer_rel} and \code{decay_rel}. Required, no default.
#' @param tracks Floating astrochronology: list with one element per continuous segment (number of hiatuses + 1),
#'   ordered from the youngest (top) to the oldest (base) segment. Each element is either (i) a list with
#'   \code{position} (numeric vector) and \code{period} (matrix positions x realisations, or a vector), (ii) a
#'   \code{\link{retrack_wt_MC}}-style object, i.e. a list whose first element is a matrix with the position in its
#'   first column and whose second element is the matrix of period realisations, or (iii) a matrix/data frame with
#'   the position in column 1 and one or more period realisations in the other columns. Periods are in position units
#'   (e.g. m) per \code{target_period}. Without hiatuses a single such object may be given. Required, no default.
#' @param hiatuses Hiatuses (data frame with the columns \code{position}, \code{mean} and \code{sd}: position of the
#'   hiatus and the mean and 1 sigma (> 0) of the normal prior on its duration in kyr, truncated at zero);
#'   \code{NULL} means no hiatus. \code{Default=NULL}.
#' @param position_type Direction of the positions (character): \code{"depth"} (positions increase downwards, i.e.
#'   with age) or \code{"height"} (positions increase upwards). \code{Default=c("depth", "height")}, i.e.
#'   \code{"depth"}.
#' @param age_unit Unit of \code{dates$age}, \code{dates$sd} and of the output ages (character): \code{"Ma"} or
#'   \code{"ka"} (= kyr). \code{Default=c("Ma", "ka")}, i.e. \code{"Ma"}.
#' @param target_period Period in kyr of the tracked astronomical cycle (numeric, > 0); the floating duration of an
#'   interval is the number of tracked cycles times \code{target_period}. \code{Default=405}.
#' @param proxy_tracks Single-proxy tracks used to compute \code{tau} (stretch prior) and \code{eps}
#'   (within-interval spread) from the between-proxy spread (list with one element per segment, each a named list of
#'   per-proxy tracks given as a data frame or matrix with the position in column 1 and the period in column 2);
#'   \code{NULL} uses the spread of the realisations in \code{tracks}. \code{Default=NULL}.
#' @param tau Widths (sd in log units) of the stretch priors that override the computed ones (numeric vector with
#'   one value per interval between consecutive dates, or a single value used for all intervals); \code{NULL}
#'   computes them (see \code{tau_method}). \code{Default=NULL}.
#' @param tau_method How tau is computed from \code{proxy_tracks} (character): \code{"se"} (between-proxy sd /
#'   sqrt(number of proxies)) or \code{"sd"} (between-proxy sd). Without \code{proxy_tracks}, tau is the sd of the
#'   log floating durations over the track realisations. \code{Default=c("se", "sd")}, i.e. \code{"se"}.
#' @param eps Within-interval relative spread per cycle used for the gamma (Dirichlet) bridge between the dates: a
#'   single number (>= 0), or a data frame with the columns \code{position} and \code{eps}; \code{NULL} computes it
#'   from \code{proxy_tracks} (between-proxy standard error per cycle) or, without them, from the spread of the track
#'   realisations. \code{Default=NULL}.
#' @param tau0_sd Standard deviation (log units) of the half-normal prior on the common additional slack tau0 of all
#'   stretch factors (numeric, >= 0); \code{0} fixes tau0 at 0. \code{Default=0.03}.
#' @param outer_sd Standard deviation (log units) of the prior tying the stretch factors of the two outer intervals
#'   (not bracketed by dates) to those of their neighbours (numeric, > 0). \code{Default=0.2}.
#' @param anchor \code{id} of the date used as the age anchor (a single value of the type of \code{dates$id});
#'   \code{NULL} uses the most precise date (smallest \code{sd}). \code{Default=NULL}.
#' @param date_model Likelihood model of the dates (character): \code{"mixture"} (each date is a mixture of a
#'   Gaussian and an outlier component, with an estimated outlier probability; see \code{outlier_type},
#'   \code{outlier_prior}, \code{outlier_k}), \code{"two_stage"} (the robust two-stage screening of version 1: a
#'   Student-t stage flags dates that move more than \code{outlier_z} sd, flagged dates enter the final Gaussian
#'   stage with \code{outlier_k} x sd) or \code{"gaussian"} (plain Gaussian likelihood, no outlier treatment).
#'   \code{Default=c("mixture", "two_stage", "gaussian")}, i.e. \code{"mixture"}.
#' @param two_stage Deprecated, kept for compatibility with version 1 (logical or \code{NULL}): \code{TRUE} is
#'   \code{date_model = "two_stage"}, \code{FALSE} is \code{date_model = "gaussian"}; if given it overrides
#'   \code{date_model}. \code{NULL} uses \code{date_model}. \code{Default=NULL}.
#' @param outlier_type Outlier component of the mixture (character): \code{"scaled"} (Gaussian around the date
#'   with \code{outlier_k} x sd), \code{"old"} (one-sided: the true age may be younger than the date, i.e. the date
#'   too old, e.g. inherited or recycled zircon) or \code{"young"} (one-sided: the date may be too young, e.g.
#'   residual Pb loss). Only used with \code{date_model = "mixture"}.
#'   \code{Default=c("scaled", "old", "young")}, i.e. \code{"scaled"}.
#' @param outlier_prior Shape parameters of the beta prior on the outlier probability p (numeric vector of length 2,
#'   both > 0: \code{c(a, b)}, prior mean a / (a + b)). Only used with \code{date_model = "mixture"}.
#'   \code{Default=c(1, 19)} (prior mean 0.05).
#' @param nu Degrees of freedom of the Student-t likelihood of the dates in stage 1 of
#'   \code{date_model = "two_stage"} (numeric, > 0); smaller values are more robust to outliers. \code{Default=4}.
#' @param outlier_z Threshold for flagging a date in \code{date_model = "two_stage"}: a date is flagged when its
#'   stage-1 shift exceeds \code{outlier_z} times its sd (numeric, > 0). \code{Default=3}.
#' @param outlier_k Width of the outlier component in multiples of the sd of the date (\code{date_model =
#'   "mixture"}), or the factor by which the sd of a flagged date is multiplied in the final stage of
#'   \code{date_model = "two_stage"} (numeric, >= 1). \code{Default=4}.
#' @param inflate_mswd Multiply the sd of dates with an MSWD > 1 (column \code{dates$mswd}) by sqrt(MSWD) (logical,
#'   TRUE or FALSE); without that column it has no effect. \code{Default=TRUE}.
#' @param tracer_rel Relative tracer-calibration uncertainty (1 sd, as a fraction of the age) of the dates (numeric,
#'   >= 0); a column \code{dates$tracer_rel} overrides it per date. It sets the prior of the offsets between
#'   laboratory groups and the tracer envelope \code{sd2_tracer_kyr}. \code{Default=1.5e-4} (0.03 \% at 2 sd,
#'   EARTHTIME ET535/ET2535; Condon et al., 2015).
#' @param decay_rel Relative 238U decay-constant uncertainty (1 sd, as a fraction of the age) (numeric, >= 0); only
#'   used for the total envelope \code{sd2_total_kyr}, it does not affect durations. \code{Default=5.4e-4} (0.108 \% at
#'   2 sd; Jaffey et al., 1971).
#' @param depth_offset_sd Standard deviation (position units) of the common depth offset of each depth group
#'   (numeric, > 0: one value for all groups, or a vector named by the labels in \code{dates$depth_group}). Required
#'   when \code{dates$depth_group} is given, ignored otherwise. \code{Default=NULL}.
#' @param n_iter Number of MCMC iterations of the final stage (integer, > 0; stage 1 of \code{"two_stage"} uses
#'   half). \code{Default=3e5}.
#' @param thin Thinning interval of the MCMC chain: every \code{thin}-th state after the burn-in is kept (integer,
#'   >= 1). \code{Default=50}.
#' @param burn_frac Fraction of the iterations discarded as burn-in, during which the proposal widths are adapted
#'   (numeric, between 0 and 1). \code{Default=0.2}.
#' @param n_tracks Number of track realisations drawn per segment for the floating astrochronology (integer, > 0).
#'   \code{Default=1000}.
#' @param loo Predict every date (or the dates in \code{loo_ids}) from all other dates and the astrochronology
#'   (leave-one-out), which tests the consistency between dates and cycles (logical, TRUE or FALSE); this runs one
#'   extra MCMC per date, with the date model refitted each time. \code{Default=TRUE}.
#' @param loo_iter Number of MCMC iterations per leave-one-out run (integer, > 0); \code{NULL} uses
#'   \code{round(n_iter / 3)}. \code{Default=NULL}.
#' @param predict_positions Positions at which the age model is evaluated (numeric vector, same unit and type as
#'   \code{dates$position}; must span all dated positions); \code{NULL} uses 501 equidistant positions spanning the
#'   dates and the tracks. \code{Default=NULL}.
#' @param max_curves Maximum number of age-depth curves stored in the result (integer, > 0); larger values need more
#'   memory. \code{Default=5000}.
#' @param n_bridge Number of bridge realisations per posterior draw used for the percentiles (integer, > 0);
#'   \code{NULL} uses enough for at least 10000 curves, at most 20 per draw. \code{Default=NULL}.
#' @param date_loglik Own likelihoods of the dates that replace the Gaussian date likelihood, e.g. from
#'   single-zircon dates with \code{\link{upb_grain_loglik}} (list with one function per date, in the order of the
#'   rows of \code{dates}; each takes an age in \code{age_unit} and returns its log-likelihood). The functions are
#'   normalised numerically and combined with the outlier mixture; \code{dates$age} and \code{dates$sd} then only set
#'   the integration grid. Cannot be combined with \code{date_model = "two_stage"} (\code{"mixture"} is used, with a
#'   warning). \code{NULL} uses the Gaussian likelihood. \code{Default=NULL}.
#' @param bchron Also run the U-Pb-only \code{\link{modBchron}} age-depth model with the same dates on the same
#'   prediction grid, for comparison (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param bchron_split With hiatuses, run one \code{\link{modBchron}} model per segment (dates within the segment) and
#'   estimate the hiatus duration from the difference between the models at the hiatus, as in Cramer et al. (2026)
#'   (logical, TRUE or FALSE); \code{FALSE} runs one model through all dates. \code{Default=FALSE}.
#' @param bchron_args Further arguments passed to \code{\link{modBchron}} (named list, e.g. \code{MC}, \code{burn},
#'   \code{truncateUp}, \code{truncateDown}). \code{Default=list()}.
#' @param temper Inverse temperatures for parallel tempering of the final stage (numeric vector, first value 1,
#'   strictly decreasing, all > 0; Swendsen and Wang, 1986), e.g. \code{exp(seq(0, log(0.01), length.out = 8))}.
#'   Only the astronomical part of the posterior (stretch-factor priors, ties of the outer intervals and hiatus
#'   priors) is tempered; the dates, positions and tau0 enter every chain at full weight (partial tempering, Kjellson
#'   et al., 2026). Hot chains can change the number of cycles between dates, so alternative solutions that a single
#'   chain would miss are found and weighted; samples are taken from the chain at 1. Costs one MCMC per temperature.
#'   \code{NULL} runs one chain. \code{Default=NULL}.
#' @param loo_ids \code{id}s of the dates for the leave-one-out test (vector of values of \code{dates$id});
#'   \code{NULL} tests all dates. Only used with \code{loo = TRUE}. \code{Default=NULL}.
#' @param misfit How the posterior departure of an interval's duration from its posterior mean (the misfit between
#'   the floating time scale and the dates) is placed inside the interval (character; one of \code{"window"},
#'   \code{"localised"}, \code{"stretch"}). \code{"window"}: at one random position,
#'   spread uniformly in time over one contiguous window of \code{misfit_window} kyr of floating time (one cycle with a
#'   different accumulation rate); a negative misfit compresses the window (the window is widened where needed so that
#'   no step loses more than 70 \% of its time). \code{"localised"}: at one horizon of unknown position (one step of the
#'   prediction grid; a miscounted cycle or a condensed level). \code{"stretch"}: spread uniformly over the interval
#'   (rigid stretch, the behaviour of earlier versions). With \code{"stretch"} the age between two dates is a fixed
#'   interpolation between two uncertain end ages, whose variance (1 - w)^2 V_i + w^2 V_i+1 + 2 w (1 - w) C is convex and
#'   can be smaller between the dates than at either date. With \code{"localised"} the variance becomes (1 - w) V_i + w
#'   V_i+1 (w = fraction of the interval in floating time) plus the per-cycle bridge variance: linear between the dated
#'   horizons, so it never has an interior minimum. With \code{"window"} the variance is the rigid-stretch variance plus
#'   2 (E[F] - w) Cov(a_i, delta) + (E[F^2] - w^2) E[delta^2] (F = share of the window below w), which tends to the
#'   localised variance for a short window and stays within a few kyr of it for windows of ~100 kyr, while the stored
#'   curves keep time everywhere (no empty steps, no one-step jumps). No option loosens a constraint: the posterior, the
#'   mean age model and the interval durations between dates are the same. Where a hiatus inside an interval (an
#'   uncertainty step at a known position) or Monte Carlo noise would still create a dip, the 2 sd and the half-widths
#'   of the 95 \% band are replaced by their smallest unimodal majorant within each dated interval (conservative;
#'   \code{"window"} and \code{"localised"} only; the amount is reported in \code{summary$guard_raise_kyr}).
#'   \code{Default=c("window", "localised", "stretch")}, i.e. \code{"window"}.
#' @param misfit_window Length in kyr of floating time of the window used for \code{misfit = "window"} (numeric,
#'   single value > 0). \code{Default=100} (about one short-eccentricity cycle).
#' @param band How the 95 \% band and the median in \code{summary} are computed (character): \code{"analytic"} mean
#'   +/- 1.96 x the analytic sd, median = mean (no Monte Carlo noise), or \code{"curves"} the 2.5, 50 and 97.5
#'   percentiles of the age-depth curves. Both are returned (\code{median_mc}, \code{lower_mc}, \code{upper_mc} hold
#'   the percentiles). \code{Default=c("analytic", "curves")}, i.e. \code{"analytic"}.
#' @param seed Seed of the random number generator (integer), so that the MCMC results can be reproduced;
#'   \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param verbose Print progress messages (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param genplot Plot the result with \code{plot()} (\code{plot.astroMCBchron}) (logical, TRUE or FALSE).
#'   \code{Default=FALSE}.
#'
#' @return An object of class \code{astroMCBchron}: a list with
#' \describe{
#'   \item{model}{matrix of age-depth curves (rows = \code{predictPositions}, columns =
#'   realisations), in \code{age_unit}; the same layout as the \code{model} element of
#'   \code{\link{modBchron}} output.}
#'   \item{predictPositions}{positions of the rows of \code{model}.}
#'   \item{summary}{data.frame: position, mean, median, lower and upper end of the 95 \% band (see \code{band}), the
#'   internal 2 sd (\code{sd2_kyr}), the 2 sd including the tracer term (\code{sd2_tracer_kyr}) and including the
#'   decay constant (\code{sd2_total_kyr}) with the corresponding 95 \% ranges (\code{lower_tracer},
#'   \code{upper_tracer}, \code{lower_total}, \code{upper_total}), the 2.5, 50 and 97.5 percentiles of the curves
#'   (\code{median_mc}, \code{lower_mc}, \code{upper_mc}) and the amount by which the guard raised the 2 sd
#'   (\code{guard_raise_kyr}).}
#'   \item{dates}{data.frame: date, effective date (mean and 2 sd of its likelihood), posterior
#'   age and 2 sd, shift (kyr and sd), posterior outlier probability and flag (> 0.5),
#'   leave-one-out prediction and z-score.}
#'   \item{outlier_prob, group_offsets, depth_offsets}{posterior draws of the outlier
#'   probability (\code{date_model = "mixture"} only), the laboratory offsets (kyr) and the depth offsets.}
#'   \item{stretch}{matrix of posterior stretch factors (draws x intervals).}
#'   \item{hiatus}{matrix of posterior hiatus durations (kyr).}
#'   \item{tau, tau0, eps}{stretch-prior widths, posterior of the common slack and the
#'   within-interval spread profile.}
#'   \item{stage1, loo}{stage-1 screening and leave-one-out tables.}
#'   \item{bchron, comparison}{the \code{modBchron} model(s) and a comparison table of the
#'   2 sd of both models (if \code{bchron = TRUE}).}
#'   \item{posterior, tempering, acceptance, call, settings}{raw posterior states, parallel-tempering diagnostics
#'   (swap acceptance, round trips; \code{NULL} without \code{temper}), acceptance rates and settings.}
#' }
#'
#' @references
#' Condon, D.J., Schoene, B., McLean, N.M., Bowring, S.A., Parrish, R.R., 2015. Metrology and traceability of U-Pb
#' isotope dilution geochronology (EARTHTIME Tracer Calibration Part I). Geochimica et Cosmochimica Acta 164,
#' 464-480.
#'
#' Cramer, B.D., Schmitz, M.D., Trayler, R.B., Stolfus, B.M., Peate, D.W., McAdams,
#' N.E.B., Arts, M., Da Silva, A.-C., 2026. Quantifying uncertainty in the structure and
#' duration of Silurian biogeochemical events. Scientific Reports.
#'
#' Haslett, J., Parnell, A., 2008. A simple monotone process with application to
#' radiocarbon-dated depth chronologies. Journal of the Royal Statistical Society C 57,
#' 399-418.
#'
#' Jaffey, A.H., Flynn, K.F., Glendenin, L.E., Bentley, W.C., Essling, A.M., 1971. Precision measurement of
#' half-lives and specific activities of 235U and 238U. Physical Review C 4, 1889-1906.
#'
#' Swendsen, R.H., Wang, J.-S., 1986. Replica Monte Carlo simulation of spin-glasses. Physical Review Letters 57,
#' 2607-2609.
#'
#' Trayler, R.B., Schmitz, M.D., Cuitino, J.I., Kohn, M.J., Bargo, M.S., Kay, R.F.,
#' Stromberg, C.A.E., Vizcaino, S.F., 2020. An improved approach to age-modeling in deep
#' time: Implications for the Santa Cruz Formation, Argentina. GSA Bulletin 132, 233-244.
#'
#' @examples
#' \donttest{
#' ## synthetic example: accumulation rate varying between 1 and 3 cm/kyr (so that interpolating
#' ## between the dates is wrong), a 300-kyr hiatus at 60 m, five dates with 1 sigma = 80 kyr and a
#' ## floating 405-kyr astrochronology (the tracked period in m, 40 realisations with 2 % noise)
#' set.seed(1)
#' rate <- function(z) 0.02 * (1 + 0.5 * sin(2 * pi * z / 40))                  # m/kyr
#' zz <- seq(0, 100, 0.05)
#' age_z <- 400000 + c(0, cumsum(diff(zz) / rate(zz[-1]))) + ifelse(zz > 60, 300, 0)  # kyr
#' true_age <- function(z) approx(zz, age_z, z)$y
#' z_dates <- c(10, 35, 55, 70, 95)
#' dates <- data.frame(id = paste0("ash", 1:5), position = z_dates,
#'                     age = true_age(z_dates) / 1000 + rnorm(5, 0, 0.08), sd = 0.08,
#'                     thickness = 0.5)
#' trk <- function(z) {
#'   period <- matrix(NA_real_, length(z), 40)
#'   for (i in 1:40) {
#'     period[, i] <- 405 * rate(z) * exp(rnorm(1, 0, 0.02))
#'   }
#'   list(position = z, period = period)
#' }
#' fit <- astroMCBchron(dates, tracks = list(trk(seq(0, 60, 0.25)), trk(seq(60, 100, 0.25))),
#'                      hiatuses = data.frame(position = 60, mean = 280, sd = 60),
#'                      n_iter = 5000, loo = FALSE, bchron = TRUE,
#'                      bchron_args = list(MC = 1000, burn = 300))
#' fit
#' plot(fit)
#' # error against the truth: astrochronology + dates against the dates alone
#' s <- fit$summary
#' mean(abs(s$median - true_age(s$position) / 1000)) * 1000                     # kyr
#' }
#' @importFrom stats approx dnorm dt median quantile rgamma rnorm runif sd var setNames
#' @export
astroMCBchron <- function(dates, tracks, hiatuses = NULL, position_type = c("depth", "height"),
    age_unit = c("Ma", "ka"), target_period = 405,
    proxy_tracks = NULL, tau = NULL, tau_method = c("se", "sd"), eps = NULL,
    tau0_sd = 0.03, outer_sd = 0.2, anchor = NULL,
    date_model = c("mixture", "two_stage", "gaussian"), two_stage = NULL,
    outlier_type = c("scaled", "old", "young"), outlier_prior = c(1, 19),
    nu = 4, outlier_z = 3, outlier_k = 4, inflate_mswd = TRUE,
    tracer_rel = 1.5e-4, decay_rel = 5.4e-4, depth_offset_sd = NULL,
    n_iter = 3e5, thin = 50, burn_frac = 0.2, n_tracks = 1000,
    loo = TRUE, loo_iter = NULL, predict_positions = NULL,
    max_curves = 5000, n_bridge = NULL, date_loglik = NULL,
    bchron = FALSE, bchron_split = FALSE, bchron_args = list(),
    temper = NULL, loo_ids = NULL,
    misfit = c("window", "localised", "stretch"), misfit_window = 100, band = c("analytic", "curves"),
    seed = NULL, verbose = TRUE, genplot = FALSE) {
  cl <- match.call()
  misfit <- match.arg(misfit)
  band <- match.arg(band)
  if (!is.numeric(misfit_window) || length(misfit_window) != 1 || !is.finite(misfit_window) || misfit_window <= 0) {
    stop("misfit_window must be one positive number (kyr)")
  }
  position_type <- match.arg(position_type)
  age_unit <- match.arg(age_unit)
  tau_method <- match.arg(tau_method)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sc <- if (age_unit == "Ma") { 1000 } else { 1 }                                 # input unit -> kyr
  sgn <- if (position_type == "depth") { 1 } else { -1 }                           # internal z: increases with age
  msg <- function(...) if (verbose) { cat(sprintf(...)) }
  # Dates
  dates <- as.data.frame(dates)
  for (nm in c("age", "sd", "position")) {
    if (is.null(dates[[nm]])) {
      stop("dates needs a column '", nm, "'")
    }
  }
  if (is.null(dates$id)) {
    dates$id <- paste0("date_", seq_len(nrow(dates)))
  }
  if (is.null(dates$thickness)) {
    dates$thickness <- 0
  }
  if (nrow(dates) < 2) {
    stop("at least two dates are needed")
  }
  if (any(!is.finite(dates$age) | !is.finite(dates$sd) | dates$sd <= 0)) {
    stop("dates: ages and sds must be finite and sds > 0")
  }
  ord <- order(sgn * dates$position)
  dates <- dates[ord, ]
  rownames(dates) <- NULL
  ## date model: "mixture" (default; outlier mixture with an estimated outlier probability), "two_stage" (screening of
  ## version 1, kept to reproduce earlier results) or "gaussian" (no outlier treatment); two_stage = TRUE/FALSE is accepted
  ## for backward compatibility
  date_model <- if (!is.null(two_stage)) {
    if (isTRUE(two_stage)) {
      "two_stage"
    } else {
      "gaussian"
    }
  } else {
    match.arg(date_model)
  }
  outlier_type <- match.arg(outlier_type)
  if (!is.null(date_loglik)) {
    date_loglik <- date_loglik[ord]
    if (date_model == "two_stage") {
      warning("date_loglik cannot be combined with the two-stage screening; date_model = 'mixture' is used")
      date_model <- "mixture"
    }
  }
  if (length(outlier_prior) != 2 || any(outlier_prior <= 0)) {
    stop("outlier_prior must be the two shape parameters of a beta prior")
  }
  nb <- nrow(dates)
  zb <- sgn * dates$position
  th <- abs(dates$thickness)
  ## over-dispersed weighted means: inflate the uncertainty by sqrt(MSWD) when MSWD > 1
  if (inflate_mswd && !is.null(dates$mswd)) {
    f <- sqrt(pmax(1, ifelse(is.finite(dates$mswd), dates$mswd, 1)))
    if (any(f > 1)) {
      msg("MSWD > 1: uncertainty of %s inflated by %s\n", paste(dates$id[f > 1], collapse = ", "), paste(round(f[f > 1], 2), collapse = ", "))
    }
    dates$sd <- dates$sd * f
  }
  age_k <- dates$age * sc
  sd_k <- dates$sd * sc
  if (any(diff(zb) <= 0)) {
    stop("two dates share the same position; combine them first")
  }
  ia <- if (is.null(anchor)) { which.min(sd_k) } else { match(anchor, dates$id) }
  if (is.na(ia)) {
    stop("anchor id not found in dates$id")
  }
  ## laboratory / tracer groups: one systematic offset per group relative to the reference group (the group of the anchor)
  grp <- if (is.null(dates$group)) { rep("1", nb) } else { as.character(dates$group) }
  grp[is.na(grp)] <- "1"
  glev <- unique(grp)
  gi <- match(grp, glev)
  ng <- length(glev)
  gref <- gi[ia]
  trel <- if (!is.null(dates$tracer_rel)) { dates$tracer_rel } else { rep(tracer_rel, nb) }
  trel[!is.finite(trel)] <- tracer_rel
  sig_tr_g <- numeric(ng)                                                          # 1 sd tracer term per group (kyr)
  for (g in seq_len(ng)) {
    sig_tr_g[g] <- mean(trel[gi == g] * age_k[gi == g])
  }
  sig_doff <- numeric(ng)
  for (g in seq_len(ng)) {
    sig_doff[g] <- if (g == gref) { 0 } else { sqrt(sig_tr_g[g]^2 + sig_tr_g[gref]^2) }
  }
  ## correlated depth offsets: dates that share a depth group are shifted together relative to the depth scale of the tracks
  dg_raw <- if (is.null(dates$depth_group)) { rep(NA_character_, nb) } else { as.character(dates$depth_group) }
  dlev <- unique(dg_raw[!is.na(dg_raw)])
  di <- match(dg_raw, dlev)
  nd <- length(dlev)
  if (nd) {
    if (is.null(depth_offset_sd)) {
      stop("dates$depth_group given: depth_offset_sd (sd of the depth offset per group, position units) is needed")
    }
    osd <- if (!is.null(names(depth_offset_sd))) { depth_offset_sd[dlev] } else { rep(depth_offset_sd, length.out = nd) }
    if (any(!is.finite(osd) | osd <= 0)) {
      stop("depth_offset_sd must be > 0 for every depth group")
    }
  } else {
    osd <- numeric(0)
  }
  # Hiatuses and segments
  if (!is.null(hiatuses) && nrow(hiatuses)) {
    hiatuses <- as.data.frame(hiatuses)
    hiatuses <- hiatuses[order(sgn * hiatuses$position), , drop = FALSE]
    zh <- sgn * hiatuses$position
    Hm <- hiatuses$mean
    Hs <- hiatuses$sd
    if (any(Hs <= 0)) {
      stop("hiatuses$sd must be > 0")
    }
  } else {
    zh <- numeric(0)
    Hm <- Hs <- numeric(0)
  }
  nh <- length(zh)
  nseg <- nh + 1
  single <- function(x) is.matrix(x) || is.data.frame(x) || (is.list(x) && !is.null(x$position)) ||
      (is.list(x) && length(x) == 2 && (is.matrix(x[[1]]) || is.data.frame(x[[1]])) && is.matrix(x[[2]]) && nrow(as.matrix(x[[1]])) == nrow(x[[2]]) && nseg == 1)
  if (single(tracks)) {
    tracks <- list(tracks)
  }
  if (length(tracks) != nseg) {
    stop(sprintf("tracks must have one element per segment: %d hiatus(es) -> %d segments, got %d", nh, nseg, length(tracks)))
  }
  TRK <- vector("list", length(tracks))
  names(TRK) <- names(tracks)
  for (s in seq_along(tracks)) {
    TRK[[s]] <- .amc_track(tracks[[s]])
  }
  for (s in seq_len(nseg)) {
    TRK[[s]]$z <- sgn * TRK[[s]]$position
  }
  # Grid and knots
  trk_z <- vector("list", length(TRK))
  names(trk_z) <- names(TRK)
  n_trk_z <- numeric(length(TRK))
  for (s in seq_along(TRK)) {
    trk_z[[s]] <- TRK[[s]]$z
    n_trk_z[s] <- length(TRK[[s]]$z)
  }
  zall <- c(zb - th, zb + th, unlist(trk_z))
  zpred <- if (is.null(predict_positions)) { seq(min(zall), max(zall), length.out = 501) } else { sort(unique(sgn * predict_positions)) }
  ztop <- min(c(zall, zpred))
  zbot <- max(c(zall, zpred))
  zg <- sort(unique(zpred))
  if (min(zg) > min(zb) || max(zg) < max(zb)) {
    stop("predict_positions must span all dated positions")
  }
  ## memory guard (as in WaverideR 0.6.0): stop before allocating the stored curves and floating-time matrices if they do not fit
  if (exists("wr_ram_available", mode = "function")) {
    est_bytes <- 8 * (length(zg) * max_curves * 3 + sum(n_trk_z) * n_tracks * 3)
    free_bytes <- tryCatch(wr_ram_available(), error = function(e) NA)
    if (!is.na(free_bytes) && est_bytes > 0.7 * free_bytes) {
      stop(sprintf(paste0("astroMCBchron needs about %.1f GB for max_curves = %g, n_tracks = %g and %d prediction positions, ",
          "but only %.1f GB is free. Reduce max_curves, n_tracks or the number of predict_positions."),
          est_bytes / 2^30, max_curves, n_tracks, length(zg), free_bytes / 2^30), call. = FALSE)
    }
  }
  kn <- c(ztop, zb, zbot)
  J <- nb + 1
  seg_of <- function(z) findInterval(z, zh, left.open = TRUE) + 1                  # z exactly at a hiatus belongs to the upper segment
  bounds <- c(ztop, zh, zbot)
  # Floating time per segment (realisations)
  NTR <- n_tracks
  ft_seg <- function(s) {
    t <- TRK[[s]]
    P <- t$period
    k <- t$z >= bounds[s] - 1e-9 & t$z <= bounds[s + 1] + 1e-9
    if (sum(k) < 2) {
      k <- rep(TRUE, length(t$z))
    }                               # segment shorter than the track spacing: use all
    zs <- sort(unique(c(bounds[s], t$z[k], bounds[s + 1])))
    Pz <- .amc_interp(t$z, P, zs)
    dr <- if (ncol(P) >= NTR) { sample.int(ncol(P), NTR) } else { sample.int(ncol(P), NTR, replace = TRUE) }
    Tz <- .amc_cumfloat(zs, Pz[, dr, drop = FALSE], target_period)
    list(z = zs, T = Tz, draws = dr, mean = .amc_cumfloat(zs, Pz, target_period))
  }
  FT <- vector("list", nseg)
  for (s in seq_len(nseg)) {
    FT[[s]] <- ft_seg(s)
  }
  off <- matrix(0, nseg, NTR)
  if (nseg > 1) {
    for (s in 2:nseg) {
      off[s, ] <- off[s - 1, ] + FT[[s - 1]]$T[nrow(FT[[s - 1]]$T), ]
    }
  }
  Fall <- function(z) {
    out <- matrix(0, length(z), NTR)
    s <- seg_of(z)
    for (q in unique(s)) {
      u <- s == q
      out[u, ] <- sweep(.amc_interp(FT[[q]]$z, FT[[q]]$T, z[u]), 2, off[q, ], "+")
    }
    out
  }
  Fk <- function(z, k) {
    out <- numeric(length(z))
    s <- seg_of(z)
    for (q in unique(s)) {
      u <- s == q
      out[u] <- off[q, k] + .amc_interp(FT[[q]]$z, FT[[q]]$T[, k, drop = FALSE], z[u])[, 1]
    }
    out
  }
  Fkn <- Fall(kn)
  D <- matrix(t(diff(Fkn)), nrow = NTR)                          # NTR x J floating durations per interval
  ivl <- function(z) pmin(pmax(findInterval(z, kn, rightmost.closed = TRUE), 1), J)
  loc_bed <- function(i) {
    m <- if (th[i] > 0) { 41 } else { 1 }
    dz <- if (m > 1) { seq(-th[i], th[i], length.out = m) } else { 0 }
    zl <- pmin(pmax(zb[i] + dz, ztop), zbot)
    jl <- ivl(zl)
    list(dz = dz, j = jl, f = t(Fall(zl) - Fkn[jl, , drop = FALSE]), hz = if (nh) { outer(zl, zh, ">") * 1 } else { matrix(0, m, 0) })
  }
  LOC <- vector("list", nb)
  for (i in seq_len(nb)) {
    LOC[[i]] <- loc_bed(i)
  }
  nz_loc <- numeric(nb)                                                           # number of candidate positions within each bed
  for (i in seq_len(nb)) {
    nz_loc[i] <- length(LOC[[i]]$dz)
  }
  # Tau (stretch prior) and eps (within-interval spread)
  PT <- NULL
  if (!is.null(proxy_tracks)) {
    if (length(proxy_tracks) != nseg) {
      stop("proxy_tracks must have one element (named list of proxy tracks) per segment")
    }
    pt_one <- function(x, s) {
      x <- as.matrix(x)
      z <- sgn * x[, 1]
      o <- order(z)
      z <- z[o]
      p <- x[o, 2]
      k <- is.finite(p) & p > 0
      zs <- sort(unique(c(bounds[s], z[k & z >= bounds[s] & z <= bounds[s + 1]], bounds[s + 1])))
      list(z = zs, T = .amc_cumfloat(zs, .amc_interp(z[k], matrix(p[k]), zs), target_period)[, 1])
    }
    PT <- vector("list", nseg)
    for (s in seq_len(nseg)) {
      pxs <- proxy_tracks[[s]]
      pts_s <- vector("list", length(pxs))
      names(pts_s) <- names(pxs)
      for (r in seq_along(pxs)) {
        pts_s[[r]] <- pt_one(pxs[[r]], s)
      }
      PT[s] <- list(pts_s)
    }
  }
  pfloat <- function(nm, z) {
    s <- seg_of(z)
    out <- numeric(length(z))
    tot <- 0
    for (q in seq_len(nseg)) {
      x <- PT[[q]][[nm]]
      if (is.null(x)) {
        if (any(s == q)) {
          return(NA)
        } else {
          next
        }
      }
      u <- s == q
      out[u] <- tot + approx(x$z, x$T, z[u], rule = 2)$y
      tot <- tot + x$T[length(x$T)]
    }
    out
  }
  inner <- if (J > 2) { 2:(J - 1) } else { integer(0) }
  ilab <- sprintf("%s-%s", format(sgn * kn[-length(kn)], nsmall = 2), format(sgn * kn[-1], nsmall = 2))
  if (!is.null(tau)) {
    if (length(tau) == 1) {
      tau <- rep(tau, length(inner))
    }
    if (length(tau) != length(inner)) {
      stop(sprintf("tau must have %d values (one per interval between dates)", length(inner)))
    }
    tau_v <- tau
    tau_src <- "user"
  } else if (!is.null(PT)) {
    tau_proxy <- function(j) {
      sj <- unique(seg_of(c(kn[j] + 1e-9, kn[j + 1])))
      sj <- seq(min(sj), max(sj))
      nms <- names(PT[[sj[1]]])
      for (q in sj[-1]) {
        nms <- intersect(nms, names(PT[[q]]))
      }
      if (length(nms) < 2) {
        return(NA)
      }
      dl <- vector("list", length(nms))
      names(dl) <- nms
      for (r in seq_along(nms)) {
        dl[r] <- list(diff(pfloat(nms[r], c(kn[j], kn[j + 1]))))
      }
      d <- if (all(lengths(dl) == 1L)) { unlist(dl) } else { dl }                  # as sapply() would simplify
      d <- d[is.finite(d) & d > 0]
      if (length(d) < 2) {
        return(NA)
      }
      sd(log(d)) / if (tau_method == "se") { sqrt(length(d)) } else { 1 }
    }
    tau_l <- vector("list", length(inner))
    for (r in seq_along(inner)) {
      tau_l[r] <- list(tau_proxy(inner[r]))
    }
    tau_v <- unlist(tau_l)                                                         # as sapply() would simplify (all length 1)
    if (any(!is.finite(tau_v))) {
      msg("tau: fewer than two proxies in some intervals; using the track-realisation spread there\n")
      ld_in <- log(D[, inner, drop = FALSE])
      sd_tr <- numeric(ncol(ld_in))
      for (r in seq_len(ncol(ld_in))) {
        sd_tr[r] <- sd(ld_in[, r])
      }
      tau_v[!is.finite(tau_v)] <- sd_tr[!is.finite(tau_v)]
    }
    tau_src <- paste0("proxy ", tau_method)
  } else {
    ld_in <- log(pmax(D[, inner, drop = FALSE], 1e-9))
    tau_v <- numeric(ncol(ld_in))
    for (r in seq_len(ncol(ld_in))) {
      tau_v[r] <- sd(ld_in[, r])
    }
    tau_src <- "track realisations"
  }
  names(tau_v) <- ilab[inner]
  ## eps per cycle -> continuous profile on zg
  EPS <- .amc_eps(eps, FT, PT, zg, zh, bounds, target_period, sgn)
  cz <- EPS$profile^2 * target_period
  # Likelihood and prior
  ## Normal component of every date: a Gaussian (age, sd), or the user's date_loglik normalised numerically on a grid.
  ## mu_e and sd_e are the mean and sd of that component; they define the outlier component, the leave-one-out z-score and
  ## the effective date reported in the output.
  if (!is.null(date_loglik)) {
    norm_one <- function(i) {
      g <- age_k[i] + seq(-15, 15, length.out = 3001) * sd_k[i]
      v <- date_loglik[[i]](g / sc)
      v[!is.finite(v)] <- -1e6
      m <- max(v)
      w <- exp(v - m)
      dx <- g[2] - g[1]
      lz <- m + log(sum(w) * dx)
      w <- w / sum(w)
      mu <- sum(w * g)
      list(lz = lz, mu = mu, sd = sqrt(sum(w * (g - mu)^2)))
    }
    LN <- vector("list", nb)
    for (i in seq_len(nb)) {
      LN[[i]] <- norm_one(i)
    }
    mu_e <- numeric(nb)
    sd_e <- numeric(nb)
    lz_e <- numeric(nb)
    for (i in seq_len(nb)) {
      mu_e[i] <- LN[[i]][["mu"]]
      sd_e[i] <- LN[[i]][["sd"]]
      lz_e[i] <- LN[[i]][["lz"]]
    }
    lnorm_comp <- function(i, a) {
      v <- date_loglik[[i]](a / sc)
      if (is.finite(v)) {
        v - lz_e[i]
      } else {
        -1e6
      }
    }
  } else {
    mu_e <- age_k
    sd_e <- sd_k
    lnorm_comp <- function(i, a) dnorm(a, age_k[i], sd_k[i], log = TRUE)
  }
  ## outlier component: "scaled" = Gaussian with outlier_k x sd around the date; "old" = the date may be too old
  ## (inherited or recycled zircon: true age younger than the date, half-normal with scale outlier_k x sd);
  ## "young" = the date may be too young (residual Pb loss)
  lout_comp <- function(i, a) {
    s <- outlier_k * sd_e[i]
    switch(outlier_type, scaled = dnorm(a, mu_e[i], s, log = TRUE),
        old = if (a <= mu_e[i]) { log(2) + dnorm(a, mu_e[i], s, log = TRUE) } else { -Inf },
        young = if (a >= mu_e[i]) { log(2) + dnorm(a, mu_e[i], s, log = TRUE) } else { -Inf })
  }
  LF_mix <- function(i, a, p) {
    ln <- log1p(-p$po) + lnorm_comp(i, a)
    lo <- log(p$po) + lout_comp(i, a)
    m <- max(ln, lo)
    if (!is.finite(m)) {
      return(-1e6)
    }
    m + log(exp(ln - m) + exp(lo - m))
  }
  LF_gauss <- function(i, a, p) lnorm_comp(i, a)
  resp_out <- function(i, a, po) {
    ln <- log1p(-po) + lnorm_comp(i, a)
    lo <- log(po) + lout_comp(i, a)
    if (!is.finite(lo)) {
      0
    } else {
      1 / (1 + exp(ln - lo))
    }
  }
  logL_fun <- function(type, mult) function(i, a, p) {
    x <- (a - mu_e[i]) / (sd_e[i] * mult[i])            # two-stage screening (version 1)
    if (type[i] == "t") {
      dt(x, nu, log = TRUE) - dt(0, nu, log = TRUE)
    } else {
      -0.5 * x^2
    }
  }
  logL_date <- logL_fun(rep("t", nb), rep(1, nb))
  LF_main <- switch(date_model, mixture = LF_mix, gaussian = LF_gauss, two_stage = NULL)
  ## position of a dated bed: nominal position + position within the bed (+/- thickness) + offset of its depth group
  pos_bed <- function(p, i) zb[i] + LOC[[i]]$dz[p$m[i]] + if (!is.na(di[i])) { p$o[di[i]] } else { 0 }
  ages_beds <- function(p) {
    cs <- c(0, cumsum(p$s * D[p$k, ]))
    G <- numeric(nb)
    for (i in seq_len(nb)) {
      L <- LOC[[i]]
      m <- p$m[i]
      if (is.na(di[i])) {
        j <- L$j[m]
        G[i] <- cs[j] + p$s[j] * L$f[p$k, m] + if (nh) { sum(p$H * L$hz[m, ]) } else { 0 }
      } else {
        z <- min(max(pos_bed(p, i), ztop), zbot)
        j <- ivl(z)
        G[i] <- cs[j] + p$s[j] * (Fk(z, p$k) - Fkn[j, p$k]) + if (nh) { sum(p$H * (z > zh)) } else { 0 }
      }
    }
    p$A + G - G[ia]
  }
  dated <- function(p) ages_beds(p) + p$d[gi]                                     # model age in the frame of each laboratory group
  Gz <- function(z, p) {
    F <- Fk(z, p$k)
    Fn <- Fk(kn, p$k)
    j <- ivl(z)
    cs <- c(0, cumsum(p$s * diff(Fn)))
    cs[j] + p$s[j] * (F - Fn[j]) + if (nh) { as.vector(outer(z, zh, ">") %*% p$H) } else { 0 }
  }
  ages_z <- function(z, p) {
    za <- min(max(pos_bed(p, ia), ztop), zbot)
    g <- Gz(c(z, za), p)
    p$A + g[seq_along(z)] - g[length(g)]
  }
  t0max <- max(tau0_sd, 1e-12)
  lprior_h <- function(p) {                                                      # astronomical part (tempered with `temper`)
    if (any(p$s <= 0) || (nh && any(p$H <= 0)) || p$tau0 < 0) {
      return(-Inf)
    }
    ls <- log(p$s)
    lp <- 0
    if (length(inner)) {
      lp <- lp + sum(dnorm(ls[inner], 0, sqrt(tau_v^2 + p$tau0^2 + 1e-12), log = TRUE))
    }
    lp <- lp + dnorm(ls[1] - ls[2], 0, outer_sd, log = TRUE) + dnorm(ls[J] - ls[J - 1], 0, outer_sd, log = TRUE)
    if (nh) {
      lp <- lp + sum(dnorm(p$H, Hm, Hs, log = TRUE))
    }
    lp
  }
  ## non-astronomical priors: tau0, outlier probability (beta prior, sampled on the logit scale: + log Jacobian),
  ## laboratory offsets and depth offsets
  lprior_e <- function(p) {
    lp <- if (tau0_sd > 0) { dnorm(p$tau0, 0, t0max, log = TRUE) + log(max(p$tau0, 1e-12)) } else { 0 }
    if (date_model == "mixture") {
      if (p$po <= 0 || p$po >= 1) {
        return(-Inf)
      }
      lp <- lp + stats::dbeta(p$po, outlier_prior[1], outlier_prior[2], log = TRUE) + log(p$po) + log1p(-p$po)
    }
    if (ng > 1) {
      lp <- lp + sum(dnorm(p$d[-gref], 0, sig_doff[-gref], log = TRUE))
    }
    if (nd) {
      lp <- lp + sum(dnorm(p$o, 0, osd, log = TRUE))
    }
    lp
  }
  lprior <- function(p) {
    h <- lprior_h(p)
    if (!is.finite(h)) {
      return(-Inf)
    }
    h + lprior_e(p)
  }
  llik <- function(p, use, LF) {
    a <- dated(p)
    w <- which(use)
    ll <- numeric(length(w))
    for (r in seq_along(w)) {
      ll[r] <- LF(w[r], a[w[r]], p)
    }
    sum(ll)
  }
  m_init <- integer(nb)                                                           # start in the middle of each bed
  for (i in seq_len(nb)) {
    m_init[i] <- as.integer((length(LOC[[i]]$dz) + 1) / 2)
  }
  p_init <- function() list(A = age_k[ia], s = rep(1, J), tau0 = if (tau0_sd > 0) { tau0_sd } else { 0 }, H = Hm,
      m = m_init, k = 1L,
      po = if (date_model == "mixture") { outlier_prior[1] / sum(outlier_prior) } else { 0 }, d = rep(0, ng), o = rep(0, nd))
  ## proposals for the date-model parameters (outlier probability, laboratory offsets, depth offsets)
  prop_extra <- function(getp, tryf, st) {
    if (date_model == "mixture") {
      p <- getp()
      q <- p
      lg <- stats::qlogis(p$po) + rnorm(1, 0, st$po)
      q$po <- stats::plogis(lg)
      tryf(q, "po")
    }
    if (ng > 1) {
      for (g in setdiff(seq_len(ng), gref)) {
        p <- getp()
        q <- p
        q$d[g] <- p$d[g] + rnorm(1, 0, st$d[g])
        tryf(q, "d")
      }
    }
    if (nd) {
      for (g in seq_len(nd)) {
        p <- getp()
        q <- p
        q$o[g] <- p$o[g] + rnorm(1, 0, st$o[g])
        tryf(q, "o")
      }
    }
  }
  # MCMC (Metropolis-within-Gibbs, adaptive during burn-in)
  run_mcmc <- function(use, niter, LF) {
    p <- p_init()
    lp <- lprior(p) + llik(p, use, LF)
    step <- list(A = min(sd_k), s = 0.05, t = 0.8, H = if (nh) { Hs } else { numeric(0) }, m = pmax(1, round(nz_loc / 6)), po = 0.5, d = pmax(sig_doff, 1e-6) / 2, o = osd / 2)
    nacc <- c(A = 0, s = 0, t = 0, H = 0, z = 0, k = 0, po = 0, d = 0, o = 0)
    win <- nacc
    ntry <- nacc
    wtry <- nacc
    nburn <- floor(burn_frac * niter)
    keep <- vector("list", floor((niter - nburn) / thin))
    nk <- 0
    movable <- which(nz_loc > 1)
    try_q <- function(q, b) {
      ntry[b] <<- ntry[b] + 1
      wtry[b] <<- wtry[b] + 1
      lq <- lprior(q)
      if (!is.finite(lq)) {
        return(FALSE)
      }
      lq <- lq + llik(q, use, LF)
      if (log(runif(1)) < lq - lp) {
        p <<- q
        lp <<- lq
        nacc[b] <<- nacc[b] + 1
        win[b] <<- win[b] + 1
        TRUE
      } else {
        FALSE
      }
    }
    for (it in seq_len(niter)) {
      q <- p
      q$A <- p$A + rnorm(1, 0, step$A)
      try_q(q, "A")
      for (j in seq_len(J)) {
        q <- p
        q$s[j] <- p$s[j] * exp(rnorm(1, 0, step$s))
        try_q(q, "s")
      }
      if (tau0_sd > 0) {
        q <- p
        q$tau0 <- p$tau0 * exp(rnorm(1, 0, step$t))
        try_q(q, "t")
      }
      if (nh) {
        for (h in seq_len(nh)) {
          q <- p
          q$H[h] <- p$H[h] + rnorm(1, 0, step$H[h])
          try_q(q, "H")
        }
      }
      if (length(movable)) {
        i <- movable[sample.int(length(movable), 1)]
        q <- p
        q$m[i] <- p$m[i] + sample(c(-1L, 1L), 1) * sample.int(step$m[i], 1)
        if (q$m[i] >= 1 && q$m[i] <= length(LOC[[i]]$dz)) {
          try_q(q, "z")
        }
      }
      if (NTR > 1) {
        q <- p
        q$k <- sample.int(NTR, 1)
        try_q(q, "k")
      }
      prop_extra(function() p, try_q, step)
      if (it <= nburn && it %% 200 == 0) {                                       # adapt the step sizes towards ~30% acceptance
        r <- win / pmax(wtry, 1)
        step$A <- step$A * ifelse(r["A"] > 0.3, 1.25, 0.8)
        step$s <- step$s * ifelse(r["s"] > 0.3, 1.25, 0.8)
        step$t <- step$t * ifelse(r["t"] > 0.3, 1.25, 0.8)
        if (nh) {
          step$H <- step$H * ifelse(r["H"] > 0.3, 1.25, 0.8)
        }
        step$po <- step$po * ifelse(r["po"] > 0.3, 1.25, 0.8)
        step$d <- step$d * ifelse(r["d"] > 0.3, 1.25, 0.8)
        step$o <- step$o * ifelse(r["o"] > 0.3, 1.25, 0.8)
        win[] <- 0
        wtry[] <- 0
      }
      if (it > nburn && (it - nburn) %% thin == 0) {
        nk <- nk + 1
        keep[[nk]] <- p
      }
    }
    if (nk == 0) {
      keep <- list(p)
      nk <- 1
    }
    keep <- keep[seq_len(nk)]
    attr(keep, "acc") <- round(nacc / pmax(ntry, 1), 2)
    keep
  }
  # Parallel tempering (astronomical part only)
  run_pt <- function(use, niter, LF, betas) {
    S <- length(betas)
    init <- p_init()
    ev <- function(p) {
      h <- lprior_h(p)
      if (!is.finite(h)) {
        return(c(h = -Inf, e = -Inf, l = -Inf))
      }
      c(h = h, e = lprior_e(p), l = llik(p, use, LF))
    }
    P <- rep(list(init), S)
    E <- vector("list", S)
    for (ci in seq_len(S)) {
      E[[ci]] <- ev(P[[ci]])
    }
    st <- vector("list", S)
    for (ci in seq_len(S)) {
      st[[ci]] <- list(A = min(sd_k), s = 0.05, t = 0.8, H = if (nh) { Hs / sqrt(betas[ci]) } else { numeric(0) },
          m = pmax(1, round(nz_loc / 6)), po = 0.5, d = pmax(sig_doff, 1e-6) / 2, o = osd / 2)
    }
    win <- matrix(0, S, 4, dimnames = list(NULL, c("A", "s", "t", "H")))
    wtry <- win
    nburn <- floor(burn_frac * niter)
    keep <- vector("list", floor((niter - nburn) / thin))
    nk <- 0
    hot <- vector("list", floor((niter - nburn) / thin) * S)
    nh2 <- 0
    sw_try <- sw_acc <- numeric(S - 1)
    lab <- seq_len(S)
    dir <- rep(0L, S)
    trips <- 0
    movable <- which(nz_loc > 1)
    for (it in seq_len(niter)) {
      for (c in seq_len(S)) {
        b <- betas[c]
        p <- P[[c]]
        e <- E[[c]]
        tgt <- function(x) x[["e"]] + b * x[["h"]] + x[["l"]]
        tryq <- function(q, blk) {
          if (!is.na(blk)) {
            wtry[c, blk] <<- wtry[c, blk] + 1
          }
          x <- ev(q)
          if (!is.finite(x[["h"]])) {
            return(invisible(FALSE))
          }
          if (log(runif(1)) < tgt(x) - tgt(e)) {
            p <<- q
            e <<- x
            if (!is.na(blk)) {
              win[c, blk] <<- win[c, blk] + 1
            }
          }
        }
        q <- p
        q$A <- p$A + rnorm(1, 0, st[[c]]$A)
        tryq(q, "A")
        for (j in seq_len(J)) {
          q <- p
          q$s[j] <- p$s[j] * exp(rnorm(1, 0, st[[c]]$s))
          tryq(q, "s")
        }
        if (tau0_sd > 0) {
          q <- p
          q$tau0 <- p$tau0 * exp(rnorm(1, 0, st[[c]]$t))
          tryq(q, "t")
        }
        if (nh) {
          for (h in seq_len(nh)) {
            q <- p
            q$H[h] <- p$H[h] + rnorm(1, 0, st[[c]]$H[h])
            tryq(q, "H")
          }
        }
        if (length(movable)) {
          i <- movable[sample.int(length(movable), 1)]
          q <- p
          q$m[i] <- p$m[i] + sample(c(-1L, 1L), 1) * sample.int(st[[c]]$m[i], 1)
          if (q$m[i] >= 1 && q$m[i] <= length(LOC[[i]]$dz)) {
            tryq(q, NA)
          }
        }
        if (NTR > 1) {
          q <- p
          q$k <- sample.int(NTR, 1)
          tryq(q, NA)
        }
        prop_extra(function() p, function(q, b) tryq(q, NA), st[[c]])
        P[[c]] <- p
        E[[c]] <- e
        if (it <= nburn && it %% 200 == 0) {
          r <- win[c, ] / pmax(wtry[c, ], 1)
          f <- ifelse(r > 0.3, 1.25, 0.8)
          st[[c]]$A <- st[[c]]$A * f[["A"]]
          st[[c]]$s <- st[[c]]$s * f[["s"]]
          st[[c]]$t <- st[[c]]$t * f[["t"]]
          if (nh) {
            st[[c]]$H <- st[[c]]$H * f[["H"]]
          }
          win[c, ] <- 0
          wtry[c, ] <- 0
        }
      }
      for (c in seq(1 + it %% 2, S - 1, by = 2)) {
        if (S > 1 && c < S) {          # non-reversible even/odd swaps (Syed et al., 2022)
          sw_try[c] <- sw_try[c] + 1
          la <- (betas[c] - betas[c + 1]) * (E[[c + 1]][["h"]] - E[[c]][["h"]])
          if (log(runif(1)) < la) {
            tmp <- P[[c]]
            P[[c]] <- P[[c + 1]]
            P[[c + 1]] <- tmp
            tmp <- E[[c]]
            E[[c]] <- E[[c + 1]]
            E[[c + 1]] <- tmp
            tl <- lab[c]
            lab[c] <- lab[c + 1]
            lab[c + 1] <- tl
            sw_acc[c] <- sw_acc[c] + 1
          }
        }
      }
      ## round trips of the replica that sits in the cold chain: cold -> hottest -> cold
      for (c in c(1, S)) {
        r <- lab[c]
        if (c == S && dir[r] == 1L) {
          dir[r] <- 2L
        }
        if (c == 1 && dir[r] == 2L) {
          trips <- trips + 1
          dir[r] <- 1L
        }
        if (c == 1 && dir[r] == 0L) {
          dir[r] <- 1L
        }
      }
      if (it > nburn && (it - nburn) %% thin == 0) {
        nk <- nk + 1
        keep[[nk]] <- P[[1]]
        for (c in seq_len(S)) {
          nh2 <- nh2 + 1
          hot[[nh2]] <- c(beta = betas[c], logpost = sum(E[[c]]), H = if (nh) { P[[c]]$H } else { NULL }, s = P[[c]]$s, A = P[[c]]$A)
        }
      }
    }
    if (nk == 0) {
      keep <- list(P[[1]])
      nk <- 1
    }
    keep <- keep[seq_len(nk)]
    attr(keep, "acc") <- round(colSums(win) / pmax(colSums(wtry), 1), 2)
    attr(keep, "pt") <- list(betas = betas, swap_acceptance = round(sw_acc / pmax(sw_try, 1), 3), round_trips = trips,
        states = as.data.frame(do.call(rbind, hot[seq_len(nh2)])))
    keep
  }
  # Sampling: mixture / gaussian, or the two-stage screening of version 1
  stage1 <- NULL
  outlier <- rep(FALSE, nb)
  t0 <- Sys.time()
  msg("astroMCBchron: %d dates, %d hiatus(es), %d intervals, %d track realisations; tau from %s; date model %s%s%s\n", nb, nh, J, NTR, tau_src, date_model,
      if (date_model == "mixture") { sprintf(" (%s outliers, beta(%g, %g) prior on the outlier probability)", outlier_type, outlier_prior[1], outlier_prior[2]) } else { "" },
      if (ng > 1) { sprintf("; %d laboratory groups", ng) } else { "" })
  screen <- function(use, niter) {                                                # stage 1 of the two-stage screening (version 1)
    p1 <- run_mcmc(use, niter, logL_date)
    P1 <- matrix(NA_real_, length(p1), nb)
    for (d in seq_along(p1)) {
      P1[d, ] <- dated(p1[[d]])
    }
    sh <- colMeans(P1) - mu_e
    list(P1 = P1, shift = sh, z = sh / sd_e, flag = abs(sh / sd_e) > outlier_z & use)
  }
  if (date_model == "two_stage") {
    s1 <- screen(rep(TRUE, nb), round(n_iter / 2))
    outlier <- s1$flag
    stage1 <- data.frame(id = dates$id, position = dates$position, date = dates$age, date_sd = dates$sd, stage1 = colMeans(s1$P1) / sc,
        shift_kyr = round(s1$shift, 1), shift_sd = round(s1$z, 2), outlier = outlier)
    msg("stage 1 (Student-t, nu = %g): largest shift %.0f kyr (%.1f sd); %d outlier(s)%s\n", nu, s1$shift[which.max(abs(s1$z))], max(abs(s1$z)),
        sum(outlier), if (any(outlier)) { paste0(": ", paste(dates$id[outlier], collapse = ", ")) } else { "" })
    LF_main <- logL_fun(rep("n", nb), ifelse(outlier, outlier_k, 1))
  }
  if (!is.null(temper)) {
    if (abs(temper[1] - 1) > 1e-12 || any(diff(temper) >= 0)) {
      stop("temper must start at 1 and decrease")
    }
    post <- run_pt(rep(TRUE, nb), n_iter, LF_main, temper)
    msg("parallel tempering (%d chains): swap acceptance %s; round trips %d\n", length(temper), paste(attr(post, "pt")$swap_acceptance, collapse = " "), attr(post, "pt")$round_trips)
  } else {
    post <- run_mcmc(rep(TRUE, nb), n_iter, LF_main)
  }
  msg("MCMC: %d iterations, %d draws, %.1f min; acceptance %s\n", n_iter, length(post), as.numeric(difftime(Sys.time(), t0, units = "mins")),
      paste(names(attr(post, "acc")), attr(post, "acc"), collapse = " "))
  # Posterior summaries
  npost <- length(post)
  PA <- matrix(NA_real_, npost, nb)                                               # model age of each dated bed in the frame of its laboratory
  S <- matrix(NA_real_, npost, J)
  HH <- if (nh) { matrix(NA_real_, npost, nh) } else { matrix(numeric(0), npost, 0) }
  T0 <- numeric(npost)
  PO <- numeric(npost)
  GO <- if (ng > 1) { matrix(NA_real_, npost, ng) } else { NULL }
  DO <- if (nd) { matrix(NA_real_, npost, nd) } else { NULL }
  for (d in seq_len(npost)) {
    p <- post[[d]]
    PA[d, ] <- dated(p)
    S[d, ] <- p$s
    if (nh) {
      HH[d, ] <- p$H
    }
    T0[d] <- p[["tau0"]]
    PO[d] <- p[["po"]]
    if (ng > 1) {
      GO[d, ] <- p$d
    }
    if (nd) {
      DO[d, ] <- p$o
    }
  }
  colnames(S) <- ilab
  if (nh == 1) {
    HH <- matrix(HH, ncol = 1)
  }
  if (nh) {
    colnames(HH) <- sprintf("hiatus_%g", sgn * zh)
  }
  if (ng > 1) {
    colnames(GO) <- glev
  }
  if (nd) {
    colnames(DO) <- dlev
  }
  p_out <- if (date_model == "mixture") {
    R_out <- matrix(NA_real_, npost, nb)
    for (r in seq_len(npost)) {
      for (i in seq_len(nb)) {
        R_out[r, i] <- resp_out(i, PA[r, i], PO[r])
      }
    }
    colMeans(R_out)
  } else {
    rep(NA_real_, nb)
  }
  if (date_model == "mixture") {
    outlier <- p_out > 0.5
  }
  if (date_model == "mixture") {
    msg("outlier probability: %.3f [%.3f-%.3f]; posterior outlier probability per date: max %.2f (%s)\n", median(PO), quantile(PO, .025), quantile(PO, .975),
        max(p_out), dates$id[which.max(p_out)])
  }
  # Age-depth ensemble with the gamma bridge
  hstep <- integer(nh)
  for (h in seq_len(nh)) {
    hstep[h] <- which(zg > zh[h])[1]
  }
  kb <- findInterval(zb[if (any(!outlier)) { !outlier } else { TRUE }], zg, rightmost.closed = TRUE)
  incr <- function(p) {
    a <- ages_z(zg, p)
    inc <- c(0, pmax(diff(a), 1e-9))
    if (nh) {
      for (h in seq_len(nh)) {
        if (!is.na(hstep[h])) {
          inc[hstep[h]] <- max(inc[hstep[h]] - p$H[h], 1e-9)
        }
      }
    }
    list(a = a, inc = inc)
  }
  segs <- c(1, kb, length(zg))
  bridge_ages <- function(p, I) {
    inc <- I$inc
    new <- inc
    for (m in seq_len(length(segs) - 1)) {
      if (segs[m + 1] <= segs[m]) {
        next
      }
      ii <- (segs[m] + 1):segs[m + 1]
      g <- rgamma(length(ii), shape = inc[ii] / cz[ii], scale = cz[ii])
      if (m == 1 || m == length(segs) - 1) {
        new[ii] <- g
      } else if (misfit == "stretch") {
        new[ii] <- g / sum(g) * sum(inc[ii])                                    # misfit spread uniformly (rigid stretch)
      } else if (misfit == "window") {
        ## misfit "window": the departure delta is spread uniformly in time over one contiguous window of
        ## misfit_window kyr of floating time at a random position inside the interval (one cycle with a different
        ## accumulation rate). A negative delta compresses the window (rate x (1 + delta / L)); the window is widened
        ## where needed so that no step loses more than 70 % of its time (no empty steps)
        Db <- Dbar[m]
        base <- g / sum(g) * Db
        delta <- sum(inc[ii]) - Db
        L <- min(max(misfit_window, if (delta < 0) { -delta / 0.7 } else { 0 }), Db)
        u <- runif(1, 0, Db - L)
        T0 <- c(0, cumsum(base))
        ov <- pmax(pmin(T0[-1], u + L) - pmax(T0[-length(T0)], u), 0)
        new[ii] <- pmax(base + delta * ov / L, 0)
      } else {
        ## misfit "localised": the floating-model shape (with its per-cycle gamma jitter) carries the posterior-mean
        ## duration Dbar; the draw's departure delta from it sits at ONE horizon of unknown position (uniform in
        ## time) - an extra (delta > 0) or spurious (delta < 0) slice of time, e.g. a miscounted cycle or condensation
        Db <- Dbar[m]
        base <- g / sum(g) * Db
        delta <- sum(inc[ii]) - Db
        if (delta >= 0) {
          j <- sample.int(length(ii), 1, prob = base)
          base[j] <- base[j] + delta
        } else {
          T0 <- c(0, cumsum(base))
          rem <- min(-delta, Db * (1 - 1e-9))
          u <- runif(1, 0, Db)
          cut <- if (u + rem <= Db) { pmin(pmax(T0 - u, 0), rem) } else { pmin(pmax(T0 - u, 0), Db - u) + pmin(T0, u + rem - Db) }
          base <- pmax(diff(T0 - cut), 0)
        }
        new[ii] <- base
      }
    }
    if (nh) {
      for (h in seq_len(nh)) {
        if (!is.na(hstep[h])) {
          new[hstep[h]] <- new[hstep[h]] + p$H[h]
        }
      }
    }
    out <- numeric(length(zg))
    k1 <- kb[1]
    out[k1] <- I$a[k1]
    if (k1 < length(zg)) {
      for (i in (k1 + 1):length(zg)) {
        out[i] <- out[i - 1] + new[i]
      }
    }
    if (k1 > 1) {
      for (i in (k1 - 1):1) {
        out[i] <- out[i + 1] - new[i + 1]
      }
    }
    out
  }
  bridge_var <- function(I) {
    inc <- I$inc
    v <- inc * cz
    out <- numeric(length(zg))
    for (m in seq_len(length(segs) - 1)) {
      if (segs[m + 1] <= segs[m]) {
        next
      }
      ii <- (segs[m] + 1):segs[m + 1]
      if (m == 1) {
        out[segs[m]:(segs[m + 1] - 1)] <- rev(cumsum(rev(v[ii])))
        next
      }
      if (m == length(segs) - 1) {
        out[ii] <- cumsum(v[ii])
        next
      }
      mu <- cumsum(inc[ii])
      Vz <- cumsum(v[ii])
      w <- mu / mu[length(ii)]
      out[ii] <- pmax(Vz - 2 * w * Vz + w^2 * Vz[length(ii)], 0)
    }
    out
  }
  IN <- vector("list", length(post))
  PAR <- matrix(NA_real_, length(zg), length(post))
  BV <- PAR
  for (d in seq_along(post)) {
    IN[[d]] <- incr(post[[d]])
  }
  for (d in seq_along(post)) {
    PAR[, d] <- IN[[d]][["a"]]
    BV[, d] <- bridge_var(IN[[d]])
  }
  var_par <- numeric(length(zg))
  for (r in seq_along(zg)) {
    var_par[r] <- var(PAR[r, ])
  }
  sd_rb <- sqrt(var_par + rowMeans(BV))
  ## posterior-mean floating duration of every interval between dates (hiatuses excluded), used by misfit = "localised"
  ## and "window"
  Dbar <- rep(NA_real_, length(segs) - 1)
  for (m in seq_len(length(segs) - 1)) {
    if (segs[m + 1] <= segs[m]) {
      next
    }
    ii <- (segs[m] + 1):segs[m + 1]
    di <- numeric(length(IN))
    for (d in seq_along(IN)) {
      di[d] <- sum(IN[[d]]$inc[ii])
    }
    Dbar[m] <- mean(di)
  }
  if (misfit == "localised") {
    ## analytic variance with a localised misfit: for an interval with end ages a_i, a_i+1 and duration D = Dbar + delta,
    ## a(phi) = a_i + Dbar * phi + delta * 1(u < phi), u ~ U(0, 1) in floating time, so
    ## Var(phi) = Var(a_i) + phi * Var(D) + 2 * phi * Cov(a_i, D) = (1 - phi) * Var(a_i) + phi * Var(a_i+1)
    ## i.e. the rigid-interpolation variance plus phi * (1 - phi) * Var(D): linear between the dated horizons, never
    ## below the better-dated end and never above the worse-dated end (+ the per-cycle gamma-bridge variance)
    IM <- matrix(NA_real_, length(zg), length(IN))
    for (d in seq_along(IN)) {
      IM[, d] <- IN[[d]]$inc
    }
    inc_mean <- rowMeans(IM)
    for (m in seq_len(length(segs) - 1)) {
      if (m == 1 || m == length(segs) - 1 || segs[m + 1] <= segs[m]) {
        next
      }
      ii <- segs[m]:segs[m + 1]
      phi <- c(0, cumsum(inc_mean[ii[-1]])) / sum(inc_mean[ii[-1]])
      ## general form (also with a hiatus inside the interval): Var = Var_rigid(phi) + phi (1 - phi) E[delta^2];
      ## without a hiatus this is exactly the linear interpolation of the variances at the two dated horizons
      dd <- colSums(IM[ii[-1], , drop = FALSE]) - Dbar[m]
      var_par[ii] <- var_par[ii] + phi * (1 - phi) * mean(dd^2)
    }
    sd_rb <- sqrt(var_par + rowMeans(BV))
  }
  if (misfit == "window") {
    ## analytic variance with a window misfit: a(phi) = a_i + Dbar phi + delta F(phi), F = share of the window below phi,
    ## window [u, u + f] with u ~ U(0, 1 - f), f = misfit_window / Dbar. Relative to the rigid stretch (F = phi):
    ## Var = Var_rigid + 2 (E[F] - phi) Cov(a_i, delta) + (E[F^2] - phi^2) E[delta^2]
    IM <- matrix(NA_real_, length(zg), length(IN))
    for (d in seq_along(IN)) {
      IM[, d] <- IN[[d]]$inc
    }
    inc_mean <- rowMeans(IM)
    ug <- (seq_len(2000) - 0.5) / 2000
    for (m in seq_len(length(segs) - 1)) {
      if (m == 1 || m == length(segs) - 1 || segs[m + 1] <= segs[m]) {
        next
      }
      ii <- segs[m]:segs[m + 1]
      phi <- c(0, cumsum(inc_mean[ii[-1]])) / sum(inc_mean[ii[-1]])
      f <- min(misfit_window / Dbar[m], 1)
      u <- ug * (1 - f)
      EF <- numeric(length(phi))
      EF2 <- numeric(length(phi))
      for (r in seq_along(phi)) {
        Fr <- pmin(pmax((phi[r] - u) / f, 0), 1)
        EF[r] <- mean(Fr)
        EF2[r] <- mean(Fr^2)
      }
      dd <- colSums(IM[ii[-1], , drop = FALSE]) - Dbar[m]
      ai <- PAR[segs[m], ]
      cv <- mean((ai - mean(ai)) * dd)
      var_par[ii] <- var_par[ii] + 2 * (EF - phi) * cv + (EF2 - phi^2) * mean(dd^2)
    }
    sd_rb <- sqrt(pmax(var_par, 0) + rowMeans(BV))
  }
  ndr <- length(post)
  NBQ <- if (is.null(n_bridge)) { max(1L, min(20L, as.integer(ceiling(10000 / ndr)))) } else { as.integer(n_bridge) }
  NBS <- max(1L, min(NBQ, as.integer(floor(max_curves / ndr))))
  AGQ <- matrix(NA_real_, length(zg), ndr * NBQ)
  for (d in seq_len(ndr)) {
    for (b in seq_len(NBQ)) {
      AGQ[, (d - 1) * NBQ + b] <- bridge_ages(post[[d]], IN[[d]])
    }
  }
  qq <- matrix(NA_real_, 3, nrow(AGQ), dimnames = list(names(quantile(0, c(.025, .5, .975))), NULL))
  for (r in seq_len(nrow(AGQ))) {
    qq[, r] <- quantile(AGQ[r, ], c(.025, .5, .975))
  }
  guard_raise <- rep(0, length(zg))
  if (misfit %in% c("localised", "window")) {
    sd_pre <- sd_rb
    ## guard: between two consecutive dated horizons the uncertainty may never have a local minimum. The localised
    ## misfit makes the variance linear there; a hiatus inside an interval (uncertainty step at a known position) or
    ## Monte Carlo noise in the quantiles could still create a dip, so the 2 sd and the half-widths of the 95 % band
    ## are replaced by their smallest unimodal majorant within every dated interval (conservative: never reduced)
    qc <- function(x) pmin(cummax(x), rev(cummax(rev(x))))
    for (m in seq_len(length(segs) - 1)) {
      if (m == 1 || m == length(segs) - 1 || segs[m + 1] <= segs[m]) {
        next
      }
      ii <- segs[m]:segs[m + 1]
      sd_rb[ii] <- qc(sd_rb[ii])
      qq[1, ii] <- qq[2, ii] - qc(qq[2, ii] - qq[1, ii])
      qq[3, ii] <- qq[2, ii] + qc(qq[3, ii] - qq[2, ii])
    }
    guard_raise <- sd_rb - sd_pre
  }
  qq_curves <- qq
  if (band == "analytic") {
    ## 95 % band from the analytic mean and sd (no Monte Carlo noise); the curve percentiles are kept as lower_mc/upper_mc
    mu_a <- rowMeans(PAR)
    qq[1, ] <- mu_a - 1.959964 * sd_rb
    qq[2, ] <- mu_a
    qq[3, ] <- mu_a + 1.959964 * sd_rb
  }
  AG <- AGQ[, as.vector(outer(seq_len(NBS), (seq_len(ndr) - 1) * NBQ, `+`)), drop = FALSE]
  rm(AGQ)
  ## systematic uncertainty, added in quadrature (a common offset of all dates of the reference laboratory; it does not
  ## affect durations): tracer calibration of the reference group, and the 238U decay constant (common to all U-Pb dates)
  mage <- rowMeans(PAR)
  s_tr <- mean(trel[gi == gref]) * mage
  s_dc <- decay_rel * mage
  wid <- function(dq, s) sqrt(dq^2 + (2 * s)^2)
  summ <- data.frame(position = sgn * zg, mean = mage / sc, median = qq[2, ] / sc, lower = qq[1, ] / sc, upper = qq[3, ] / sc,
      sd2_kyr = 2 * sd_rb, sd2_tracer_kyr = 2 * sqrt(sd_rb^2 + s_tr^2), sd2_total_kyr = 2 * sqrt(sd_rb^2 + s_tr^2 + s_dc^2),
      lower_tracer = (qq[2, ] - wid(qq[2, ] - qq[1, ], s_tr)) / sc, upper_tracer = (qq[2, ] + wid(qq[3, ] - qq[2, ], s_tr)) / sc,
      lower_total = (qq[2, ] - wid(qq[2, ] - qq[1, ], sqrt(s_tr^2 + s_dc^2))) / sc, upper_total = (qq[2, ] + wid(qq[3, ] - qq[2, ], sqrt(s_tr^2 + s_dc^2))) / sc,
      median_mc = qq_curves[2, ] / sc, lower_mc = qq_curves[1, ] / sc, upper_mc = qq_curves[3, ] / sc, guard_raise_kyr = 2 * guard_raise)
  if (misfit %in% c("localised", "window")) {
    ## the same guard for the envelopes with systematic uncertainty (added in quadrature, slowly varying with age)
    for (m in seq_len(length(segs) - 1)) {
      if (m == 1 || m == length(segs) - 1 || segs[m + 1] <= segs[m]) {
        next
      }
      ii <- segs[m]:segs[m + 1]
      summ$sd2_tracer_kyr[ii] <- qc(summ$sd2_tracer_kyr[ii])
      summ$sd2_total_kyr[ii] <- qc(summ$sd2_total_kyr[ii])
      for (nm in c("tracer", "total")) {
        lo <- paste0("lower_", nm)
        up <- paste0("upper_", nm)
        summ[[lo]][ii] <- summ$median[ii] - qc(summ$median[ii] - summ[[lo]][ii])
        summ[[up]][ii] <- summ$median[ii] + qc(summ[[up]][ii] - summ$median[ii])
      }
    }
  }
  # Leave-one-out (each date predicted without its own age;
  ## the date model is refitted without that date, so no information from the full fit is reused)
  L <- NULL
  if (loo) {
    li <- if (is.null(loo_iter)) { round(n_iter / 3) } else { loo_iter }
    lid <- if (is.null(loo_ids)) { seq_len(nb) } else { match(loo_ids, dates$id) }
    L <- matrix(NA_real_, nb, 3, dimnames = list(NULL, c("pred", "sd", "z")))
    for (i in lid) {
      use <- seq_len(nb) != i
      LFi <- if (date_model == "two_stage") {
        f1 <- screen(use, round(li / 2))$flag
        logL_fun(rep("n", nb), ifelse(f1, outlier_k, 1))
      } else {
        LF_main
      }
      pp <- if (is.null(temper)) { run_mcmc(use, li, LFi) } else { run_pt(use, li, LFi, temper) }
      a <- numeric(length(pp))
      for (d in seq_along(pp)) {
        a[d] <- dated(pp[[d]])[i]
      }
      L[i, ] <- c(mean(a), sd(a), (mean(a) - mu_e[i]) / sqrt(sd(a)^2 + sd_e[i]^2))
    }
    msg("leave-one-out: %d of %d dates predicted within 2 sd\n", sum(abs(L[, "z"]) <= 2, na.rm = TRUE), sum(is.finite(L[, "z"])))
  }
  sd_pa <- numeric(nb)
  for (i in seq_len(nb)) {
    sd_pa[i] <- sd(PA[, i])
  }
  dtab <- data.frame(id = dates$id, position = dates$position, group = grp, date = dates$age, date_2sd_kyr = 2 * sd_k,
      date_eff = mu_e / sc, date_eff_2sd_kyr = 2 * sd_e,
      posterior = colMeans(PA) / sc, posterior_2sd_kyr = 2 * sd_pa,
      shift_kyr = colMeans(PA) - mu_e, shift_sd = (colMeans(PA) - mu_e) / sd_e, p_outlier = p_out, outlier = outlier)
  if (!is.null(L)) {
    dtab$loo_pred <- L[, "pred"] / sc
    dtab$loo_2sd_kyr <- 2 * L[, "sd"]
    dtab$loo_z <- L[, "z"]
  }
  # U-Pb-only comparison model
  BC <- NULL
  cmp <- NULL
  if (bchron) {
    BC <- .amc_bchron(dates, zg, zh, sgn, sc, bchron_split, bchron_args, verbose)
    cmp <- data.frame(position = sgn * zg, astro_2sd_kyr = 2 * sd_rb, bchron_2sd_kyr = BC$sd2_kyr)
    cmp$reduction_pct <- 100 * (1 - cmp$astro_2sd_kyr / cmp$bchron_2sd_kyr)
    inr <- zg >= min(zb) & zg <= max(zb)
    msg("uncertainty reduction vs modBchron between the outermost dates: mean %.0f%% (median %.0f%%)\n", mean(cmp$reduction_pct[inr], na.rm = TRUE),
        median(cmp$reduction_pct[inr], na.rm = TRUE))
  }
  out <- list(model = AG / sc, predictPositions = sgn * zg, summary = summ, dates = dtab,
      stretch = S, hiatus = HH, tau = tau_v, tau0 = T0, outlier_prob = if (date_model == "mixture") { PO } else { NULL },
      group_offsets = GO, depth_offsets = DO, eps = EPS, stage1 = stage1,
      loo = if (is.null(L)) { NULL } else { dtab[, c("id", "position", "date", "date_eff", "loo_pred", "loo_2sd_kyr", "loo_z")] },
      bchron = BC, comparison = cmp, posterior = post, tempering = attr(post, "pt"),
      acceptance = attr(post, "acc"), call = cl,
      settings = list(misfit = misfit, misfit_window = misfit_window, band = band, position_type = position_type,
          age_unit = age_unit, target_period = target_period, hiatuses = hiatuses,
          anchor = dates$id[ia], tau_source = tau_src, n_iter = n_iter, thin = thin, n_tracks = NTR,
          bridge_per_draw = NBQ, curves_stored = ncol(AG), date_model = date_model, outlier_type = outlier_type,
          outlier_prior = outlier_prior, outlier_k = outlier_k, groups = glev, reference_group = glev[gref],
          tracer_rel_ref = mean(trel[gi == gref]), decay_rel = decay_rel, depth_groups = dlev, depth_offset_sd = osd,
          date_loglik = !is.null(date_loglik)))
  class(out) <- "astroMCBchron"
  if (genplot) {
    plot(out)
  }
  out
}

# Internal helpers
#' @keywords internal
#' @noRd
.amc_track <- function(tr) {
  if (is.list(tr) && !is.data.frame(tr) && !is.null(tr$position)) {
    pos <- tr$position
    P <- tr$period
  } else if (is.list(tr) && !is.data.frame(tr) && length(tr) >= 2 && (is.matrix(tr[[1]]) || is.data.frame(tr[[1]]))) {
    pos <- as.matrix(tr[[1]])[, 1]
    P <- tr[[2]]
  } else if (is.matrix(tr) || is.data.frame(tr)) {
    tr <- as.matrix(tr)
    pos <- tr[, 1]
    P <- tr[, -1, drop = FALSE]
  } else {
    stop("unrecognised track format")
  }
  P <- as.matrix(P)
  if (nrow(P) != length(pos) && ncol(P) == length(pos)) {
    P <- t(P)
  }
  if (nrow(P) != length(pos)) {
    stop("track: number of positions and period rows differ")
  }
  o <- order(pos)
  pos <- pos[o]
  P <- P[o, , drop = FALSE]
  bad <- !is.finite(P) | P <= 0
  if (any(bad)) {
    for (j in which(colSums(bad) > 0)) {
      k <- !bad[, j]
      if (sum(k) < 2) {
        stop("track realisation with < 2 valid periods")
      }
      P[, j] <- approx(pos[k], P[k, j], pos, rule = 2)$y
    }
  }
  list(position = pos, period = P)
}
#' @keywords internal
#' @noRd
.amc_interp <- function(x, Y, xout) {                       # column-wise linear interpolation, constant beyond the ends
  Y <- as.matrix(Y)
  o <- order(x)
  x <- x[o]
  Y <- Y[o, , drop = FALSE]
  d <- !duplicated(x)
  x <- x[d]
  Y <- Y[d, , drop = FALSE]
  if (length(x) == 1) {
    return(Y[rep(1, length(xout)), , drop = FALSE])
  }
  xo <- pmin(pmax(xout, x[1]), x[length(x)])
  i <- findInterval(xo, x, rightmost.closed = TRUE)
  i[i >= length(x)] <- length(x) - 1
  w <- (xo - x[i]) / (x[i + 1] - x[i])
  Y[i, , drop = FALSE] * (1 - w) + Y[i + 1, , drop = FALSE] * w
}
#' @keywords internal
#' @noRd
.amc_cumfloat <- function(z, P, target) {                    # floating time (kyr) from the top of the grid, trapezoid in 1/period
  P <- as.matrix(P)
  n <- length(z)
  if (n < 2) {
    return(matrix(0, 1, ncol(P)))
  }
  inc <- diff(z) * target / ((P[-1, , drop = FALSE] + P[-n, , drop = FALSE]) / 2)
  for (j in seq_len(ncol(inc))) {
    inc[, j] <- cumsum(inc[, j])
  }
  rbind(0, matrix(inc, nrow = n - 1))
}
#' @keywords internal
#' @noRd
.amc_eps <- function(eps, FT, PT, zg, zh, bounds, target, sgn) {
  nseg <- length(FT)
  if (is.numeric(eps) && length(eps) == 1) {
    return(list(source = "user (constant)", cycles = NULL, profile = rep(eps, length(zg))))
  }
  if (is.data.frame(eps)) {
    e <- approx(sgn * eps$position, eps$eps, zg, rule = 2)$y
    return(list(source = "user (profile)", cycles = NULL, profile = e))
  }
  eps_seg <- function(s) {
    f <- FT[[s]]
    tm <- rowMeans(f$mean)
    nodes <- approx(tm, f$z, seq(0, max(tm), by = target), ties = mean)$y
    nodes <- unique(nodes[is.finite(nodes)])
    if (length(nodes) < 2) {
      nodes <- c(bounds[s], bounds[s + 1])
    }
    if (!is.null(PT) && length(PT[[s]]) >= 2) {
      d <- matrix(NA_real_, length(nodes) - 1, length(PT[[s]]))
      for (r in seq_along(PT[[s]])) {
        x <- PT[[s]][[r]]
        d[, r] <- diff(approx(x$z, x$T, nodes, rule = 2)$y)
      }
      ld <- log(pmax(d, 1e-6))
      e <- numeric(nrow(ld))
      for (r in seq_len(nrow(ld))) {
        e[r] <- sd(ld[r, ])
      }
      e <- e / sqrt(ncol(d))
      src <- "proxy se"
    } else {
      d <- matrix(NA_real_, length(nodes) - 1, ncol(f$T))
      for (r in seq_len(ncol(f$T))) {
        d[, r] <- diff(approx(f$z, f$T[, r], nodes, rule = 2)$y)
      }
      ld <- log(pmax(d, 1e-6))
      e <- numeric(nrow(ld))
      for (r in seq_len(nrow(ld))) {
        e[r] <- sd(ld[r, ])
      }
      src <- "track realisations"
    }
    data.frame(segment = s, top = nodes[-length(nodes)], bottom = nodes[-1], eps = e, source = src)
  }
  cyc_list <- vector("list", nseg)
  for (s in seq_len(nseg)) {
    cyc_list[[s]] <- eps_seg(s)
  }
  cyc <- do.call(rbind, cyc_list)
  cyc$eps[!is.finite(cyc$eps)] <- median(cyc$eps, na.rm = TRUE)
  if (all(!is.finite(cyc$eps))) {
    cyc$eps <- 0.05
  }
  prof <- numeric(length(zg))
  sz <- findInterval(zg, zh, left.open = TRUE) + 1
  for (s in seq_len(nseg)) {
    k <- cyc$segment == s
    u <- sz == s
    if (!any(u)) {
      next
    }
    m <- (cyc$top[k] + cyc$bottom[k]) / 2
    prof[u] <- if (sum(k) > 1) { approx(m, cyc$eps[k], zg[u], rule = 2)$y } else { cyc$eps[k] }
  }
  cyc$top <- sgn * cyc$top
  cyc$bottom <- sgn * cyc$bottom
  list(source = unique(cyc$source), cycles = cyc, profile = prof)
}
#' @keywords internal
#' @noRd
.amc_bchron <- function(dates, zg, zh, sgn, sc, split, args, verbose) {
  run1 <- function(d, zp) {
    a <- list(ages = d$age * sc, ageSds = d$sd * sc, positions = -(sgn * d$position), positionThicknesses = abs(d$thickness), ids = d$id,
        predictPositions = -zp, MC = 10000, burn = 2000)                 # modBchron works in height (increasing upwards), kyr
    for (nm in setdiff(names(args), c("truncateUp", "truncateDown"))) {
      a[[nm]] <- args[[nm]]
    }
    if (!is.null(args$truncateUp)) {
      a$truncateUp <- args$truncateUp * sc
    }          # truncation ages in age_unit
    if (!is.null(args$truncateDown)) {
      a$truncateDown <- args$truncateDown * sc
    }
    if (verbose) {
      cat(sprintf("modBchron: %d dates, %d iterations\n", nrow(d), a$MC))
    }
    m <- if (verbose) {
      do.call(modBchron, a)
    } else {
      utils::capture.output(r <- do.call(modBchron, a))
      r
    }
    keep <- m$burn:ncol(m$model)
    list(fit = m, curves = m$model[, keep, drop = FALSE])
  }
  if (!split || !length(zh)) {
    r <- run1(dates, zg)
    C <- r$curves
    fits <- list(r$fit)
    hz <- NULL
  } else {
    s_d <- findInterval(sgn * dates$position, zh, left.open = TRUE) + 1
    s_g <- findInterval(zg, zh, left.open = TRUE) + 1
    C <- NULL
    fits <- list()
    edge <- list()
    for (s in sort(unique(s_g))) {
      d <- dates[s_d == s, ]
      if (nrow(d) < 2) {
        warning(sprintf("segment %d has < 2 dates; modBchron run with all dates there", s))
        d <- dates
      }
      zp <- sort(unique(c(zg[s_g == s], if (s > 1) { zh[s - 1] }, if (s <= length(zh)) { zh[s] })))
      r <- run1(d, zp)
      fits[[s]] <- r$fit
      n <- min(ncol(r$curves), if (is.null(C)) { Inf } else { ncol(C) })
      edge[[s]] <- list(z = zp, curves = r$curves)
      cs <- r$curves[match(zg[s_g == s], zp), seq_len(n), drop = FALSE]
      C <- if (is.null(C)) { cs } else { rbind(C[, seq_len(n), drop = FALSE], cs) }
    }
    hz_one <- function(h) {
      up <- edge[[h]]
      lo <- edge[[h + 1]]
      n <- min(ncol(up$curves), ncol(lo$curves))
      d <- lo$curves[match(zh[h], lo$z), seq_len(n)] - up$curves[match(zh[h], up$z), seq_len(n)]
      d <- d[d > 0]
      q <- quantile(d, c(.5, .025, .975))
      data.frame(position = sgn * zh[h], median_kyr = q[[1]], plus_kyr = q[[3]] - q[[1]], minus_kyr = q[[1]] - q[[2]])
    }
    hz_list <- vector("list", length(zh))
    for (h in seq_along(zh)) {
      hz_list[[h]] <- hz_one(h)
    }
    hz <- do.call(rbind, hz_list)
    if (verbose) {
      for (h in seq_len(nrow(hz))) {
        cat(sprintf("modBchron (split): hiatus at %g = %.0f +%.0f/-%.0f kyr\n", hz$position[h], hz$median_kyr[h], hz$plus_kyr[h], hz$minus_kyr[h]))
      }
    }
  }
  qq <- matrix(NA_real_, 3, nrow(C), dimnames = list(names(quantile(0, c(.025, .5, .975))), rownames(C)))
  sd_c <- numeric(nrow(C))
  for (r in seq_len(nrow(C))) {
    qq[, r] <- quantile(C[r, ], c(.025, .5, .975), na.rm = TRUE)
    sd_c[r] <- sd(C[r, ], na.rm = TRUE)
  }
  if (length(rownames(C))) {
    names(sd_c) <- rownames(C)
  }
  list(fits = fits, model = C / sc, predictPositions = sgn * zg, sd2_kyr = 2 * sd_c,
      summary = data.frame(position = sgn * zg, mean = rowMeans(C, na.rm = TRUE) / sc, median = qq[2, ] / sc, lower = qq[1, ] / sc, upper = qq[3, ] / sc,
          sd2_kyr = 2 * sd_c),
      hiatus = hz, split = split)
}

# Methods
#' @export
print.astroMCBchron <- function(x, ...) {
  s <- x$settings
  u <- s$age_unit
  cat(sprintf("astroMCBchron model: %d dates, %d hiatus(es); anchor %s; tau from %s; %d posterior draws, %d curves stored\n",
      nrow(x$dates), NCOL(x$hiatus), s$anchor, s$tau_source, length(x$posterior), s$curves_stored))
  d <- x$dates
  show <- data.frame(id = d$id, position = d$position, date = round(d$date, 3), posterior = round(d$posterior, 3),
      `2sd_kyr` = round(d$posterior_2sd_kyr), shift_kyr = round(d$shift_kyr), shift_sd = round(d$shift_sd, 2), outlier = d$outlier, check.names = FALSE)
  if (!is.null(d$p_outlier) && any(is.finite(d$p_outlier))) {
    show$p_outlier <- round(d$p_outlier, 2)
  }
  if (!is.null(d$group) && length(unique(d$group)) > 1) {
    show$group <- d$group
  }
  if (!is.null(d$loo_z)) {
    show$loo_z <- round(d$loo_z, 2)
  }
  dm <- if (is.null(s$date_model)) { "two_stage" } else { s$date_model }
  cat(sprintf("date model: %s%s\n", dm, if (dm == "mixture") {
    sprintf(" (%s outliers); outlier probability %.3f [%.3f-%.3f]", s$outlier_type,
        median(x$outlier_prob), quantile(x$outlier_prob, .025), quantile(x$outlier_prob, .975))
  } else {
    ""
  }))
  cat(sprintf("\nDates (%s) and posterior ages:\n", u))
  print(show, row.names = FALSE)
  cat("\nStretch factors (median [95%]):\n")
  s_med <- s_lo <- s_up <- numeric(ncol(x$stretch))
  for (j in seq_len(ncol(x$stretch))) {
    s_med[j] <- median(x$stretch[, j])
    s_lo[j] <- quantile(x$stretch[, j], .025)
    s_up[j] <- quantile(x$stretch[, j], .975)
  }
  if (length(colnames(x$stretch))) {
    names(s_med) <- names(s_lo) <- names(s_up) <- colnames(x$stretch)
  }
  print(data.frame(interval = colnames(x$stretch), median = round(s_med, 2), lower = round(s_lo, 2),
      upper = round(s_up, 2)), row.names = FALSE)
  if (NCOL(x$hiatus)) {
    for (h in seq_len(ncol(x$hiatus))) {
      cat(sprintf("Hiatus at %s: %.0f +/- %.0f kyr (2 sd)\n", sub("hiatus_", "", colnames(x$hiatus)[h]),
          mean(x$hiatus[, h]), 2 * sd(x$hiatus[, h])))
    }
  }
  cat(sprintf("tau0: %.3f [%.3f-%.3f]\n", median(x$tau0), quantile(x$tau0, .025), quantile(x$tau0, .975)))
  if (!is.null(x$group_offsets)) {
    for (g in setdiff(colnames(x$group_offsets), s$reference_group)) {
      cat(sprintf("laboratory offset of group %s relative to %s: %.0f [%.0f, %.0f] kyr\n", g, s$reference_group, median(x$group_offsets[, g]),
          quantile(x$group_offsets[, g], .025), quantile(x$group_offsets[, g], .975)))
    }
  }
  if (!is.null(x$depth_offsets)) {
    for (g in colnames(x$depth_offsets)) {
      cat(sprintf("depth offset of group %s: %.3f [%.3f, %.3f]\n", g, median(x$depth_offsets[, g]), quantile(x$depth_offsets[, g], .025), quantile(x$depth_offsets[, g], .975)))
    }
  }
  if (!is.null(x$summary$sd2_total_kyr)) {
    cat(sprintf("2 sd (mean over the record): internal %.0f kyr; + tracer %.0f kyr; + decay constant %.0f kyr\n",
        mean(x$summary$sd2_kyr), mean(x$summary$sd2_tracer_kyr), mean(x$summary$sd2_total_kyr)))
  }
  if (!is.null(x$comparison)) {
    k <- x$comparison$position >= min(x$dates$position) & x$comparison$position <= max(x$dates$position)
    cat(sprintf("2 sd reduction vs modBchron between the outermost dates: mean %.0f%%, median %.0f%%\n", mean(x$comparison$reduction_pct[k], na.rm = TRUE),
        median(x$comparison$reduction_pct[k], na.rm = TRUE)))
  }
  invisible(x)
}
#' @export
summary.astroMCBchron <- function(object, ...) object$summary

#' Predict ages at new positions from an astroMCBchron model
#'
#' @param object an \code{astroMCBchron} object.
#' @param newPositions positions (same type and unit as the model).
#' @param side for positions exactly at a hiatus: "upper" (younger side, default) or "lower".
#' @param ... unused.
#' @return list with \code{raw} (draws x positions, in age units) and \code{HDI}
#'   (position, 2.5, 50 and 97.5 percentiles).
#' @export
predict.astroMCBchron <- function(object, newPositions, side = c("upper", "lower"), ...) {
  side <- match.arg(side)
  sgn <- if (object$settings$position_type == "depth") { 1 } else { -1 }
  zg <- sgn * object$predictPositions
  zn <- sgn * newPositions
  zh <- if (NCOL(object$hiatus)) { sgn * object$settings$hiatuses$position } else { numeric(0) }
  sg <- findInterval(zg, zh, left.open = TRUE)
  sn <- findInterval(zn, zh, left.open = side == "upper")
  if (any(zn < min(zg) | zn > max(zg))) {
    stop("newPositions outside the model range")
  }
  raw <- matrix(NA_real_, ncol(object$model), length(zn))
  for (q in unique(sn)) {
    u <- sn == q
    g <- sg == q
    raw[, u] <- t(.amc_interp(zg[g], object$model[g, , drop = FALSE], zn[u]))
  }
  hq <- matrix(NA_real_, ncol(raw), 3, dimnames = list(NULL, names(quantile(0, c(.025, .5, .975)))))
  for (j in seq_len(ncol(raw))) {
    hq[j, ] <- quantile(raw[, j], c(.025, .5, .975), na.rm = TRUE)
  }
  list(raw = raw, HDI = data.frame(position = newPositions, hq, check.names = FALSE))
}

#' Plot an astroMCBchron model
#'
#' Three panels: (A) age model (95\% range) with the dates and, if present, the modBchron
#' model; (B) 2 sd age uncertainty of both models; (C) posterior and leave-one-out ages
#' minus the dates. Oldest strata/ages on the left and bottom.
#' @param x an \code{astroMCBchron} object.
#' @param ... unused.
#' @export
plot.astroMCBchron <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  s <- x$summary
  d <- x$dates
  u <- x$settings$age_unit
  dep <- x$settings$position_type == "depth"
  xl <- if (dep) { rev(range(s$position)) } else { range(s$position) }                # oldest strata on the left
  graphics::layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(1.3, 1))
  graphics::par(mar = c(4, 4.5, 1.5, 1), mgp = c(2.4, .7, 0))
  yl <- rev(range(c(s$lower, s$upper, d$date)))                               # oldest ages at the bottom
  plot(NA, xlim = xl, ylim = yl, xlab = sprintf("Position (%s)", x$settings$position_type), ylab = sprintf("Age (%s)", u))
  if (!is.null(x$bchron)) {
    b <- x$bchron$summary
    graphics::polygon(c(b$position, rev(b$position)), c(b$lower, rev(b$upper)), col = grDevices::adjustcolor("grey60", .35), border = NA)
    graphics::lines(b$position, b$median, lty = 2, col = "grey30")
  }
  graphics::polygon(c(s$position, rev(s$position)), c(s$lower, rev(s$upper)), col = grDevices::adjustcolor("steelblue", .45), border = NA)
  graphics::lines(s$position, s$median, col = "navy", lwd = 1.6)
  if (!is.null(s$lower_total)) {
    graphics::lines(s$position, s$lower_total, lty = 3, col = "navy")
    graphics::lines(s$position, s$upper_total, lty = 3, col = "navy")
  }
  k2 <- d$date_2sd_kyr / (if (u == "Ma") { 1000 } else { 1 })
  graphics::arrows(d$position, d$date - k2, d$position, d$date + k2,
      angle = 90, code = 3, length = .03, col = "firebrick", lwd = 1.4)
  graphics::points(d$position, d$date, pch = 21, bg = ifelse(d$outlier, "white", "firebrick"), cex = .9)
  if (!is.null(d$p_outlier) && any(d$p_outlier > 0.2, na.rm = TRUE)) {
    k <- which(d$p_outlier > 0.2)
    graphics::text(d$position[k], d$date[k], sprintf("p = %.2f", d$p_outlier[k]), pos = 3, cex = .65, col = "firebrick")
  }
  if (NCOL(x$hiatus)) {
    graphics::abline(v = x$settings$hiatuses$position, lty = 4, col = "darkorange")
  }
  graphics::legend(if (dep) { "topleft" } else { "topright" }, bty = "n", cex = .8, legend = c("astroMCBchron (median, 95%)", if (!is.null(x$bchron)) { "modBchron (median, 95%)" }, "dates (2 sd)"),
      col = c("navy", if (!is.null(x$bchron)) { "grey30" }, "firebrick"), lty = c(1, if (!is.null(x$bchron)) { 2 }, NA), pch = c(NA, if (!is.null(x$bchron)) { NA }, 21),
      pt.bg = "firebrick", lwd = c(1.6, if (!is.null(x$bchron)) { 1 }, NA))
  graphics::mtext("A", 3, .2, adj = 0, font = 2)
  ym <- max(c(s$sd2_kyr, s$sd2_tracer_kyr, if (!is.null(x$bchron)) { x$bchron$summary$sd2_kyr }, d$date_2sd_kyr), na.rm = TRUE)
  plot(s$position, s$sd2_kyr, type = "l", col = "navy", lwd = 1.6, xlim = xl, ylim = c(0, ym * 1.05), xlab = sprintf("Position (%s)", x$settings$position_type), ylab = "2 sd (kyr)")
  if (!is.null(x$bchron)) {
    graphics::lines(x$bchron$summary$position, x$bchron$summary$sd2_kyr, lty = 2, col = "grey30")
  }
  if (!is.null(s$sd2_tracer_kyr)) {
    graphics::lines(s$position, s$sd2_tracer_kyr, lty = 3, col = "navy")
    graphics::legend("topleft", bty = "n", cex = .7, legend = c("internal", "+ tracer (dotted)"), lty = c(1, 3), col = "navy")
  }
  graphics::points(d$position, d$date_2sd_kyr, pch = 23, bg = "firebrick", cex = .8)
  graphics::mtext("B", 3, .2, adj = 0, font = 2)
  o <- order(d$position, decreasing = dep)
  dd <- d[o, ]
  n <- nrow(dd)
  xx <- seq_len(n)
  graphics::par(mar = c(7.5, 4.5, 1.5, 1))
  yr <- range(c(dd$shift_kyr - dd$posterior_2sd_kyr, dd$shift_kyr + dd$posterior_2sd_kyr, -dd$date_2sd_kyr, dd$date_2sd_kyr,
      if (!is.null(dd$loo_pred)) { c((dd$loo_pred - dd$date) * (if (u == "Ma") { 1000 } else { 1 }) + c(-1, 1) %o% dd$loo_2sd_kyr) }), na.rm = TRUE)
  plot(NA, xlim = c(.5, n + .5), ylim = yr, xaxt = "n", xlab = "", ylab = "Model age - date (kyr)")
  graphics::axis(1, xx, labels = dd$id, las = 2, cex.axis = .7)
  graphics::rect(xx - .35, -dd$date_2sd_kyr, xx + .35, dd$date_2sd_kyr, col = grDevices::adjustcolor("firebrick", .15), border = NA)
  graphics::abline(h = 0, col = "grey50")
  graphics::arrows(xx - .1, dd$shift_kyr - dd$posterior_2sd_kyr, xx - .1, dd$shift_kyr + dd$posterior_2sd_kyr, angle = 90, code = 3, length = .02, col = "navy")
  graphics::points(xx - .1, dd$shift_kyr, pch = 19, col = "navy", cex = .7)
  if (!is.null(dd$loo_pred)) {
    lp <- (dd$loo_pred - dd$date) * (if (u == "Ma") { 1000 } else { 1 })
    graphics::arrows(xx + .1, lp - dd$loo_2sd_kyr, xx + .1, lp + dd$loo_2sd_kyr, angle = 90, code = 3, length = .02, col = "darkgreen")
    graphics::points(xx + .1, lp, pch = 17, col = "darkgreen", cex = .7)
  }
  graphics::mtext("C", 3, .2, adj = 0, font = 2)
  invisible(x)
}
