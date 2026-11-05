#' @title Test sedimentation rates with noise-corrected TimeOpt, likelihood and COCO scores
#'
#' @description Scans a grid of constant sedimentation rates for one stratigraphic record and
#' scores every rate in three ways: the TimeOpt envelope and spectral fits (Meyers, 2015), a
#' likelihood-based evidence score with an information-criterion penalty (in the spirit of
#' Malinverno and Meyers, 2024) and the COCO correlation coefficient (Li et al., 2018). Every
#' score is expressed against noise analysed at the same rate, rates that cannot be tested are
#' masked, astronomical target periods can be taken at the age of the record from AstroGeo22
#' (Farhat et al., 2022), the TimeOpt and COCO scores are combined into one score, and the result
#' is repeated for several noise models.
#'
#' @details
#' \strong{Why a same-rate noise correction.} The raw TimeOpt r-squared values (and COCO
#' correlations) that pure noise reaches depend on the trial rate: at high rates the record spans
#' a short time, fewer cycles fit, and the regression has fewer independent data points per
#' fitted sinusoid; for red noise the target frequencies also move to the strong low-frequency end
#' of the spectrum. Li et al. (2018) describe the same effect for COCO (fewer astronomical terms
#' resolved at high rates). Here every statistic is compared with the same statistic of noise
#' series analysed at the same rate:
#' \deqn{z(u) = (S(u) - \mu_{noise}(u)) / \sigma_{noise}(u)}
#' (r-squared values on a logit scale, COCO correlations on a Fisher z scale; the noise mean and
#' standard deviation are smoothed over neighbouring rates). The family-wise probability of a
#' score is the fraction of noise series whose maximum z over all testable rates is at least as
#' large (as in the full-grid test of Meyers and Malinverno, 2018), with a floor of
#' 1/(n_null + 1).
#'
#' \strong{TimeOpt scores.} For each rate the record is calibrated to time, the carrier band
#' (precession for fit 1, short eccentricity for fit 2) is bandpassed with
#' \code{\link[astrochron:taner]{taner}}, its amplitude envelope is taken with
#' \code{\link[astrochron:hilbert]{hilbert}} and regressed on the modulating periods
#' (eccentricity for fit 1, the 405-kyr term for fit 2), and the record itself is regressed on the
#' target periods. r2_envelope, r2_power and r2_opt follow Meyers (2015); z_opt combines the
#' z-scores of the envelope and power fits (Stouffer's method).
#' 
#' \strong{Significance spectrum (main TimeOpt score).} S_opt refers the envelope and the power fit
#' each to the noise distribution at the same rate (the noise series of that rate and of the
#' \code{null_smooth} neighbouring rates pooled; an exponential tail above the pooled 75 \% quantile
#' resolves p below 1 / n_null), combines the two p-values (Fisher) and calibrates the combination
#' the same way; S_opt = -log10 of that p. The family-wise p (p_S) uses the minimum-p method
#' (Westfall and Young, 1993): the share of noise series whose smallest p over the testable rates is
#' at least as small. Because every rate is referred to its own noise distribution, noise gives a
#' flat baseline near zero over the whole range (no drift to fast rates as in raw r2_opt, and no
#' inflated low-rate scatter as in z_opt), and the family-wise levels are horizontal lines.
#' 
#' \strong{Tuned-band peaks and joint score.} TimeOpt's envelope fit has little power when the
#' record spans few modulating cycles and the rate changes (in a Kosov-like test, about 2 Myr with
#' +-17 \% rate changes, the fit-2 envelope reached p < 0.05 at the true rate in 15-30 \% of records
#' with a known signal). With \code{tuned_peaks = TRUE} a second score is computed at every rate: the
#' record is minimally tuned on the 405-kyr cycle (fit 2) or on the short eccentricity (fit 1) at that
#' rate (as \code{tune = "minimal"}), and the bands not used for the tuning (short eccentricity,
#' obliquity, precession for fit 2; obliquity and precession for fit 1) are tested for spectral peaks:
#' the mean over these bands of the largest log ratio of the periodogram to a smooth continuum
#' (cubic in log frequency fitted to the log periodogram with the peaks left out) within
#' \code{peak_tol} of each period. Every noise series is tuned and scored in exactly the same way, so
#' the regularity that the tuning imposes is part of the null and the test is not circular; the tuned
#' band itself is not scored. S_peaks is -log10 of the same-rate p with the minimum-p family-wise
#' p_peaks, as for S_opt. In the Kosov-like test this score reached p < 0.05 in 90-100 \% of records
#' with a known signal. S_joint combines S_opt and S_peaks (Fisher) and is calibrated again on the
#' same noise series (p_joint); use it as the main summary when both are available.
#'
#' \strong{Likelihood score.} This score is complementary to TimeOpt and follows ideas of
#' \code{\link[astrochron:timeOptB]{timeOptB}} (Malinverno and Meyers, 2024): the residuals of the
#' spectral fit are modelled as an AR(2) process (the likelihood is that of the AR(2) innovations),
#' and the envelope fit uses an effective number of independent points
#' N_eff = N / tau, with tau the first lag at which the autocorrelation of the envelope reaches
#' zero. At every rate the astronomical model is compared with a model without astronomical
#' terms, and an information-criterion penalty is applied for the number of resolvable terms:
#' \deqn{B = 2 (\log L_{astro} - \log L_{noise}) - k \log N}
#' (BIC, \code{penalty = "BIC"}) or \eqn{- 2k} (AIC). Only target terms that the record can
#' resolve at that rate enter the likelihood model: terms with fewer than \code{min_samples} samples
#' per cycle or with periods longer than the record are dropped, and terms closer in frequency
#' than the Rayleigh resolution (1 / duration) are merged into one. B > 0 means the astronomical
#' model is preferred at that rate; \code{weight_lik} gives rate weights proportional to
#' exp(B / 2). The likelihood score is a model comparison per rate, not a full Bayesian posterior:
#' the envelope differs between rates, so the weights are indicative. Use
#' \code{\link[astrochron:timeOptB]{timeOptB}} for posterior distributions of the sedimentation
#' rate and the precession frequency.
#'
#' \strong{COCO.} The correlation coefficient between the AR(1)-corrected periodogram of the
#' record and an astronomical target periodogram (\code{\link{coco_classic}},
#' \code{\link{coco_target_from_periods}}; Li et al., 2018).
#'
#' \strong{Combined scores.} z_comb is the sum of z_opt and z_coco, standardised again against the
#' same noise series, so the correlation between the TimeOpt and COCO scores is accounted for.
#' z_cons (consensus) is the smaller of the two z-scores, standardised in the same way: it is only
#' high where both methods support the rate, which suppresses peaks that only one method shows
#' (e.g. the COCO alias at ~4.3 times the true rate, where eccentricity cycles are read as
#' precession and short eccentricity). With \code{combine = "all"} the likelihood score is added to
#' both.
#'
#' \strong{Noise models.} AR(1), AR(2), power law (f^-beta) and bending power law
#' (1 / (1 + (f / f_b)^beta); Vaughan et al., 2011), fitted to the periodogram of the record
#' (Whittle likelihood). \code{noise_model = "auto"} uses the model with the lowest AIC.
#' \code{sensitivity = TRUE} repeats the test for every noise model.
#'
#' \strong{Noise fitted without the cycles} (\code{null_fit}). A noise model fitted to the whole
#' periodogram of the record also fits the astronomical cycles in it, and its shape near the
#' tested frequencies is set mainly by the many high-frequency periodogram values. Both distort
#' the background at the frequencies that matter: cycles that are present bend the fitted
#' spectrum (e.g. an AR(1) background is fitted as AR(2)), so the null partly contains the cycles
#' it is tested against and genuine cycles look less significant. Therefore, by default
#' (\code{null_fit = "local"}), the noise model is fitted separately for every trial rate
#' (Whittle likelihood):
#' \itemize{
#'   \item leaving out the frequency bands where that rate puts the astronomical terms (405-kyr
#'   eccentricity, short eccentricity, obliquity and precession, each band widened by
#'   \code{mask_pad} Rayleigh frequencies 1 / duration), as in robust background estimates that
#'   exclude spectral peaks (Mann and Lees, 1996; Meyers, 2012);
#'   \item using only frequencies up to \code{local_factor} times the highest target frequency of
#'   the spectral fit (at least 40 frequencies), so the background is estimated where the tested
#'   cycles are;
#'   \item carrying the uncertainty of the fitted parameters (Whittle Hessian) into the noise
#'   series: every noise series draws its parameters from that uncertainty, so a background that
#'   is poorly constrained by few frequencies gives a wider null instead of false alarms.
#' }
#' The noise series for a rate come from the model fitted at that rate; the same random numbers
#' are used at all rates, so each noise series is one coherent "noise world" over the rate grid
#' and the family-wise probability keeps its meaning. The model choice
#' (\code{noise_model = "auto"}) uses the AIC summed over the testable rates.
#' \code{null_fit = "masked"} leaves out the bands but uses all other frequencies;
#' \code{null_fit = "full"} fits one noise model to the whole periodogram (earlier behaviour).
#'
#' \strong{Minimal tuning} (\code{tune = "minimal"}). A constant trial rate smears the cycles
#' when the true rate varies, so the TimeOpt fit can favour a wrong rate. With minimal tuning the
#' record is first tuned at every trial (mean) rate with the highest-order cycle only: for fit 2
#' the 405-kyr cycle, for fit 1 the short eccentricity (tied at the harmonic mean of its periods,
#' ~110 kyr). The cycle is bandpassed in depth around the period the trial rate predicts
#' (+- \code{tune_width}), every cycle of its (monotonic) Hilbert phase is tied to the cycle
#' period, time is interpolated linearly between the ties and extrapolated at the trial rate
#' beyond them, and TimeOpt, the likelihood score and COCO are computed on the tuned series. The
#' tuned cycle is left out of the spectral fit (\code{tune_fit_tuned = FALSE}) because the tuning
#' imposes its regularity. Every noise series is tuned in exactly the same way at every rate, so
#' the gain in fit that tuning gives to any series is part of the null and the probabilities
#' stay calibrated. \code{tuning} returns the minimal-tuning age model at the rate selected by
#' z_opt. The trial rate is then the mean rate of a tuned age model, so the test asks at which
#' scale a 405-kyr (or short-eccentricity) tuning gives the best astronomical fit.
#'
#' \strong{Envelope and band corrections (fit 2).} Three settings protect the modulation fit
#' of TimeOpt (fit 2 in particular) against artefacts found on real records:
#' \itemize{
#'   \item \code{edge_trim}: the Taner bandpass and the Hilbert transform ring at the ends of a
#'   finite record, and the ringing extends over about one impulse length of the filter
#'   (1 / bandwidth, ~230 kyr for the short-eccentricity band). With only a few 405-kyr cycles
#'   these end values steer the envelope regression; they are removed before the fit (and from
#'   the duration used for the untestable mask).
#'   \item \code{band_guard}: the default fit-2 band (1/130.7 - 0.001 to 1/94.9 + 0.001
#'   cycles/kyr, i.e. 87-151 kyr) contains the s4-s6 obliquity-modulation term (151 kyr); the
#'   long-period edge is moved halfway to it (~140 kyr).
#'   \item \code{obliquity_mod}: the s3-s6 and s4-s6 terms are part of the spectral model, so
#'   their power counts at the rate where they belong. Otherwise a rate ~1.2-1.3 times too high,
#'   at which the short-eccentricity targets fall on the obliquity-modulation cycles, can win.
#' }
#' With \code{band_guard = FALSE, obliquity_mod = FALSE, edge_trim = 0} the r-squared values are
#' those of astrochron::timeOpt.
#'
#' \strong{Age-dependent targets.} With \code{age_ma}, the target periods are computed from the
#' AstroGeo22 precession frequency (Farhat et al., 2022; \code{\link{astrogeo22}}) and the La2004
#' secular frequencies (Laskar et al., 2004; held at their present values): eccentricity
#' g2-g5, g3-g2, g4-g2, g3-g5, g4-g5; precession k+g5, k+g2, k+g4, k+g3; obliquity k+s3,
#' k+s4, k+s6. \code{target_sensitivity = TRUE} repeats the scores with the AstroGeo22 minimum
#' and maximum precession frequency.
#'
#' \strong{Untestable rates} (\code{untestable}): 1 = fewer than \code{min_samples} samples per
#' shortest carrier period; 2 = the record spans fewer than \code{min_cycles} of the longest
#' modulating period.
#'
#' @param data Evenly spaced stratigraphic record (matrix or data frame with two columns: depth in metres in
#' column 1, proxy value in column 2; further columns are ignored). The rows are sorted by depth; the record
#' is detrended (see \code{detrend}) and standardised before the analysis. Required, no default
#' (\code{Default=NULL} gives an error).
#' @param sedmin,sedmax Lowest and highest trial sedimentation rate of the rate grid, in cm/kyr (numeric,
#' > 0, \code{sedmin} < \code{sedmax}). \code{Default=1} and \code{Default=20}.
#' @param numsed Number of trial sedimentation rates in the grid (integer, >= 1; 1 tests \code{sedmin} only).
#' \code{Default=100}.
#' @param linLog Spacing of the rate grid (numeric): \code{1} logarithmic, \code{0} linear. \code{Default=1}.
#' @param fit Which TimeOpt fit is used (numeric): \code{1} precession carrier modulated by eccentricity,
#' \code{2} short-eccentricity carrier modulated by the 405-kyr eccentricity. \code{Default=1}.
#' @param age_ma Age of the record in Ma (numeric, 0 to 4429), used to compute the precession and obliquity
#' target periods from the AstroGeo22 precession frequency (\code{\link{astrogeo22}}); \code{NULL} uses the
#' present-day precession frequency. Values given in \code{targetE}, \code{targetP} or \code{targetO}
#' override the computed periods. \code{Default=NULL}.
#' @param targetE,targetP,targetO Your own eccentricity, precession and obliquity target periods in kyr
#' (numeric vectors, any order; they are sorted from longest to shortest). For \code{fit = 2} the longest
#' value of \code{targetE} is the modulating (405-kyr) period and the others are the carrier periods.
#' \code{NULL} uses the La2004/AstroGeo22 periods (see \code{age_ma}). \code{Default=NULL}.
#' @param obliquity Include the obliquity periods in the spectral (power) fit and in the COCO target
#' (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param fitModPwr Include the modulating eccentricity periods in the spectral (power) fit (logical, TRUE or
#' FALSE), as in \code{astrochron::timeOpt}. \code{Default=TRUE}.
#' @param ftol Half-width added outside the carrier frequencies to define the Taner bandpass, in cycles/kyr
#' (numeric, >= 0); \code{NULL} uses 0.005 for the precession carrier (fit 1) and 0.001 for the
#' short-eccentricity carrier (fit 2). \code{Default=NULL}.
#' @param roll Roll-off rate of the Taner bandpass filter used for the carrier band (numeric, > 0; see
#' \code{astrochron::taner}). \code{Default=10^7}.
#' @param detrend Remove a linear trend from the proxy before the analysis (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param penalty Information-criterion penalty of the likelihood score (character): \code{"BIC"}
#' (k log N) or \code{"AIC"} (2k). \code{Default=c("BIC", "AIC")}, i.e. \code{"BIC"}.
#' @param coco Compute the COCO correlation score (and the scores that combine it) (logical, TRUE or FALSE).
#' \code{Default=TRUE}.
#' @param combine Which scores enter the combined scores z_comb and z_cons (character):
#' \code{"timeopt_coco"} (TimeOpt and COCO) or \code{"all"} (TimeOpt, COCO and the likelihood score).
#' \code{Default=c("timeopt_coco", "all")}, i.e. \code{"timeopt_coco"}.
#' @param noise_model Noise model of the null series (character): \code{"auto"} (the model with the lowest AIC),
#' \code{"ar1"}, \code{"ar2"}, \code{"powerlaw"} or \code{"bending"} (bending power law). With
#' \code{sensitivity = TRUE} the other models are also run, but the main result uses this one.
#' \code{Default=c("auto", "ar1", "ar2", "powerlaw", "bending")}, i.e. \code{"auto"}.
#' @param null_fit How the noise model is fitted to the periodogram (character): \code{"local"} at every rate
#' without the astronomical bands, near the tested frequencies and with its parameter uncertainty (see Details);
#' \code{"masked"} at every rate without the bands but with all other frequencies; \code{"full"} once, to the
#' whole periodogram. \code{Default=c("local", "masked", "full")}, i.e. \code{"local"}.
#' @param mask_pad Widening of every masked astronomical band on each side, in Rayleigh frequencies
#' (1 / duration) (numeric, >= 0); used when \code{null_fit} is \code{"local"} or \code{"masked"}.
#' \code{Default=2}.
#' @param local_factor Upper frequency limit of the \code{"local"} noise fit, in multiples of the highest target
#' frequency of the spectral fit (numeric, > 1; at least 40 frequencies are always used). \code{Default=4}.
#' @param tune Tuning applied at every trial rate before scoring (character): \code{"none"} (constant rate) or
#' \code{"minimal"} (tune with the 405-kyr cycle for fit 2, the short eccentricity for fit 1; see Details).
#' \code{Default=c("none", "minimal")}, i.e. \code{"none"}.
#' @param tune_width Relative half-width of the bandpass around the period predicted by the trial rate that
#' extracts the tuning cycle (numeric, between 0 and 1; also the allowed deviation of the local rate); used by
#' \code{tune = "minimal"} and by the tuned-band peak score. \code{Default=0.3}.
#' @param tune_fit_tuned Keep the tuned cycle in the spectral fit after minimal tuning (logical, TRUE or FALSE);
#' \code{FALSE} leaves it out because the tuning imposes its regularity. \code{Default=FALSE}.
#' @param band_guard Keep every other astronomical term out of the carrier bandpass (logical, TRUE or FALSE):
#' each band edge is moved at most halfway to the nearest term outside the carrier set (for fit 2 the ~151-kyr
#' s4-s6 obliquity-modulation term next to the 131-kyr short eccentricity). \code{Default=TRUE}.
#' @param obliquity_mod Include the obliquity-modulation terms s3-s6 (~173 kyr) and s4-s6 (~151 kyr) in the
#' spectral fit (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param edge_trim Time in kyr removed at both ends of the envelope before the envelope fit (numeric, >= 0, or
#' the character \code{"auto"}): \code{"auto"} removes one impulse length of the carrier filter
#' (1 / bandwidth), \code{0} removes nothing (as \code{astrochron::timeOpt}). \code{Default="auto"}.
#' @param tuned_peaks Compute the tuned-band peak score S_peaks and the joint score S_joint (see Details)
#' (logical, TRUE or FALSE). \code{Default=TRUE}.
#' @param peak_tol Relative tolerance around each tested period in the tuned-band peak score (numeric, > 0;
#' 0.05 = periods within 5\% of the target). \code{Default=0.05}.
#' @param keep_null Return the same-rate p of every noise series (matrix rates x noise series) for S_opt, S_peaks
#' and S_joint in \code{null_p}, as needed by \code{\link{astro_sedrate_multi}} (logical, TRUE or FALSE).
#' \code{Default=FALSE}.
#' @param sensitivity Repeat the test for all four noise models and report them in \code{sensitivity} (logical,
#' TRUE or FALSE); \code{FALSE} runs only the chosen model, which is faster. \code{Default=TRUE}.
#' @param target_sensitivity Repeat the scores with the AstroGeo22 minimum and maximum precession frequency
#' (logical, TRUE or FALSE; only used when \code{age_ma} is given). \code{Default=TRUE}.
#' @param n_null Number of noise series simulated per noise model (integer, > 0); the smallest family-wise
#' probability that can be reached is 1/(n_null + 1). \code{Default=200}.
#' @param null_smooth Half-width in grid steps of the smoothing of the noise mean and standard deviation over
#' neighbouring rates, and of the pooling of the noise series of neighbouring rates for S_opt (integer, >= 0).
#' \code{Default=2}.
#' @param alpha False-alarm level (numeric, between 0 and 1); stored in the \code{settings} of the result.
#' \code{Default=0.05}.
#' @param min_samples Minimum number of samples per shortest carrier period; rates with fewer are marked
#' untestable (code 1) (numeric, > 0). \code{Default=3}.
#' @param min_cycles Minimum number of cycles of the longest modulating period that the record must span at a
#' rate; rates with fewer are marked untestable (code 2) (numeric, > 0). \code{Default=2}.
#' @param run_multicore Simulate and score the noise series in parallel on several cores (logical, TRUE or
#' FALSE). \code{Default=FALSE}.
#' @param n_cores Number of cores used for the noise series (integer); \code{NULL} uses the option
#' \code{waverider.n_cores} or, with \code{run_multicore = TRUE}, all but one detected core (at most
#' \code{getOption("waverider.max_cores", 4)}). \code{Default=NULL}.
#' @param verbose Print progress messages (logical, TRUE or FALSE). \code{Default=FALSE}.
#'
#' @return A list of class \code{astro_sedrate_test} with
#' \describe{
#'   \item{scores}{data frame per rate: sedrate, untestable code, r2_envelope, r2_power, r2_opt,
#'   z_envelope, z_power, z_opt, p_opt, B_spectral, B_envelope, B_total, weight_lik, z_lik, p_lik,
#'   coco_rho, z_coco, p_coco, z_comb, p_comb, z_cons, p_cons, S_opt, p_S, S_peaks, p_peaks,
#'   S_joint, p_joint (p = family-wise probability; S = -log10 of the same-rate p)}
#'   \item{best}{data frame: best rate, value and family-wise p of each score (raw r2_opt included)}
#'   \item{null_p}{with \code{keep_null = TRUE}: same-rate p of every noise series (rates x noise series)
#'   for S_opt, S_peaks and S_joint; otherwise NULL}
#'   \item{significance_levels}{matrix (rows TimeOpt, peaks, joint): S values of the family-wise
#'   5 \% and 1 \% levels}
#'   \item{tuning}{for \code{tune = "minimal"}: depth, time (kyr), sedimentation rate (cm/kyr)
#'   and tie points of the minimal-tuning age model at the rate selected by S_opt}
#'   \item{corroboration}{value and family-wise p of every score at the rate selected by S_opt, and
#'   the ratio of each score's own best rate to that rate (to spot aliases)}
#'   \item{sensitivity}{best rate, z and p of every score for every noise model}
#'   \item{target_sensitivity}{best rates with the AstroGeo22 minimum / nominal / maximum
#'   precession frequency}
#'   \item{noise_fits}{fitted noise models with their AIC (summed over the testable rates for
#'   \code{null_fit = "local"} or \code{"masked"})}
#'   \item{noise_fits_rate}{for \code{null_fit = "local"} or \code{"masked"}: AIC and parameters of every noise
#'   model at every testable rate, and the fraction of periodogram frequencies kept}
#'   \item{targets}{the target periods used}
#'   \item{settings}{the settings}
#' }
#'
#' @references
#' Meyers, S.R. (2015). The evaluation of eccentricity-related amplitude modulation and bundling
#' in paleoclimate data: An inverse approach for astrochronologic testing and time scale
#' optimization. Paleoceanography 30, 1625-1640 <doi:10.1002/2015PA002850>
#'
#' Westfall, P.H. and Young, S.S. (1993). Resampling-Based Multiple Testing: Examples and Methods
#' for p-Value Adjustment. Wiley, New York, ISBN 978-0-471-55761-6.
#'
#' Meyers, S.R. and Malinverno, A. (2018). Proterozoic Milankovitch cycles and the history of the
#' solar system. PNAS 115, 6363-6368 <doi:10.1073/pnas.1717689115>
#'
#' Malinverno, A. and Meyers, S.R. (2024). Bayesian estimation of past astronomical frequencies,
#' lunar distance, and length of day from sediment cycles. Geochemistry, Geophysics, Geosystems 25,
#' e2023GC011176 <doi:10.1029/2023GC011176>
#'
#' Li, M., Kump, L.R., Hinnov, L.A. and Mann, M.E. (2018). Tracking variable sedimentation rates
#' and astronomical forcing in Phanerozoic paleoclimate proxy series with evolutionary correlation
#' coefficients and hypothesis testing. Earth and Planetary Science Letters 501, 165-179
#' <doi:10.1016/j.epsl.2018.08.041>
#'
#' Farhat, M., Auclair-Desrotour, P., Boue, G. and Laskar, J. (2022). The resonant tidal
#' evolution of the Earth-Moon distance. Astronomy and Astrophysics 665, L1
#' <doi:10.1051/0004-6361/202243445>
#'
#' Mann, M.E. and Lees, J.M. (1996). Robust estimation of background noise and signal detection
#' in climatic time series. Climatic Change 33, 409-445 <doi:10.1007/BF00142586>
#'
#' Meyers, S.R. (2012). Seeing red in cyclic stratigraphy: Spectral noise estimation for
#' astrochronology. Paleoceanography 27, PA3228 <doi:10.1029/2012PA002307>
#'
#' Vaughan, S., Bailey, R.J. and Smith, D.G. (2011). Detecting cycles in stratigraphic data:
#' Spectral analysis in the presence of red noise. Paleoceanography 26, PA4211
#' <doi:10.1029/2011PA002195>
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' res <- astro_sedrate_test(data = mag[1:400, ], sedmin = 0.3, sedmax = 1.5, numsed = 25,
#'     fit = 1, n_null = 10, sensitivity = FALSE, target_sensitivity = FALSE)
#' res$best
#' }
#' @export
#' @importFrom stats lm residuals sd median optimize optim fft acf ar.yw arima.sim cor approx
#' @importFrom stats rnorm quantile
#' @importFrom astrochron taner hilbert
astro_sedrate_test <- function(data = NULL,
    sedmin = 1,
    sedmax = 20,
    numsed = 100,
    linLog = 1,
    fit = 1,
    age_ma = NULL,
    targetE = NULL,
    targetP = NULL,
    targetO = NULL,
    obliquity = TRUE,
    fitModPwr = TRUE,
    ftol = NULL,
    roll = 10^7,
    detrend = TRUE,
    penalty = c("BIC", "AIC"),
    coco = TRUE,
    combine = c("timeopt_coco", "all"),
    noise_model = c("auto", "ar1", "ar2", "powerlaw", "bending"),
    null_fit = c("local", "masked", "full"),
    mask_pad = 2,
    local_factor = 4,
    tune = c("none", "minimal"),
    tune_width = 0.3,
    tune_fit_tuned = FALSE,
    band_guard = TRUE,
    obliquity_mod = TRUE,
    edge_trim = "auto",
    tuned_peaks = TRUE,
    peak_tol = 0.05,
    keep_null = FALSE,
    sensitivity = TRUE,
    target_sensitivity = TRUE,
    n_null = 200,
    null_smooth = 2,
    alpha = 0.05,
    min_samples = 3,
    min_cycles = 2,
    run_multicore = FALSE,
    n_cores = NULL,
    verbose = FALSE) {
  penalty <- match.arg(penalty)
  combine <- match.arg(combine)
  noise_model <- match.arg(noise_model)
  null_fit <- match.arg(null_fit)
  tune <- match.arg(tune)
  if (!(fit %in% c(1, 2))) {
    stop("astro_sedrate_test: fit must be 1 or 2", call. = FALSE)
  }
  dat <- unname(as.matrix(data))[, 1:2]
  dat <- dat[order(dat[, 1]), , drop = FALSE]
  dd <- diff(dat[, 1])
  dx <- stats::median(dd)
  if (!all(abs(dd - dx) <= 1e-3 * abs(dx))) {
    stop("astro_sedrate_test: the depth values must be evenly spaced", call. = FALSE)
  }
  n <- nrow(dat)
  y <- dat[, 2]
  if (detrend) {
    y <- stats::residuals(stats::lm(y ~ dat[, 1]))
  }
  y <- (y - mean(y)) / stats::sd(y)
  sedrate <- wr_sedrate_grid(sedmin, sedmax, numsed, linLog)
  numsed <- length(sedrate)

  tg <- wr_ast_targets(age_ma, targetE, targetP, targetO)
  cfg <- wr_ast_cfg(tg, fit, obliquity, fitModPwr, ftol, roll, penalty, min_samples, coco,
      band_guard, obliquity_mod, edge_trim)
  cfg <- wr_ast_add_tune(cfg, tg, tune, tune_width, tune_fit_tuned)
  code <- wr_ast_testable(sedrate, dx, n, cfg, min_samples, min_cycles)
  testable <- code == 0L
  if (!any(testable)) {
    stop("astro_sedrate_test: no testable rate in sedmin-sedmax (too few samples per carrier cycle, ",
        "or the record spans fewer than min_cycles modulating cycles; see Details)", call. = FALSE)
  }

  # Observed scores
  obs <- wr_ast_eval(y, dx, sedrate, cfg)

  # Noise models
  fits_rate <- NULL
  nfr <- NULL
  if (null_fit == "full") {
    nfit <- wr_ast_noise_fit(y)
  } else {
    # noise model fitted at every testable rate without the astronomical bands of that rate
    fits_rate <- vector("list", numsed)
    for (j in which(testable)) {
      keep <- wr_ast_band_mask(n, dx * 100 / sedrate[j], tg, mask_pad,
          fmax = if (null_fit == "local") { local_factor / min(cfg$spec) } else { Inf })
      fits_rate[[j]] <- wr_ast_noise_fit(y, keep)
      f1 <- fits_rate[[j]]
      for (q in seq_len(nrow(f1))) {
        pq <- c(f1$param[[q]], NA, NA)[1:2]
        nfr <- rbind(nfr, data.frame(sedrate = sedrate[j], model = f1$model[q], AIC = f1$AIC[q],
            par1 = pq[1], par2 = pq[2], kept = mean(keep)))
      }
    }
    nfit <- fits_rate[[which(testable)[1]]][, c("model", "npar")]
    for (q in seq_len(nrow(nfit))) {
      sel <- nfr$model == nfit$model[q]
      nfit$AIC[q] <- sum(nfr$AIC[sel])
      nfit$description[q] <- sprintf("%s fitted per rate without astronomical bands; median %s",
          nfit$model[q], wr_ast_par_text(nfit$model[q],
          c(stats::median(nfr$par1[sel]), stats::median(nfr$par2[sel]))))
    }
  }
  chosen <- noise_model
  if (noise_model == "auto") {
    chosen <- nfit$model[which.min(nfit$AIC)]
  }
  models <- chosen
  if (sensitivity) {
    models <- unique(c(chosen, nfit$model))
  }
  pk <- NULL
  obs_pk <- rep(NA_real_, numsed)
  if (tuned_peaks) {
    pk <- list(cfg = wr_ast_add_tune(cfg, tg, "minimal", tune_width, FALSE), fam = wr_ast_peak_families(tg, fit),
        tol = peak_tol)
    for (j in which(testable)) {
      obs_pk[j] <- wr_ast_tuned_peaks(y, dx, sedrate[j], pk)
    }
  }
  runs <- list()
  sims_main <- NULL
  for (mdl in models) {
    if (verbose) {
      message("astro_sedrate_test: ", n_null, " noise series (", mdl, ")")
    }
    if (null_fit == "full") {
      par_m <- nfit$param[[mdl]]
      pk_m <- if (mdl == chosen) { pk } else { NULL }
      sims <- wr_lapply(seq_len(n_null), function(k) {
        xs <- wr_ast_noise_series(n, mdl, par_m)
        o <- wr_ast_eval(xs, dx, sedrate, cfg)
        o$peaks <- rep(NA_real_, numsed)
        if (!is.null(pk_m)) {
          for (j in which(testable)) {
            o$peaks[j] <- wr_ast_tuned_peaks(xs, dx, sedrate[j], pk_m)
          }
        }
        o
      }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
    } else {
      par_r <- vector("list", numsed)
      for (j in which(testable)) {
        par_r[[j]] <- list(par = fits_rate[[j]]$param[[mdl]], V = fits_rate[[j]]$theta_cov[[mdl]])
      }
      sims <- wr_lapply(seq_len(n_null), function(k) {
        wr_ast_eval_rates(n, dx, sedrate, cfg, testable, mdl, par_r, pk = if (mdl == chosen) { pk } else { NULL })
      }, n_cores = n_cores, run_multicore = run_multicore, combine = "list", progress = verbose)
    }
    runs[[mdl]] <- wr_ast_standardise(obs, sims, testable, null_smooth, combine, cfg$coco)
    if (mdl == chosen) {
      sims_main <- sims
    }
  }
  main <- runs[[chosen]]
  pkres <- wr_ast_peaks_joint(obs, obs_pk, sims_main, testable, null_smooth, tuned_peaks)

  scores <- data.frame(sedrate = sedrate, untestable = code,
      r2_envelope = obs$r2_env, r2_power = obs$r2_pwr, r2_opt = obs$r2_env * obs$r2_pwr,
      z_envelope = main$z$env, z_power = main$z$pwr, z_opt = main$z$opt, p_opt = main$p$opt,
      B_spectral = obs$B_spec, B_envelope = obs$B_env, B_total = obs$B_spec + obs$B_env,
      z_lik = main$z$lik, p_lik = main$p$lik,
      coco_rho = obs$coco, z_coco = main$z$coco, p_coco = main$p$coco,
      z_comb = main$z$comb, p_comb = main$p$comb, z_cons = main$z$cons, p_cons = main$p$cons,
      S_opt = main$z$sig, p_S = main$p$sig, S_peaks = pkres$S_peaks, p_peaks = pkres$p_peaks,
      S_joint = pkres$S_joint, p_joint = pkres$p_joint)
  bt <- scores$B_total
  bt[!testable] <- NA
  w <- rep(NA_real_, numsed)
  if (any(is.finite(bt))) {
    w <- exp((bt - max(bt, na.rm = TRUE)) / 2)
    w[!is.finite(w)] <- 0
    w <- w / sum(w)
  }
  scores$weight_lik <- w

  best <- wr_ast_best(scores, testable)
  tuning <- NULL
  if (tune == "minimal") {
    # minimal-tuning age model at the rate selected by the noise-corrected TimeOpt score
    sb <- best$best_sedrate[best$score == "TimeOpt significance S_opt"]
    if (length(sb) && is.finite(sb)) {
      tu <- wr_ast_min_tune(y, dx, sb, cfg)
      if (!is.null(tu)) {
        zz <- tu$model$depth
        ti <- tu$model$time
        sr <- c(NA_real_, diff(zz) * 100 / diff(ti))
        sr[1] <- sr[2]
        tuning <- data.frame(depth = dat[1, 1] + zz, time_kyr = ti, sedrate = sr,
            tie = tu$model$tie)
      }
    }
  }
  corroboration <- wr_ast_corroborate(scores, testable)
  sens <- NULL
  for (mdl in names(runs)) {
    sc <- scores
    sc$z_opt <- runs[[mdl]]$z$opt
    sc$p_opt <- runs[[mdl]]$p$opt
    sc$z_lik <- runs[[mdl]]$z$lik
    sc$p_lik <- runs[[mdl]]$p$lik
    sc$z_coco <- runs[[mdl]]$z$coco
    sc$p_coco <- runs[[mdl]]$p$coco
    sc$z_comb <- runs[[mdl]]$z$comb
    sc$p_comb <- runs[[mdl]]$p$comb
    sc$z_cons <- runs[[mdl]]$z$cons
    sc$p_cons <- runs[[mdl]]$p$cons
    sc$S_opt <- runs[[mdl]]$z$sig
    sc$p_S <- runs[[mdl]]$p$sig
    b <- wr_ast_best(sc, testable)
    b <- b[b$score != "r2_opt (raw)", ]
    sens <- rbind(sens, data.frame(noise_model = mdl, AIC = nfit$AIC[nfit$model == mdl], b,
        row.names = NULL))
  }

  tsens <- NULL
  if (target_sensitivity && !is.null(age_ma)) {
    for (lab in c("p_min", "p_nominal", "p_max")) {
      tg2 <- wr_ast_targets(age_ma, targetE, targetP, targetO, which_p = lab)
      cfg2 <- wr_ast_cfg(tg2, fit, obliquity, fitModPwr, ftol, roll, penalty, min_samples, coco,
          band_guard, obliquity_mod, edge_trim)
      cfg2 <- wr_ast_add_tune(cfg2, tg2, tune, tune_width, tune_fit_tuned)
      o2 <- wr_ast_eval(y, dx, sedrate, cfg2)
      r2o <- o2$r2_env * o2$r2_pwr
      r2o[!testable] <- NA
      b2 <- o2$B_spec + o2$B_env
      b2[!testable] <- NA
      cc <- o2$coco
      cc[!testable] <- NA
      tsens <- rbind(tsens, data.frame(precession = lab, p_arcsec = tg2$p_arcsec,
          best_r2_opt = wr_ast_argmax(sedrate, r2o), best_B_total = wr_ast_argmax(sedrate, b2),
          best_coco = wr_ast_argmax(sedrate, cc)))
    }
  }

  null_p <- NULL
  if (keep_null) {
    null_p <- list(S_opt = main$sig$P, S_peaks = pkres$P_peaks, S_joint = pkres$P_joint)
  }
  out <- list(scores = scores, best = best, null_p = null_p, significance_levels = rbind(TimeOpt = main$levels, peaks = pkres$lev_peaks,
      joint = pkres$lev_joint), tuning = tuning, corroboration = corroboration, sensitivity = sens,
      target_sensitivity = tsens,
      noise_fits = nfit[, c("model", "AIC", "description")], noise_fits_rate = nfr,
      noise_model = chosen,
      targets = tg, settings = list(fit = fit, dx = dx, n = n, sedmin = sedmin,
      sedmax = sedmax, numsed = numsed, linLog = linLog, penalty = penalty, combine = combine,
      null_fit = null_fit, mask_pad = mask_pad, tune = tune, tune_width = tune_width,
      n_null = n_null, alpha = alpha, min_samples = min_samples, min_cycles = min_cycles,
      age_ma = age_ma, flow = cfg$flow, fhigh = cfg$fhigh, roll = roll),
      data = dat)
  class(out) <- "astro_sedrate_test"
  out
}

# Internal helpers

# target periods (kyr): AstroGeo22 precession frequency at age_ma (Farhat et al. 2022) and
# present-day La2004 g/s frequencies; which_p selects the nominal, minimum or maximum precession
wr_ast_targets <- function(age_ma = NULL, targetE = NULL, targetP = NULL, targetO = NULL,
    which_p = "p_nominal") {
  if (is.null(age_ma)) {
    p <- 50.467718
  } else {
    a <- astrogeo22(age_ma)
    p <- switch(which_p, p_min = a$p_min_arcsec, p_max = a$p_max_arcsec, a$p_arcsec)
  }
  f <- system.file("extdata", "La04gisi.csv", package = "WaverideR")
  gs <- readLines(f, warn = FALSE, encoding = "UTF-8")
  gs <- gs[!grepl("%", gs) & nzchar(trimws(gs))]
  v <- as.numeric(strsplit(gs[1], ",")[[1]])[1:18]
  g <- v[1:9]
  s <- v[10:18]
  per <- function(fr) {
    360 * 3600 / abs(fr) / 1000
  }
  E <- per(c(g[2] - g[5], g[3] - g[2], g[4] - g[2], g[3] - g[5], g[4] - g[5]))
  P <- per(p + c(g[5], g[2], g[4], g[3]))
  O <- per(p + c(s[3], s[4], s[6]))
  # obliquity-modulation terms s3-s6 (~173 kyr) and s4-s6 (~151 kyr)
  Omod <- per(c(s[3] - s[6], s[4] - s[6]))
  if (!is.null(targetE)) {
    E <- targetE
  }
  if (!is.null(targetP)) {
    P <- targetP
  }
  if (!is.null(targetO)) {
    O <- targetO
  }
  list(E = sort(E, decreasing = TRUE), P = sort(P, decreasing = TRUE),
      O = sort(O, decreasing = TRUE), Omod = sort(Omod, decreasing = TRUE), p_arcsec = p,
      age_ma = age_ma)
}

wr_ast_cfg <- function(tg, fit, obliquity, fitModPwr, ftol, roll, penalty, min_samples, coco,
    band_guard = FALSE, obliquity_mod = FALSE, edge_trim = 0) {
  if (fit == 1) {
    carrier <- tg$P
    modul <- tg$E
    if (is.null(ftol)) {
      ftol <- 0.005
    }
    spec <- tg$P
    if (fitModPwr) {
      spec <- c(tg$E, spec)
    }
  } else {
    carrier <- tg$E[-1]
    modul <- tg$E[1]
    if (is.null(ftol)) {
      ftol <- 0.001
    }
    spec <- tg$E
    if (!fitModPwr) {
      spec <- tg$E[-1]
    }
  }
  if (obliquity) {
    spec <- c(spec, tg$O)
  }
  if (obliquity_mod && length(tg$Omod)) {
    spec <- c(spec, tg$Omod)
  }
  coco_orbit <- c(tg$E, tg$P)
  if (obliquity) {
    coco_orbit <- c(coco_orbit, tg$O)
  }
  coco_target <- NULL
  if (coco) {
    # sinusoid amplitudes as in period2spectrum.m: 1 eccentricity, 0.8 obliquity, 0.6 precession
    wts <- c(rep(1, length(tg$E)), rep(0.6, length(tg$P)))
    if (obliquity) {
      wts <- c(wts, rep(0.8, length(tg$O)))
    }
    ord <- order(coco_orbit, decreasing = TRUE)
    coco_target <- coco_target_from_periods(orbit = coco_orbit[ord], weights = wts[ord])
  }
  flow <- 1 / max(carrier) - ftol
  fhigh <- 1 / min(carrier) + ftol
  if (band_guard) {
    # keep every other astronomical term out of the carrier band: each band edge at most halfway
    # to the nearest term outside the carrier set (e.g. the ~151/173-kyr obliquity modulation next
    # to the short eccentricity)
    others <- setdiff(round(c(tg$E, tg$O, tg$Omod, tg$P), 6), round(carrier, 6))
    longer <- others[others > max(carrier)]
    shorter <- others[others < min(carrier)]
    if (length(longer)) {
      flow <- max(flow, 0.5 * (1 / max(carrier) + 1 / min(longer)))
    }
    if (length(shorter)) {
      fhigh <- min(fhigh, 0.5 * (1 / min(carrier) + 1 / max(shorter)))
    }
  }
  # edge trim of the envelope (kyr at each end): one impulse length of the carrier filter
  trim <- 0
  if (identical(edge_trim, "auto")) {
    trim <- 1 / (fhigh - flow)
  } else if (is.numeric(edge_trim)) {
    trim <- edge_trim
  }
  list(fit = fit, carrier = carrier, modul = modul, spec = sort(unique(spec), decreasing = TRUE),
      flow = flow, fhigh = fhigh, roll = roll, trim = trim,
      penalty = penalty, min_samples = min_samples, coco = coco, coco_target = coco_target)
}

# minimal-tuning settings: fit 2 ties every 405-kyr cycle, fit 1 every short-eccentricity cycle
# (at the harmonic mean of the short-eccentricity periods); the tuned cycle is left out of the
# spectral fit unless tune_fit_tuned = TRUE, because its regularity is imposed by the tuning
wr_ast_add_tune <- function(cfg, tg, tune, tune_width, tune_fit_tuned) {
  cfg$tune <- tune == "minimal"
  cfg$tune_width <- tune_width
  if (cfg$fit == 2) {
    cfg$tune_period <- tg$E[1]
    cfg$tune_band <- c(tg$E[1], tg$E[1])
  } else {
    se <- tg$E[-1]
    cfg$tune_period <- 1 / mean(1 / se)
    cfg$tune_band <- c(max(se), min(se))
  }
  sp <- cfg$spec
  if (!tune_fit_tuned) {
    drop <- sp <= cfg$tune_band[1] * (1 + 1e-9) & sp >= cfg$tune_band[2] * (1 - 1e-9)
    sp <- sp[!drop]
  }
  cfg$spec_tuned <- sp
  cfg
}

# unwrapped phase of the analytic signal (FFT Hilbert transform)
wr_ast_phase <- function(x) {
  n <- length(x)
  h <- numeric(n)
  h[1] <- 1
  if (n %% 2 == 0) {
    h[n / 2 + 1] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[2:((n + 1) / 2)] <- 2
  }
  a <- stats::fft(stats::fft(x) * h, inverse = TRUE) / n
  ph <- Arg(a)
  d <- diff(ph)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(ph[1], ph[1] + cumsum(d))
}

# Minimal tuning at a trial (mean) rate: the tuning cycle is bandpassed in depth around the period
# the rate predicts (+- tune_width), every cycle (each 2 pi of the monotonic Hilbert phase) is tied
# to tune_period, time is interpolated linearly between ties and extrapolated at the trial rate
# beyond the first and last tie. Returns the tuned series resampled at the trial-rate time step,
# or NULL if fewer than three ties are found.
wr_ast_min_tune <- function(y, dx, sedrate, cfg) {
  n <- length(y)
  z <- (0:(n - 1)) * dx
  w <- cfg$tune_width
  fl <- 100 / (sedrate * cfg$tune_band[1] * (1 + w))
  fh <- 100 / (sedrate * cfg$tune_band[2] * max(0.05, 1 - w))
  if (fh >= 0.5 / dx) {
    return(NULL)
  }
  bp <- tryCatch(astrochron::taner(cbind(z, y), padfac = 2, flow = fl, fhigh = fh, roll = 10^3,
      demean = TRUE, detrend = FALSE, addmean = FALSE, output = 1, genplot = FALSE, check = FALSE,
      verbose = FALSE), error = function(e) { NULL })
  if (is.null(bp)) {
    return(NULL)
  }
  ph <- cummax(wr_ast_phase(bp[, 2]))
  k <- ceiling(ph[1] / (2 * pi)):floor(ph[n] / (2 * pi))
  if (length(k) < 3) {
    return(NULL)
  }
  tz <- stats::approx(ph, z, xout = 2 * pi * k, ties = min)$y
  ok <- is.finite(tz) & !duplicated(tz)
  tz <- tz[ok]
  if (length(tz) < 3) {
    return(NULL)
  }
  tt <- (seq_along(tz) - 1) * cfg$tune_period
  tim <- stats::approx(tz, tt, xout = z)$y
  lo <- z < tz[1]
  hi <- z > tz[length(tz)]
  tim[lo] <- tt[1] - (tz[1] - z[lo]) * 100 / sedrate
  tim[hi] <- tt[length(tt)] + (z[hi] - tz[length(tz)]) * 100 / sedrate
  dt <- dx * 100 / sedrate
  tg <- seq(tim[1], tim[n], by = dt)
  if (length(tg) < 10) {
    return(NULL)
  }
  yt <- stats::approx(tim, y, xout = tg)$y
  yt <- yt - mean(yt)
  tie <- rep(FALSE, n)
  for (q in seq_along(tz)) {
    tie[which.min(abs(z - tz[q]))] <- TRUE
  }
  list(y = yt, dt = dt, model = data.frame(depth = z, time = tim - tim[1], tie = tie))
}

# 0 = testable, 1 = too few samples per shortest carrier period, 2 = record shorter than
# min_cycles x longest modulating period
wr_ast_testable <- function(sedrate, dx, n, cfg, min_samples, min_cycles) {
  code <- integer(length(sedrate))
  for (j in seq_along(sedrate)) {
    dt <- dx * 100 / sedrate[j]
    dur <- (n - 1) * dt - 2 * max(0, cfg$trim)
    if (dur < min_cycles * max(cfg$modul)) {
      code[j] <- 2L
    }
    if (min(cfg$carrier) < min_samples * dt) {
      code[j] <- 1L
    }
  }
  code
}

# design matrix (intercept, cos, sin) for periods (kyr) at times tt
wr_ast_design <- function(tt, periods) {
  X <- matrix(1, length(tt), 1 + 2 * length(periods))
  for (i in seq_along(periods)) {
    X[, 2 * i] <- cos(2 * pi * tt / periods[i])
    X[, 2 * i + 1] <- sin(2 * pi * tt / periods[i])
  }
  X
}

# periods that a record of duration dur with sample spacing dt can resolve: drop periods longer
# than the record or with fewer than min_samples samples per cycle, merge periods closer in
# frequency than 1 / dur
wr_ast_resolvable <- function(periods, dur, dt, min_samples) {
  p <- periods[periods <= dur & periods >= min_samples * dt]
  if (length(p) < 2) {
    return(p)
  }
  fr <- sort(1 / p)
  groups <- list(fr[1])
  for (i in 2:length(fr)) {
    last <- groups[[length(groups)]]
    if (fr[i] - max(last) < 1 / dur) {
      groups[[length(groups)]] <- c(last, fr[i])
    } else {
      groups[[length(groups) + 1]] <- fr[i]
    }
  }
  out <- numeric(length(groups))
  for (i in seq_along(groups)) {
    out[i] <- 1 / mean(groups[[i]])
  }
  sort(out, decreasing = TRUE)
}

# conditional least-squares AR(2) innovation variance of a residual series
wr_ast_ar2_var <- function(e) {
  m <- length(e)
  if (m < 10) {
    return(mean(e^2))
  }
  X <- cbind(e[2:(m - 1)], e[1:(m - 2)])
  fitq <- stats::lm.fit(X, e[3:m])
  mean(fitq$residuals^2)
}

wr_ast_tau <- function(e) {
  a <- stats::acf(e, lag.max = length(e) - 2, plot = FALSE, demean = TRUE)$acf[-1]
  z <- which(a <= 0)
  if (length(z)) {
    max(1, z[1])
  } else {
    length(e)
  }
}

# all scores of one (standardised) series at every rate
wr_ast_eval <- function(y, dx, sedrate, cfg) {
  if (isTRUE(cfg$tune)) {
    # minimal tuning at every trial rate, then the scores of the tuned series (time in kyr)
    ns <- length(sedrate)
    out <- list(r2_env = rep(NA_real_, ns), r2_pwr = rep(NA_real_, ns),
        B_spec = rep(NA_real_, ns), B_env = rep(NA_real_, ns), coco = rep(NA_real_, ns))
    cf <- cfg
    cf$tune <- FALSE
    cf$spec <- cfg$spec_tuned
    for (j in seq_len(ns)) {
      tu <- wr_ast_min_tune(y, dx, sedrate[j], cfg)
      if (is.null(tu)) {
        next
      }
      o <- wr_ast_eval(tu$y, tu$dt, 100, cf)
      for (nm in names(out)) {
        out[[nm]][j] <- o[[nm]]
      }
    }
    return(out)
  }
  n <- length(y)
  ns <- length(sedrate)
  r2_env <- rep(NA_real_, ns)
  r2_pwr <- r2_env
  B_spec <- r2_env
  B_env <- r2_env
  pen <- function(k, nn) {
    if (cfg$penalty == "BIC") {
      k * log(nn)
    } else {
      2 * k
    }
  }
  v0 <- wr_ast_ar2_var(y - mean(y))
  for (j in seq_len(ns)) {
    dt <- dx * 100 / sedrate[j]
    tt <- (0:(n - 1)) * dt
    dur <- (n - 1) * dt
    # spectral (power) fit, as TimeOpt
    X <- wr_ast_design(tt, cfg$spec)
    fq <- stats::lm.fit(X, y)
    fitted <- y - fq$residuals
    r2_pwr[j] <- suppressWarnings(stats::cor(y, fitted))^2
    # likelihood with resolvable terms and AR(2) residuals
    pr <- wr_ast_resolvable(cfg$spec, dur, dt, cfg$min_samples)
    if (length(pr)) {
      Xr <- wr_ast_design(tt, pr)
      er <- stats::lm.fit(Xr, y)$residuals
      v1 <- wr_ast_ar2_var(er)
      B_spec[j] <- (n - 2) * log(v0 / v1) - pen(2 * length(pr), n - 2)
    } else {
      B_spec[j] <- 0
    }
    # envelope fit
    bp <- tryCatch(astrochron::taner(cbind(tt, y), padfac = 2, flow = cfg$flow,
        fhigh = cfg$fhigh, roll = cfg$roll, demean = TRUE, detrend = FALSE, addmean = FALSE,
        output = 1, genplot = FALSE, check = FALSE, verbose = FALSE), error = function(e) { NULL })
    if (is.null(bp)) {
      next
    }
    hil <- tryCatch(astrochron::hilbert(bp, padfac = 2, demean = TRUE, detrend = FALSE,
        addmean = FALSE, genplot = FALSE, check = FALSE, verbose = FALSE),
        error = function(e) { NULL })
    if (is.null(hil)) {
      next
    }
    env <- hil[, 2]
    tte <- tt
    if (isTRUE(cfg$trim > 0)) {
      keep <- tt >= cfg$trim & tt <= dur - cfg$trim
      if (sum(keep) < 10) {
        next
      }
      env <- env[keep]
      tte <- tt[keep]
    }
    Xe <- wr_ast_design(tte, cfg$modul)
    fe <- stats::lm.fit(Xe, env)
    r2_env[j] <- suppressWarnings(stats::cor(env, env - fe$residuals))^2
    me <- wr_ast_resolvable(cfg$modul, max(tte) - min(tte), dt, cfg$min_samples)
    e0 <- env - mean(env)
    tau <- wr_ast_tau(e0)
    neff <- max(3, length(env) / tau)
    if (length(me)) {
      e1 <- stats::lm.fit(wr_ast_design(tte, me), env)$residuals
      B_env[j] <- neff * log(mean(e0^2) / mean(e1^2)) - pen(2 * length(me), neff)
    } else {
      B_env[j] <- 0
    }
  }
  coco_rho <- rep(NA_real_, ns)
  if (cfg$coco) {
    cr <- tryCatch(coco_classic(cbind((0:(n - 1)) * dx, y), cfg$coco_target, sr = sedrate,
        nsim = 0, verbose = FALSE)$result$rho, error = function(e) { NULL })
    if (!is.null(cr) && length(cr) == ns) {
      coco_rho <- cr
    }
  }
  list(r2_env = r2_env, r2_pwr = r2_pwr, B_spec = B_spec, B_env = B_env, coco = coco_rho)
}

# Noise models

# Whittle fits of the four noise models to the periodogram of y; keep (logical over the Fourier
# frequencies 1..m) restricts the fit to the frequencies outside the astronomical bands. Without
# keep the AR(2) model is the Yule-Walker estimate; with keep it is a Whittle fit (stationary
# parametrisation through the partial autocorrelations).
wr_ast_noise_fit <- function(y, keep = NULL) {
  n <- length(y)
  m <- floor((n - 1) / 2)
  pg <- Mod(stats::fft(y))[2:(m + 1)]^2 / n
  fr <- (1:m) / n
  if (is.null(keep)) {
    keep <- rep(TRUE, m)
  }
  pk <- pg[keep]
  fk <- fr[keep]
  mk <- length(pk)
  whittle <- function(g) {
    cc <- mean(pk / g)
    -sum(log(cc * g)) - mk
  }
  # AR(1)
  oa <- stats::optimize(function(r) { whittle(wr_ast_shape("ar1", r, fk)) },
      c(0, 0.999), maximum = TRUE)
  # AR(2)
  a2 <- stats::ar.yw(y, aic = FALSE, order.max = 2, demean = TRUE)$ar
  if (mk < m) {
    p2 <- max(-0.98, min(0.98, a2[2]))
    p1 <- max(-0.98, min(0.98, a2[1] / (1 - p2)))
    o2 <- stats::optim(atanh(c(p1, p2)), function(q) {
      pp <- tanh(q)
      -whittle(wr_ast_shape("ar2", c(pp[1] * (1 - pp[2]), pp[2]), fk))
    }, method = "Nelder-Mead")
    pp <- tanh(o2$par)
    a2 <- c(pp[1] * (1 - pp[2]), pp[2])
  }
  ll2 <- whittle(wr_ast_shape("ar2", a2, fk))
  # power law
  op <- stats::optimize(function(b) { whittle(wr_ast_shape("powerlaw", b, fk)) }, c(0, 4),
      maximum = TRUE)
  # bending power law
  ob <- stats::optim(c(log(fr[max(2, round(m / 20))]), 2), function(q) {
    -whittle(wr_ast_shape("bending", c(exp(q[1]), q[2]), fk))
  }, method = "L-BFGS-B", lower = c(log(fr[1]), 0.1), upper = c(log(fr[m]), 6))
  out <- data.frame(model = c("ar1", "ar2", "powerlaw", "bending"),
      loglik = c(oa$objective, ll2, op$objective, -ob$value),
      npar = c(2, 3, 2, 3))
  out$AIC <- -2 * out$loglik + 2 * out$npar
  out$param <- I(list(oa$maximum, a2, op$maximum, c(exp(ob$par[1]), ob$par[2])))
  names(out$param) <- out$model
  # covariance of the parameters on a transformed (unbounded) scale from the Whittle Hessian,
  # used to carry the uncertainty of the noise fit into the noise series
  cv <- vector("list", 4)
  for (q in 1:4) {
    th <- wr_ast_par_to_theta(out$model[q], out$param[[q]])
    nll <- function(v) {
      -whittle(wr_ast_shape(out$model[q], wr_ast_theta_to_par(out$model[q], v), fk))
    }
    H <- tryCatch(stats::optimHess(th, nll), error = function(e) { NULL })
    V <- NULL
    if (!is.null(H) && all(is.finite(H))) {
      V <- tryCatch(solve(H), error = function(e) { NULL })
      if (!is.null(V) && (any(!is.finite(V)) || any(eigen(V, symmetric = TRUE,
          only.values = TRUE)$values <= 0))) {
        V <- NULL
      }
    }
    cv[q] <- list(V)
  }
  out$theta_cov <- I(cv)
  names(out$theta_cov) <- out$model
  out$description <- character(4)
  for (q in 1:4) {
    out$description[q] <- wr_ast_par_text(out$model[q], out$param[[q]])
  }
  out
}

wr_ast_par_to_theta <- function(model, par) {
  if (model == "ar1") {
    atanh(min(par[1], 0.999))
  } else if (model == "ar2") {
    p2 <- max(-0.999, min(0.999, par[2]))
    atanh(c(max(-0.999, min(0.999, par[1] / (1 - p2))), p2))
  } else if (model == "powerlaw") {
    par[1]
  } else {
    c(log(par[1]), par[2])
  }
}

wr_ast_theta_to_par <- function(model, th) {
  if (model == "ar1") {
    tanh(th[1])
  } else if (model == "ar2") {
    pp <- tanh(th)
    c(pp[1] * (1 - pp[2]), pp[2])
  } else if (model == "powerlaw") {
    th[1]
  } else {
    c(exp(th[1]), max(0.05, th[2]))
  }
}

# noise parameters for one noise world: fitted values plus a draw from their uncertainty (u is a
# standard normal vector shared by all rates of the noise world)
wr_ast_draw_par <- function(model, par, V, u) {
  if (is.null(V)) {
    return(par)
  }
  th <- wr_ast_par_to_theta(model, par)
  L <- tryCatch(t(chol(V)), error = function(e) { NULL })
  if (is.null(L)) {
    return(par)
  }
  # step limited to +-1.5 on the transformed scale (poorly constrained fits stay plausible)
  st <- pmax(-1.5, pmin(1.5, as.numeric(L %*% u[seq_along(th)])))
  wr_ast_theta_to_par(model, th + st)
}

wr_ast_par_text <- function(model, par) {
  if (model == "ar1") {
    sprintf("AR(1) rho = %.3f", par[1])
  } else if (model == "ar2") {
    sprintf("AR(2) phi = %.3f, %.3f", par[1], par[2])
  } else if (model == "powerlaw") {
    sprintf("power law beta = %.2f", par[1])
  } else {
    sprintf("bending power law beta = %.2f, f_b = %.4f cycles/sample", par[2], par[1])
  }
}

# spectral shape of a noise model at frequencies f (cycles/sample)
wr_ast_shape <- function(model, par, f) {
  if (model == "ar1") {
    1 / (1 + par[1]^2 - 2 * par[1] * cos(2 * pi * f))
  } else if (model == "ar2") {
    h <- 1 - par[1] * exp(-1i * 2 * pi * f) - par[2] * exp(-2i * 2 * pi * f)
    1 / Mod(h)^2
  } else if (model == "powerlaw") {
    f^(-par[1])
  } else {
    1 / (1 + (f / par[1])^par[2])
  }
}

# Fourier frequencies (1..m of a series of n points with spacing dt kyr) outside the bands of the
# astronomical terms: 405-kyr eccentricity, short eccentricity, obliquity and precession, each
# widened by pad / duration on both sides
wr_ast_band_mask <- function(n, dt, tg, pad = 2, fmax = Inf) {
  m <- floor((n - 1) / 2)
  fk <- (1:m) / (n * dt)
  w <- pad / (n * dt)
  bands <- list(c(1, 1) / tg$E[1])
  if (length(tg$E) > 1) {
    bands[[length(bands) + 1]] <- range(1 / tg$E[-1])
  }
  if (length(tg$O)) {
    bands[[length(bands) + 1]] <- range(1 / tg$O)
  }
  if (length(tg$P)) {
    bands[[length(bands) + 1]] <- range(1 / tg$P)
  }
  keep <- rep(TRUE, m)
  for (b in bands) {
    keep[fk >= b[1] - w & fk <= b[2] + w] <- FALSE
  }
  if (sum(keep) < 20) {
    keep <- rep(TRUE, m)
  }
  # local fit: only frequencies up to fmax (cycles/kyr), extended to at least 20 frequencies
  if (is.finite(fmax)) {
    kl <- keep & fk <= fmax
    if (sum(kl) < 40) {
      kl <- keep & cumsum(keep) <= 40
    }
    keep <- kl
  }
  keep
}

# noise series from given complex normal deviates (common random numbers across rates)
wr_ast_synth <- function(zr, n, model, par) {
  m <- 2 * n
  k <- 1:(m / 2)
  z <- complex(m)
  z[k + 1] <- zr * sqrt(wr_ast_shape(model, par, k / m))
  kk <- k[-length(k)]
  z[m - kk + 1] <- Conj(z[kk + 1])
  z[m / 2 + 1] <- complex(real = Re(z[m / 2 + 1]))
  x <- Re(stats::fft(z, inverse = TRUE))[1:n]
  (x - mean(x)) / stats::sd(x)
}

# scores of one noise world: at every testable rate a series from the noise model fitted at that
# rate, all drawn from the same random deviates
wr_ast_eval_rates <- function(n, dx, sedrate, cfg, testable, model, par_r, pk = NULL) {
  ns <- length(sedrate)
  zr <- complex(real = stats::rnorm(n), imaginary = stats::rnorm(n))
  u <- stats::rnorm(2)
  out <- list(r2_env = rep(NA_real_, ns), r2_pwr = rep(NA_real_, ns), B_spec = rep(NA_real_, ns),
      B_env = rep(NA_real_, ns), coco = rep(NA_real_, ns), peaks = rep(NA_real_, ns))
  for (j in which(testable)) {
    pj <- par_r[[j]]
    if (is.list(pj)) {
      pj <- wr_ast_draw_par(model, pj$par, pj$V, u)
    }
    x <- wr_ast_synth(zr, n, model, pj)
    o <- wr_ast_eval(x, dx, sedrate[j], cfg)
    for (nm in setdiff(names(out), "peaks")) {
      out[[nm]][j] <- o[[nm]]
    }
    if (!is.null(pk)) {
      out$peaks[j] <- wr_ast_tuned_peaks(x, dx, sedrate[j], pk)
    }
  }
  out
}

# noise series of length n with a spectrum shape g(f) (f in cycles/sample), unit variance
wr_ast_spectral_series <- function(n, gfun) {
  m <- 2 * n
  k <- 1:(m / 2)
  amp <- sqrt(gfun(k / m))
  z <- complex(m)
  z[k + 1] <- complex(real = stats::rnorm(m / 2), imaginary = stats::rnorm(m / 2)) * amp
  kk <- k[-length(k)]
  z[m - kk + 1] <- Conj(z[kk + 1])
  z[m / 2 + 1] <- complex(real = Re(z[m / 2 + 1]))
  x <- Re(stats::fft(z, inverse = TRUE))[1:n]
  (x - mean(x)) / stats::sd(x)
}

wr_ast_noise_series <- function(n, model, par) {
  if (model == "ar1") {
    x <- wr_ar1(n, par)
  } else if (model == "ar2") {
    x <- as.numeric(stats::arima.sim(list(ar = par), n = n, n.start = 500))
  } else if (model == "powerlaw") {
    x <- wr_ast_spectral_series(n, function(f) { f^(-par) })
  } else {
    x <- wr_ast_spectral_series(n, function(f) { 1 / (1 + (f / par[1])^par[2]) })
  }
  (x - mean(x)) / stats::sd(x)
}

# Standardisation against the noise series

# z-scores of the observed statistic and leave-one-out z-scores of the noise series
wr_ast_z <- function(v_obs, V, h) {
  ns <- length(v_obs)
  nsim <- ncol(V)
  S <- numeric(ns)
  SS <- numeric(ns)
  cnt <- numeric(ns)
  for (j in seq_len(ns)) {
    v <- V[j, ]
    v <- v[is.finite(v)]
    S[j] <- sum(v)
    SS[j] <- sum(v^2)
    cnt[j] <- length(v)
  }
  ok <- cnt >= 4
  mu <- rep(NA_real_, ns)
  sdv <- mu
  mu[ok] <- S[ok] / cnt[ok]
  sdv[ok] <- sqrt(pmax((SS[ok] - cnt[ok] * mu[ok]^2) / (cnt[ok] - 1), 0))
  mu_s <- wr_smooth_rates(mu, h)
  # floor for the noise standard deviation (10% of its median over the rates), so that rates
  # where the noise statistic hardly varies do not give extreme z-scores
  sd_floor <- max(1e-6, 0.1 * stats::median(sdv[ok], na.rm = TRUE))
  sd_s <- pmax(wr_smooth_rates(sdv, h), sd_floor)
  z_obs <- (v_obs - mu_s) / sd_s
  Z <- matrix(NA_real_, ns, nsim)
  for (k in seq_len(nsim)) {
    x <- V[, k]
    okk <- ok & is.finite(x)
    m_k <- rep(NA_real_, ns)
    s_k <- m_k
    m_k[okk] <- (S[okk] - x[okk]) / (cnt[okk] - 1)
    vk <- (SS[okk] - x[okk]^2 - (cnt[okk] - 1) * m_k[okk]^2) / (cnt[okk] - 2)
    s_k[okk] <- sqrt(pmax(vk, 0))
    Z[, k] <- (x - wr_smooth_rates(m_k, h)) / pmax(wr_smooth_rates(s_k, h), sd_floor)
  }
  list(z = z_obs, Z = Z)
}

# family-wise probability per rate: share of noise series whose maximum over testable rates is
# at least as large
wr_ast_pfw <- function(z, Z, testable) {
  nsim <- ncol(Z)
  mx <- rep(NA_real_, nsim)
  for (k in seq_len(nsim)) {
    v <- Z[testable, k]
    v <- v[is.finite(v)]
    if (length(v)) {
      mx[k] <- max(v)
    }
  }
  mx <- mx[is.finite(mx)]
  p <- rep(NA_real_, length(z))
  for (j in seq_along(z)) {
    if (testable[j] && is.finite(z[j])) {
      p[j] <- (1 + sum(mx >= z[j])) / (length(mx) + 1)
    }
  }
  p
}

wr_ast_lgt <- function(r) {
  r <- pmin(pmax(r, 1e-10), 1 - 1e-10)
  log(r / (1 - r))
}

# upper-tail probability of x in the noise sample 'pool': empirical below the pooled 75 % quantile,
# exponential tail fitted to the exceedances above it (p below 1 / n_null stays resolved)
wr_ast_pw_tail <- function(x, pool, q = 0.75) {
  pool <- sort(pool[is.finite(pool)])
  np <- length(pool)
  p <- rep(NA_real_, length(x))
  if (np < 10) {
    return(p)
  }
  u <- stats::quantile(pool, q, names = FALSE)
  sg <- max(mean(pool[pool > u] - u), 1e-8)
  g <- is.finite(x)
  lo <- g & x <= u
  hi <- g & x > u
  p[lo] <- (1 + np - findInterval(x[lo], pool, left.open = TRUE)) / (1 + np)
  p[hi] <- (1 - q) * exp(-(x[hi] - u) / sg)
  p
}

# significance spectrum: p of the observed statistic at every rate against the noise at that rate
# (rates within h grid steps pooled), and family-wise p by the minimum-p method (Westfall & Young
# 1993): the share of noise series whose smallest p over the testable rates is at least as small
wr_ast_sigspec <- function(obs, V, testable, h) {
  ns <- length(obs)
  nsim <- ncol(V)
  P <- matrix(NA_real_, ns, nsim)
  po <- rep(NA_real_, ns)
  idx <- which(testable)
  for (j in idx) {
    nb <- idx[abs(idx - j) <= max(1, h)]
    pool <- as.numeric(V[nb, ])
    po[j] <- wr_ast_pw_tail(obs[j], pool)
    P[j, ] <- wr_ast_pw_tail(V[j, ], pool)
  }
  minp <- rep(NA_real_, nsim)
  for (k in seq_len(nsim)) {
    v <- P[idx, k]
    v <- v[is.finite(v)]
    if (length(v)) {
      minp[k] <- min(v)
    }
  }
  minp <- minp[is.finite(minp)]
  fw <- rep(NA_real_, ns)
  for (j in idx) {
    if (is.finite(po[j])) {
      fw[j] <- (1 + sum(minp <= po[j])) / (length(minp) + 1)
    }
  }
  lev <- c(p05 = NA_real_, p01 = NA_real_)
  if (length(minp)) {
    lev <- -log10(stats::quantile(minp, c(0.05, 0.01), names = FALSE))
    names(lev) <- c("p05", "p01")
  }
  list(p = po, P = P, fw = fw, levels = lev)
}

# TimeOpt significance spectrum: envelope and power fits each referred to the noise at the same
# rate, combined (Fisher) and calibrated again; S = -log10 of the combined same-rate p
wr_ast_timeopt_sig <- function(obs, sims, testable, h) {
  ns <- length(obs$r2_env)
  nsim <- length(sims)
  Ve <- matrix(NA_real_, ns, nsim)
  Vp <- Ve
  for (k in seq_len(nsim)) {
    Ve[, k] <- wr_ast_lgt(sims[[k]]$r2_env)
    Vp[, k] <- wr_ast_lgt(sims[[k]]$r2_pwr)
  }
  se <- wr_ast_sigspec(wr_ast_lgt(obs$r2_env), Ve, testable, h)
  sp <- wr_ast_sigspec(wr_ast_lgt(obs$r2_pwr), Vp, testable, h)
  cf <- -2 * (log(se$p) + log(sp$p))
  cV <- -2 * (log(se$P) + log(sp$P))
  sc <- wr_ast_sigspec(cf, cV, testable, h)
  S <- -log10(sc$p)
  S[!testable] <- NA
  list(S = S, p = sc$p, P = sc$P, fw = sc$fw, levels = sc$levels)
}

wr_ast_standardise <- function(obs, sims, testable, h, combine, use_coco) {
  ns <- length(obs$r2_env)
  nsim <- length(sims)
  get <- function(name, f) {
    M <- matrix(NA_real_, ns, nsim)
    for (k in seq_len(nsim)) {
      M[, k] <- f(sims[[k]][[name]])
    }
    M
  }
  fisher <- function(r) {
    atanh(pmin(pmax(r, -0.999999), 0.999999))
  }
  ze <- wr_ast_z(wr_logit(obs$r2_env), get("r2_env", wr_logit), h)
  zp <- wr_ast_z(wr_logit(obs$r2_pwr), get("r2_pwr", wr_logit), h)
  lik_obs <- obs$B_spec + obs$B_env
  Lsim <- matrix(NA_real_, ns, nsim)
  for (k in seq_len(nsim)) {
    Lsim[, k] <- sims[[k]]$B_spec + sims[[k]]$B_env
  }
  zl <- wr_ast_z(lik_obs, Lsim, h)
  z_opt <- (ze$z + zp$z) / sqrt(2)
  Z_opt <- (ze$Z + zp$Z) / sqrt(2)
  out_z <- list(env = ze$z, pwr = zp$z, opt = z_opt, lik = zl$z)
  out_p <- list(opt = wr_ast_pfw(z_opt, Z_opt, testable), lik = wr_ast_pfw(zl$z, zl$Z, testable))
  for (nm in c("env", "pwr", "opt", "lik")) {
    out_z[[nm]][!testable] <- NA
  }
  out_z$coco <- rep(NA_real_, ns)
  out_p$coco <- rep(NA_real_, ns)
  C_obs <- z_opt
  C_sim <- Z_opt
  if (use_coco) {
    zc <- wr_ast_z(fisher(obs$coco), get("coco", fisher), h)
    out_z$coco <- zc$z
    out_z$coco[!testable] <- NA
    out_p$coco <- wr_ast_pfw(zc$z, zc$Z, testable)
    C_obs <- C_obs + zc$z
    C_sim <- C_sim + zc$Z
  }
  if (combine == "all") {
    C_obs <- C_obs + zl$z
    C_sim <- C_sim + zl$Z
  }
  zcomb <- wr_ast_z(C_obs, C_sim, h)
  out_z$comb <- zcomb$z
  out_z$comb[!testable] <- NA
  out_p$comb <- wr_ast_pfw(zcomb$z, zcomb$Z, testable)
  # consensus: the smallest of the component z-scores (all methods must support the rate),
  # standardised against the same statistic of the noise series
  out_z$cons <- rep(NA_real_, ns)
  out_p$cons <- rep(NA_real_, ns)
  comp_obs <- list(z_opt)
  comp_sim <- list(Z_opt)
  if (use_coco) {
    comp_obs[[length(comp_obs) + 1]] <- zc$z
    comp_sim[[length(comp_sim) + 1]] <- zc$Z
  }
  if (combine == "all") {
    comp_obs[[length(comp_obs) + 1]] <- zl$z
    comp_sim[[length(comp_sim) + 1]] <- zl$Z
  }
  if (length(comp_obs) > 1) {
    M_obs <- comp_obs[[1]]
    M_sim <- comp_sim[[1]]
    for (q in 2:length(comp_obs)) {
      M_obs <- pmin(M_obs, comp_obs[[q]])
      M_sim <- pmin(M_sim, comp_sim[[q]])
    }
    zcons <- wr_ast_z(M_obs, M_sim, h)
    out_z$cons <- zcons$z
    out_z$cons[!testable] <- NA
    out_p$cons <- wr_ast_pfw(zcons$z, zcons$Z, testable)
  }
  sg <- wr_ast_timeopt_sig(obs, sims, testable, h)
  out_z$sig <- sg$S
  out_p$sig <- sg$fw
  list(z = out_z, p = out_p, levels = sg$levels, sig = sg)
}

# support of the other scores at the rate selected by the noise-corrected TimeOpt score, and the
# ratio of each score's own best rate to that rate (aliases: e.g. ~0.41 obliquity read as short
# eccentricity, ~4.3 eccentricity read as precession)
wr_ast_corroborate <- function(sc, testable) {
  v <- sc$S_opt
  v[!testable] <- NA
  if (!any(is.finite(v))) {
    return(NULL)
  }
  j <- which.max(v)
  out <- NULL
  for (q in list(c("TimeOpt significance", "S_opt", "p_S"), c("TimeOpt z", "z_opt", "p_opt"), c("likelihood", "z_lik", "p_lik"),
      c("COCO", "z_coco", "p_coco"), c("combined", "z_comb", "p_comb"),
      c("consensus", "z_cons", "p_cons"))) {
    zq <- sc[[q[2]]]
    zq[!testable] <- NA
    own <- if (any(is.finite(zq))) { sc$sedrate[which.max(zq)] } else { NA }
    out <- rbind(out, data.frame(score = q[1], sedrate = sc$sedrate[j], z_at_rate = zq[j],
        p_fw_at_rate = sc[[q[3]]][j], own_best = own, ratio_own_best = own / sc$sedrate[j]))
  }
  out
}

wr_ast_argmax <- function(sedrate, v) {
  if (!any(is.finite(v))) {
    return(NA_real_)
  }
  sedrate[which.max(v)]
}

wr_ast_best <- function(sc, testable) {
  rows <- list(c("r2_opt (raw)", "r2_opt", NA), c("TimeOpt significance S_opt", "S_opt", "p_S"),
      c("tuned-band peaks S_peaks", "S_peaks", "p_peaks"), c("joint S_joint", "S_joint", "p_joint"),
      c("TimeOpt z_opt", "z_opt", "p_opt"),
      c("likelihood z_lik", "z_lik", "p_lik"), c("COCO z_coco", "z_coco", "p_coco"),
      c("combined z_comb", "z_comb", "p_comb"), c("consensus z_cons", "z_cons", "p_cons"))
  out <- NULL
  for (r in rows) {
    v <- sc[[r[2]]]
    v[!testable] <- NA
    if (!any(is.finite(v))) {
      out <- rbind(out, data.frame(score = r[1], best_sedrate = NA, value = NA, p_fw = NA))
      next
    }
    j <- which.max(v)
    p <- if (is.na(r[3])) { NA } else { sc[[r[3]]][j] }
    out <- rbind(out, data.frame(score = r[1], best_sedrate = sc$sedrate[j], value = v[j],
        p_fw = p))
  }
  bt <- sc$B_total
  bt[!testable] <- NA
  if (any(is.finite(bt))) {
    j <- which.max(bt)
    out <- rbind(out, data.frame(score = "likelihood B_total (raw)", best_sedrate = sc$sedrate[j],
        value = bt[j], p_fw = NA))
  }
  out
}

#' @export
print.astro_sedrate_test <- function(x, ...) {
  cat("astro_sedrate_test: fit", x$settings$fit, ", noise model", x$noise_model,
      ", n_null =", x$settings$n_null, "\n")
  if (identical(x$settings$null_fit, "masked")) {
    cat("noise fitted at every rate without the astronomical bands (mask_pad =",
        x$settings$mask_pad, ")\n")
  }
  if (!is.null(x$targets$age_ma)) {
    cat("targets at", x$targets$age_ma, "Ma (AstroGeo22 p =", round(x$targets$p_arcsec, 3),
        "arcsec/yr)\n")
  }
  cat("\nbest rates (testable rates only; p_fw = family-wise probability):\n")
  print(x$best, digits = 3, row.names = FALSE)
  if (!is.null(x$corroboration)) {
    cat("\nsupport of the other scores at the TimeOpt rate (z, family-wise p) and their own best rates:\n")
    print(x$corroboration, digits = 3, row.names = FALSE)
  }
  if (!is.null(x$sensitivity)) {
    cat("\nnoise-model sensitivity:\n")
    print(x$sensitivity, digits = 3, row.names = FALSE)
  }
  if (!is.null(x$target_sensitivity)) {
    cat("\nprecession-frequency sensitivity (AstroGeo22 min / nominal / max):\n")
    print(x$target_sensitivity, digits = 4, row.names = FALSE)
  }
  invisible(x)
}

#' @title Plot an astro_sedrate_test result
#' @description Main panel: the significance spectra (-log10 of the same-rate p) of TimeOpt (blue), the
#' tuned-band peaks (green) and the joint score (black), each with its family-wise 5 \% level as a
#' dashed horizontal line; below it the raw TimeOpt r2_opt and the noise-standardised z-scores.
#' Untestable rates are shaded.
#' @param x Result to plot: an object of class \code{astro_sedrate_test}, as returned by
#' \code{\link{astro_sedrate_test}}. Required, no default.
#' @param ... Not used; present for compatibility with the generic \code{plot()}.
#' @return A plot.
#' @export
#' @importFrom graphics par plot lines abline rect legend axis mtext layout points
plot.astro_sedrate_test <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  sc <- x$scores
  sr <- sc$sedrate
  logx <- if (isTRUE(x$settings$linLog == 1)) { "x" } else { "" }
  shade <- function(ylim) {
    un <- sc$untestable != 0
    if (any(un)) {
      ds <- c(diff(sr), diff(sr)[length(sr) - 1])
      for (j in which(un)) {
        graphics::rect(sr[j] - ds[j] / 2, ylim[1], sr[j] + ds[j] / 2, ylim[2],
            col = "grey90", border = NA)
      }
    }
  }
  cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b")
  graphics::layout(matrix(1:3, 3, 1), heights = c(2, 0.9, 1.1))
  ## significance spectra: TimeOpt, tuned-band peaks and joint, each with its family-wise 5 % level
  graphics::par(mar = c(1.5, 4.5, 2.2, 1))
  lev <- x$significance_levels
  if (is.null(dim(lev))) {
    lev <- rbind(TimeOpt = lev)
  }
  curves <- list(S_opt = cols[1], S_peaks = cols[3], S_joint = cols[4])
  labs <- c(S_opt = "TimeOpt", S_peaks = "tuned-band peaks", S_joint = "joint")
  pcol <- c(S_opt = "p_S", S_peaks = "p_peaks", S_joint = "p_joint")
  lrow <- c(S_opt = "TimeOpt", S_peaks = "peaks", S_joint = "joint")
  allv <- 3
  for (nm in names(curves)) {
    if (!is.null(sc[[nm]])) {
      allv <- c(allv, sc[[nm]][is.finite(sc[[nm]])])
    }
  }
  allv <- c(allv, lev[is.finite(lev)])
  yl <- c(0, max(allv) * 1.08)
  plot(sr, sc$S_opt, type = "n", log = logx, ylim = yl, xlab = "", las = 1,
      ylab = expression(-log[10] ~ p ~ "(same rate)"))
  shade(yl)
  shown <- character(0)
  info <- character(0)
  for (nm in names(curves)) {
    S <- sc[[nm]]
    if (is.null(S) || !any(is.finite(S))) {
      next
    }
    if (lrow[nm] %in% rownames(lev) && is.finite(lev[lrow[nm], 1])) {
      graphics::abline(h = lev[lrow[nm], 1], col = curves[[nm]], lty = 2)
    }
    graphics::lines(sr, S, lwd = if (nm == "S_joint") { 2.6 } else { 1.8 }, col = curves[[nm]])
    j <- which.max(S)
    graphics::points(sr[j], S[j], pch = 19, col = curves[[nm]])
    shown <- c(shown, nm)
    info <- c(info, sprintf("%s %.3g (p %.3g)", labs[nm], sr[j], sc[[pcol[nm]]][j]))
  }
  graphics::mtext(paste(info, collapse = " | "), side = 3, line = 0.3, cex = 0.75)
  graphics::legend("topleft", c(labs[shown], "family-wise 5 % (dashed)"),
      col = c(unlist(curves[shown]), "grey40"), lty = c(rep(1, length(shown)), 2), lwd = 2, bty = "n", cex = 0.8)
  ## raw r2_opt
  graphics::par(mar = c(1.5, 4.5, 0.5, 1))
  yl <- c(0, max(sc$r2_opt, na.rm = TRUE))
  plot(sr, sc$r2_opt, type = "n", log = logx, ylim = yl, xlab = "", ylab = "raw r2_opt", las = 1)
  shade(yl)
  graphics::lines(sr, sc$r2_opt, lwd = 1.5, col = "grey40")
  ## z-scores
  graphics::par(mar = c(4.2, 4.5, 0.5, 1))
  zs <- c(sc$z_opt, sc$z_lik, sc$z_coco, sc$z_comb, sc$z_cons)
  yl <- range(c(zs[is.finite(zs)], -1, 3))
  plot(sr, sc$z_opt, type = "n", log = logx, ylim = yl, xlab = "sedimentation rate (cm/kyr)",
      ylab = "z vs noise", las = 1)
  shade(yl)
  graphics::abline(h = 0, col = "grey70")
  graphics::lines(sr, sc$z_opt, lwd = 1.5, col = cols[1])
  graphics::lines(sr, sc$z_lik, lwd = 1.5, col = cols[2])
  if (any(is.finite(sc$z_coco))) {
    graphics::lines(sr, sc$z_coco, lwd = 1.5, col = cols[3])
  }
  graphics::lines(sr, sc$z_comb, lwd = 2, col = cols[4])
  graphics::legend("topleft", c("TimeOpt z", "likelihood", "COCO", "combined"),
      col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.8, horiz = TRUE)
  invisible(NULL)
}

# Tuned-band spectral peaks

# astronomical bands that the minimal tuning does not use: fit 2 tunes the 405-kyr cycle, so short
# eccentricity, obliquity and precession are tested; fit 1 tunes the short eccentricity, so
# obliquity and precession are tested. Members that a periodogram does not separate are grouped
# (main obliquity: the two shortest terms; precession: two pairs; short eccentricity: two pairs).
wr_ast_peak_families <- function(tg, fit) {
  o <- sort(tg$O)
  fam_o <- if (length(o) >= 3) { mean(o[1:2]) } else if (length(o)) { mean(o) } else { numeric(0) }
  p <- sort(tg$P, decreasing = TRUE)
  fam_p <- if (length(p) >= 4) { c(mean(p[1:2]), mean(p[3:4])) } else { mean(p) }
  fam <- c(fam_o, fam_p)
  if (fit == 2) {
    e <- sort(tg$E[-1], decreasing = TRUE)
    fam_e <- if (length(e) >= 4) { c(mean(e[1:2]), mean(e[3:length(e)])) } else { mean(e) }
    fam <- c(fam_e, fam)
  }
  fam
}

# smooth continuum of a periodogram: cubic in log frequency fitted to the log periodogram, with the
# points more than 1.5 log units above the fit left out (five passes), so that spectral peaks do not
# raise the continuum
wr_ast_continuum <- function(v) {
  n <- length(v)
  m <- floor((n - 1) / 2)
  f <- (1:m) / n
  I <- (Mod(stats::fft(v - mean(v)))^2 / n)[2:(m + 1)]
  lf <- log(f)
  li <- log(pmax(I, 1e-300))
  X <- cbind(1, lf, lf^2, lf^3)
  keep <- rep(TRUE, m)
  fit <- rep(mean(li), m)
  for (it in 1:5) {
    if (sum(keep) < 8) {
      break
    }
    cf <- stats::lm.fit(X[keep, , drop = FALSE], li[keep])$coefficients
    cf[!is.finite(cf)] <- 0
    fit <- as.numeric(X %*% cf)
    keep <- (li - fit) < 1.5
  }
  list(f = f, I = I, bg = exp(fit + 0.5772))
}

# whitened peak statistic of a series with spacing dt (kyr): mean over the families of the largest
# log(periodogram / continuum) within +-tol of each family period
wr_ast_comb <- function(v, dt, fam, tol) {
  cb <- wr_ast_continuum(v)
  fz <- cb$f / dt
  w <- log(cb$I / cb$bg)
  acc <- 0
  cnt <- 0
  for (P in fam) {
    f0 <- 1 / P
    k <- abs(log(fz / f0)) < log(1 + tol)
    if (any(k) && f0 > 2 / (length(v) * dt)) {
      acc <- acc + max(w[k])
      cnt <- cnt + 1
    }
  }
  if (cnt < 2) {
    return(NA_real_)
  }
  acc / cnt
}

# minimal tuning at the trial rate, then the whitened peaks of the untuned bands of the tuned series
wr_ast_tuned_peaks <- function(y, dx, sedrate, pk) {
  tu <- wr_ast_min_tune(y, dx, sedrate, pk$cfg)
  if (is.null(tu) || length(tu$y) * tu$dt < 4 * pk$cfg$tune_period) {
    return(NA_real_)
  }
  wr_ast_comb(tu$y, tu$dt, pk$fam, pk$tol)
}

# tuned-band peak score and the joint score (TimeOpt significance and tuned-band peaks combined,
# Fisher, calibrated again on the same noise series), both with minimum-p family-wise probabilities
wr_ast_peaks_joint <- function(obs, obs_pk, sims, testable, h, on) {
  ns <- length(obs_pk)
  na <- rep(NA_real_, ns)
  lev0 <- c(p05 = NA_real_, p01 = NA_real_)
  out <- list(S_peaks = na, p_peaks = na, S_joint = na, p_joint = na, lev_peaks = lev0, lev_joint = lev0, P_joint = NULL,
      P_peaks = NULL)
  if (!on || is.null(sims)) {
    return(out)
  }
  nsim <- length(sims)
  V <- matrix(NA_real_, ns, nsim)
  for (k in seq_len(nsim)) {
    if (!is.null(sims[[k]]$peaks)) {
      V[, k] <- sims[[k]]$peaks
    }
  }
  ok <- testable & is.finite(obs_pk)
  for (j in which(ok)) {
    if (sum(is.finite(V[j, ])) < nsim / 2) {
      ok[j] <- FALSE
    }
  }
  if (!any(ok)) {
    return(out)
  }
  sp <- wr_ast_sigspec(obs_pk, V, ok, h)
  so <- wr_ast_timeopt_sig(obs, sims, testable, h)
  cf <- -2 * (log(so$p) + log(sp$p))
  cV <- -2 * (log(so$P) + log(sp$P))
  sj <- wr_ast_sigspec(cf, cV, ok, h)
  out$S_peaks <- -log10(sp$p)
  out$S_peaks[!ok] <- NA
  out$p_peaks <- sp$fw
  out$P_joint <- sj$P
  out$P_peaks <- sp$P
  out$S_joint <- -log10(sj$p)
  out$S_joint[!ok] <- NA
  out$p_joint <- sj$fw
  out$lev_peaks <- sp$levels
  out$lev_joint <- sj$levels
  out
}
