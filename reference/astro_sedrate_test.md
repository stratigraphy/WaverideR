# Test sedimentation rates with noise-corrected TimeOpt, likelihood and COCO scores

Scans a grid of constant sedimentation rates for one stratigraphic
record and scores every rate in three ways: the TimeOpt envelope and
spectral fits (Meyers, 2015), a likelihood-based evidence score with an
information-criterion penalty (in the spirit of Malinverno and Meyers,
2024) and the COCO correlation coefficient (Li et al., 2018). Every
score is expressed against noise analysed at the same rate, rates that
cannot be tested are masked, astronomical target periods can be taken at
the age of the record from AstroGeo22 (Farhat et al., 2022), the TimeOpt
and COCO scores are combined into one score, and the result is repeated
for several noise models.

## Usage

``` r
astro_sedrate_test(
  data = NULL,
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
  verbose = FALSE
)
```

## Arguments

- data:

  Evenly spaced stratigraphic record (matrix or data frame with two
  columns: depth in metres in column 1, proxy value in column 2; further
  columns are ignored). The rows are sorted by depth; the record is
  detrended (see `detrend`) and standardised before the analysis.
  Required, no default (`Default=NULL` gives an error).

- sedmin, sedmax:

  Lowest and highest trial sedimentation rate of the rate grid, in
  cm/kyr (numeric, \> 0, `sedmin` \< `sedmax`). `Default=1` and
  `Default=20`.

- numsed:

  Number of trial sedimentation rates in the grid (integer, \>= 1; 1
  tests `sedmin` only). `Default=100`.

- linLog:

  Spacing of the rate grid (numeric): `1` logarithmic, `0` linear.
  `Default=1`.

- fit:

  Which TimeOpt fit is used (numeric): `1` precession carrier modulated
  by eccentricity, `2` short-eccentricity carrier modulated by the
  405-kyr eccentricity. `Default=1`.

- age_ma:

  Age of the record in Ma (numeric, 0 to 4429), used to compute the
  precession and obliquity target periods from the AstroGeo22 precession
  frequency
  ([`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md));
  `NULL` uses the present-day precession frequency. Values given in
  `targetE`, `targetP` or `targetO` override the computed periods.
  `Default=NULL`.

- targetE, targetP, targetO:

  Your own eccentricity, precession and obliquity target periods in kyr
  (numeric vectors, any order; they are sorted from longest to
  shortest). For `fit = 2` the longest value of `targetE` is the
  modulating (405-kyr) period and the others are the carrier periods.
  `NULL` uses the La2004/AstroGeo22 periods (see `age_ma`).
  `Default=NULL`.

- obliquity:

  Include the obliquity periods in the spectral (power) fit and in the
  COCO target (logical, TRUE or FALSE). `Default=TRUE`.

- fitModPwr:

  Include the modulating eccentricity periods in the spectral (power)
  fit (logical, TRUE or FALSE), as in
  [`astrochron::timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html).
  `Default=TRUE`.

- ftol:

  Half-width added outside the carrier frequencies to define the Taner
  bandpass, in cycles/kyr (numeric, \>= 0); `NULL` uses 0.005 for the
  precession carrier (fit 1) and 0.001 for the short-eccentricity
  carrier (fit 2). `Default=NULL`.

- roll:

  Roll-off rate of the Taner bandpass filter used for the carrier band
  (numeric, \> 0; see
  [`astrochron::taner`](https://rdrr.io/pkg/astrochron/man/taner.html)).
  `Default=10^7`.

- detrend:

  Remove a linear trend from the proxy before the analysis (logical,
  TRUE or FALSE). `Default=TRUE`.

- penalty:

  Information-criterion penalty of the likelihood score (character):
  `"BIC"` (k log N) or `"AIC"` (2k). `Default=c("BIC", "AIC")`, i.e.
  `"BIC"`.

- coco:

  Compute the COCO correlation score (and the scores that combine it)
  (logical, TRUE or FALSE). `Default=TRUE`.

- combine:

  Which scores enter the combined scores z_comb and z_cons (character):
  `"timeopt_coco"` (TimeOpt and COCO) or `"all"` (TimeOpt, COCO and the
  likelihood score). `Default=c("timeopt_coco", "all")`, i.e.
  `"timeopt_coco"`.

- noise_model:

  Noise model of the null series (character): `"auto"` (the model with
  the lowest AIC), `"ar1"`, `"ar2"`, `"powerlaw"` or `"bending"`
  (bending power law). With `sensitivity = TRUE` the other models are
  also run, but the main result uses this one.
  `Default=c("auto", "ar1", "ar2", "powerlaw", "bending")`, i.e.
  `"auto"`.

- null_fit:

  How the noise model is fitted to the periodogram (character):
  `"local"` at every rate without the astronomical bands, near the
  tested frequencies and with its parameter uncertainty (see Details);
  `"masked"` at every rate without the bands but with all other
  frequencies; `"full"` once, to the whole periodogram.
  `Default=c("local", "masked", "full")`, i.e. `"local"`.

- mask_pad:

  Widening of every masked astronomical band on each side, in Rayleigh
  frequencies (1 / duration) (numeric, \>= 0); used when `null_fit` is
  `"local"` or `"masked"`. `Default=2`.

- local_factor:

  Upper frequency limit of the `"local"` noise fit, in multiples of the
  highest target frequency of the spectral fit (numeric, \> 1; at least
  40 frequencies are always used). `Default=4`.

- tune:

  Tuning applied at every trial rate before scoring (character):
  `"none"` (constant rate) or `"minimal"` (tune with the 405-kyr cycle
  for fit 2, the short eccentricity for fit 1; see Details).
  `Default=c("none", "minimal")`, i.e. `"none"`.

- tune_width:

  Relative half-width of the bandpass around the period predicted by the
  trial rate that extracts the tuning cycle (numeric, between 0 and 1;
  also the allowed deviation of the local rate); used by
  `tune = "minimal"` and by the tuned-band peak score. `Default=0.3`.

- tune_fit_tuned:

  Keep the tuned cycle in the spectral fit after minimal tuning
  (logical, TRUE or FALSE); `FALSE` leaves it out because the tuning
  imposes its regularity. `Default=FALSE`.

- band_guard:

  Keep every other astronomical term out of the carrier bandpass
  (logical, TRUE or FALSE): each band edge is moved at most halfway to
  the nearest term outside the carrier set (for fit 2 the ~151-kyr s4-s6
  obliquity-modulation term next to the 131-kyr short eccentricity).
  `Default=TRUE`.

- obliquity_mod:

  Include the obliquity-modulation terms s3-s6 (~173 kyr) and s4-s6
  (~151 kyr) in the spectral fit (logical, TRUE or FALSE).
  `Default=TRUE`.

- edge_trim:

  Time in kyr removed at both ends of the envelope before the envelope
  fit (numeric, \>= 0, or the character `"auto"`): `"auto"` removes one
  impulse length of the carrier filter (1 / bandwidth), `0` removes
  nothing (as
  [`astrochron::timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html)).
  `Default="auto"`.

- tuned_peaks:

  Compute the tuned-band peak score S_peaks and the joint score S_joint
  (see Details) (logical, TRUE or FALSE). `Default=TRUE`.

- peak_tol:

  Relative tolerance around each tested period in the tuned-band peak
  score (numeric, \> 0; 0.05 = periods within 5% of the target).
  `Default=0.05`.

- keep_null:

  Return the same-rate p of every noise series (matrix rates x noise
  series) for S_opt, S_peaks and S_joint in `null_p`, as needed by
  [`astro_sedrate_multi`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_multi.md)
  (logical, TRUE or FALSE). `Default=FALSE`.

- sensitivity:

  Repeat the test for all four noise models and report them in
  `sensitivity` (logical, TRUE or FALSE); `FALSE` runs only the chosen
  model, which is faster. `Default=TRUE`.

- target_sensitivity:

  Repeat the scores with the AstroGeo22 minimum and maximum precession
  frequency (logical, TRUE or FALSE; only used when `age_ma` is given).
  `Default=TRUE`.

- n_null:

  Number of noise series simulated per noise model (integer, \> 0); the
  smallest family-wise probability that can be reached is 1/(n_null +
  1). `Default=200`.

- null_smooth:

  Half-width in grid steps of the smoothing of the noise mean and
  standard deviation over neighbouring rates, and of the pooling of the
  noise series of neighbouring rates for S_opt (integer, \>= 0).
  `Default=2`.

- alpha:

  False-alarm level (numeric, between 0 and 1); stored in the `settings`
  of the result. `Default=0.05`.

- min_samples:

  Minimum number of samples per shortest carrier period; rates with
  fewer are marked untestable (code 1) (numeric, \> 0). `Default=3`.

- min_cycles:

  Minimum number of cycles of the longest modulating period that the
  record must span at a rate; rates with fewer are marked untestable
  (code 2) (numeric, \> 0). `Default=2`.

- run_multicore:

  Simulate and score the noise series in parallel on several cores
  (logical, TRUE or FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used for the noise series (integer); `NULL` uses the
  option `waverider.n_cores` or, with `run_multicore = TRUE`, all but
  one detected core (at most `getOption("waverider.max_cores", 4)`).
  `Default=NULL`.

- verbose:

  Print progress messages (logical, TRUE or FALSE). `Default=FALSE`.

## Value

A list of class `astro_sedrate_test` with

- scores:

  data frame per rate: sedrate, untestable code, r2_envelope, r2_power,
  r2_opt, z_envelope, z_power, z_opt, p_opt, B_spectral, B_envelope,
  B_total, weight_lik, z_lik, p_lik, coco_rho, z_coco, p_coco, z_comb,
  p_comb, z_cons, p_cons, S_opt, p_S, S_peaks, p_peaks, S_joint, p_joint
  (p = family-wise probability; S = -log10 of the same-rate p)

- best:

  data frame: best rate, value and family-wise p of each score (raw
  r2_opt included)

- null_p:

  with `keep_null = TRUE`: same-rate p of every noise series (rates x
  noise series) for S_opt, S_peaks and S_joint; otherwise NULL

- significance_levels:

  matrix (rows TimeOpt, peaks, joint): S values of the family-wise 5 %
  and 1 % levels

- tuning:

  for `tune = "minimal"`: depth, time (kyr), sedimentation rate (cm/kyr)
  and tie points of the minimal-tuning age model at the rate selected by
  S_opt

- corroboration:

  value and family-wise p of every score at the rate selected by S_opt,
  and the ratio of each score's own best rate to that rate (to spot
  aliases)

- sensitivity:

  best rate, z and p of every score for every noise model

- target_sensitivity:

  best rates with the AstroGeo22 minimum / nominal / maximum precession
  frequency

- noise_fits:

  fitted noise models with their AIC (summed over the testable rates for
  `null_fit = "local"` or `"masked"`)

- noise_fits_rate:

  for `null_fit = "local"` or `"masked"`: AIC and parameters of every
  noise model at every testable rate, and the fraction of periodogram
  frequencies kept

- targets:

  the target periods used

- settings:

  the settings

## Details

**Why a same-rate noise correction.** The raw TimeOpt r-squared values
(and COCO correlations) that pure noise reaches depend on the trial
rate: at high rates the record spans a short time, fewer cycles fit, and
the regression has fewer independent data points per fitted sinusoid;
for red noise the target frequencies also move to the strong
low-frequency end of the spectrum. Li et al. (2018) describe the same
effect for COCO (fewer astronomical terms resolved at high rates). Here
every statistic is compared with the same statistic of noise series
analysed at the same rate: \$\$z(u) = (S(u) - \mu\_{noise}(u)) /
\sigma\_{noise}(u)\$\$ (r-squared values on a logit scale, COCO
correlations on a Fisher z scale; the noise mean and standard deviation
are smoothed over neighbouring rates). The family-wise probability of a
score is the fraction of noise series whose maximum z over all testable
rates is at least as large (as in the full-grid test of Meyers and
Malinverno, 2018), with a floor of 1/(n_null + 1).

**TimeOpt scores.** For each rate the record is calibrated to time, the
carrier band (precession for fit 1, short eccentricity for fit 2) is
bandpassed with
[`taner`](https://rdrr.io/pkg/astrochron/man/taner.html), its amplitude
envelope is taken with
[`hilbert`](https://rdrr.io/pkg/astrochron/man/hilbert.html) and
regressed on the modulating periods (eccentricity for fit 1, the 405-kyr
term for fit 2), and the record itself is regressed on the target
periods. r2_envelope, r2_power and r2_opt follow Meyers (2015); z_opt
combines the z-scores of the envelope and power fits (Stouffer's
method).

**Significance spectrum (main TimeOpt score).** S_opt refers the
envelope and the power fit each to the noise distribution at the same
rate (the noise series of that rate and of the `null_smooth`
neighbouring rates pooled; an exponential tail above the pooled 75 %
quantile resolves p below 1 / n_null), combines the two p-values
(Fisher) and calibrates the combination the same way; S_opt = -log10 of
that p. The family-wise p (p_S) uses the minimum-p method (Westfall and
Young, 1993): the share of noise series whose smallest p over the
testable rates is at least as small. Because every rate is referred to
its own noise distribution, noise gives a flat baseline near zero over
the whole range (no drift to fast rates as in raw r2_opt, and no
inflated low-rate scatter as in z_opt), and the family-wise levels are
horizontal lines.

**Tuned-band peaks and joint score.** TimeOpt's envelope fit has little
power when the record spans few modulating cycles and the rate changes
(in a Kosov-like test, about 2 Myr with +-17 % rate changes, the fit-2
envelope reached p \< 0.05 at the true rate in 15-30 % of records with a
known signal). With `tuned_peaks = TRUE` a second score is computed at
every rate: the record is minimally tuned on the 405-kyr cycle (fit 2)
or on the short eccentricity (fit 1) at that rate (as
`tune = "minimal"`), and the bands not used for the tuning (short
eccentricity, obliquity, precession for fit 2; obliquity and precession
for fit 1) are tested for spectral peaks: the mean over these bands of
the largest log ratio of the periodogram to a smooth continuum (cubic in
log frequency fitted to the log periodogram with the peaks left out)
within `peak_tol` of each period. Every noise series is tuned and scored
in exactly the same way, so the regularity that the tuning imposes is
part of the null and the test is not circular; the tuned band itself is
not scored. S_peaks is -log10 of the same-rate p with the minimum-p
family-wise p_peaks, as for S_opt. In the Kosov-like test this score
reached p \< 0.05 in 90-100 % of records with a known signal. S_joint
combines S_opt and S_peaks (Fisher) and is calibrated again on the same
noise series (p_joint); use it as the main summary when both are
available.

**Likelihood score.** This score is complementary to TimeOpt and follows
ideas of [`timeOptB`](https://rdrr.io/pkg/astrochron/man/timeOptB.html)
(Malinverno and Meyers, 2024): the residuals of the spectral fit are
modelled as an AR(2) process (the likelihood is that of the AR(2)
innovations), and the envelope fit uses an effective number of
independent points N_eff = N / tau, with tau the first lag at which the
autocorrelation of the envelope reaches zero. At every rate the
astronomical model is compared with a model without astronomical terms,
and an information-criterion penalty is applied for the number of
resolvable terms: \$\$B = 2 (\log L\_{astro} - \log L\_{noise}) - k \log
N\$\$ (BIC, `penalty = "BIC"`) or \\- 2k\\ (AIC). Only target terms that
the record can resolve at that rate enter the likelihood model: terms
with fewer than `min_samples` samples per cycle or with periods longer
than the record are dropped, and terms closer in frequency than the
Rayleigh resolution (1 / duration) are merged into one. B \> 0 means the
astronomical model is preferred at that rate; `weight_lik` gives rate
weights proportional to exp(B / 2). The likelihood score is a model
comparison per rate, not a full Bayesian posterior: the envelope differs
between rates, so the weights are indicative. Use
[`timeOptB`](https://rdrr.io/pkg/astrochron/man/timeOptB.html) for
posterior distributions of the sedimentation rate and the precession
frequency.

**COCO.** The correlation coefficient between the AR(1)-corrected
periodogram of the record and an astronomical target periodogram
([`coco_classic`](https://stratigraphy.github.io/WaverideR/reference/coco_classic.md),
[`coco_target_from_periods`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_periods.md);
Li et al., 2018).

**Combined scores.** z_comb is the sum of z_opt and z_coco, standardised
again against the same noise series, so the correlation between the
TimeOpt and COCO scores is accounted for. z_cons (consensus) is the
smaller of the two z-scores, standardised in the same way: it is only
high where both methods support the rate, which suppresses peaks that
only one method shows (e.g. the COCO alias at ~4.3 times the true rate,
where eccentricity cycles are read as precession and short
eccentricity). With `combine = "all"` the likelihood score is added to
both.

**Noise models.** AR(1), AR(2), power law (f^-beta) and bending power
law (1 / (1 + (f / f_b)^beta); Vaughan et al., 2011), fitted to the
periodogram of the record (Whittle likelihood). `noise_model = "auto"`
uses the model with the lowest AIC. `sensitivity = TRUE` repeats the
test for every noise model.

**Noise fitted without the cycles** (`null_fit`). A noise model fitted
to the whole periodogram of the record also fits the astronomical cycles
in it, and its shape near the tested frequencies is set mainly by the
many high-frequency periodogram values. Both distort the background at
the frequencies that matter: cycles that are present bend the fitted
spectrum (e.g. an AR(1) background is fitted as AR(2)), so the null
partly contains the cycles it is tested against and genuine cycles look
less significant. Therefore, by default (`null_fit = "local"`), the
noise model is fitted separately for every trial rate (Whittle
likelihood):

- leaving out the frequency bands where that rate puts the astronomical
  terms (405-kyr eccentricity, short eccentricity, obliquity and
  precession, each band widened by `mask_pad` Rayleigh frequencies 1 /
  duration), as in robust background estimates that exclude spectral
  peaks (Mann and Lees, 1996; Meyers, 2012);

- using only frequencies up to `local_factor` times the highest target
  frequency of the spectral fit (at least 40 frequencies), so the
  background is estimated where the tested cycles are;

- carrying the uncertainty of the fitted parameters (Whittle Hessian)
  into the noise series: every noise series draws its parameters from
  that uncertainty, so a background that is poorly constrained by few
  frequencies gives a wider null instead of false alarms.

The noise series for a rate come from the model fitted at that rate; the
same random numbers are used at all rates, so each noise series is one
coherent "noise world" over the rate grid and the family-wise
probability keeps its meaning. The model choice (`noise_model = "auto"`)
uses the AIC summed over the testable rates. `null_fit = "masked"`
leaves out the bands but uses all other frequencies; `null_fit = "full"`
fits one noise model to the whole periodogram (earlier behaviour).

**Minimal tuning** (`tune = "minimal"`). A constant trial rate smears
the cycles when the true rate varies, so the TimeOpt fit can favour a
wrong rate. With minimal tuning the record is first tuned at every trial
(mean) rate with the highest-order cycle only: for fit 2 the 405-kyr
cycle, for fit 1 the short eccentricity (tied at the harmonic mean of
its periods, ~110 kyr). The cycle is bandpassed in depth around the
period the trial rate predicts (+- `tune_width`), every cycle of its
(monotonic) Hilbert phase is tied to the cycle period, time is
interpolated linearly between the ties and extrapolated at the trial
rate beyond them, and TimeOpt, the likelihood score and COCO are
computed on the tuned series. The tuned cycle is left out of the
spectral fit (`tune_fit_tuned = FALSE`) because the tuning imposes its
regularity. Every noise series is tuned in exactly the same way at every
rate, so the gain in fit that tuning gives to any series is part of the
null and the probabilities stay calibrated. `tuning` returns the
minimal-tuning age model at the rate selected by z_opt. The trial rate
is then the mean rate of a tuned age model, so the test asks at which
scale a 405-kyr (or short-eccentricity) tuning gives the best
astronomical fit.

**Envelope and band corrections (fit 2).** Three settings protect the
modulation fit of TimeOpt (fit 2 in particular) against artefacts found
on real records:

- `edge_trim`: the Taner bandpass and the Hilbert transform ring at the
  ends of a finite record, and the ringing extends over about one
  impulse length of the filter (1 / bandwidth, ~230 kyr for the
  short-eccentricity band). With only a few 405-kyr cycles these end
  values steer the envelope regression; they are removed before the fit
  (and from the duration used for the untestable mask).

- `band_guard`: the default fit-2 band (1/130.7 - 0.001 to 1/94.9 +
  0.001 cycles/kyr, i.e. 87-151 kyr) contains the s4-s6
  obliquity-modulation term (151 kyr); the long-period edge is moved
  halfway to it (~140 kyr).

- `obliquity_mod`: the s3-s6 and s4-s6 terms are part of the spectral
  model, so their power counts at the rate where they belong. Otherwise
  a rate ~1.2-1.3 times too high, at which the short-eccentricity
  targets fall on the obliquity-modulation cycles, can win.

With `band_guard = FALSE, obliquity_mod = FALSE, edge_trim = 0` the
r-squared values are those of astrochron::timeOpt.

**Age-dependent targets.** With `age_ma`, the target periods are
computed from the AstroGeo22 precession frequency (Farhat et al., 2022;
[`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md))
and the La2004 secular frequencies (Laskar et al., 2004; held at their
present values): eccentricity g2-g5, g3-g2, g4-g2, g3-g5, g4-g5;
precession k+g5, k+g2, k+g4, k+g3; obliquity k+s3, k+s4, k+s6.
`target_sensitivity = TRUE` repeats the scores with the AstroGeo22
minimum and maximum precession frequency.

**Untestable rates** (`untestable`): 1 = fewer than `min_samples`
samples per shortest carrier period; 2 = the record spans fewer than
`min_cycles` of the longest modulating period.

## References

Meyers, S.R. (2015). The evaluation of eccentricity-related amplitude
modulation and bundling in paleoclimate data: An inverse approach for
astrochronologic testing and time scale optimization. Paleoceanography
30, 1625-1640 \<doi:10.1002/2015PA002850\>

Westfall, P.H. and Young, S.S. (1993). Resampling-Based Multiple
Testing: Examples and Methods for p-Value Adjustment. Wiley, New York,
ISBN 978-0-471-55761-6.

Meyers, S.R. and Malinverno, A. (2018). Proterozoic Milankovitch cycles
and the history of the solar system. PNAS 115, 6363-6368
\<doi:10.1073/pnas.1717689115\>

Malinverno, A. and Meyers, S.R. (2024). Bayesian estimation of past
astronomical frequencies, lunar distance, and length of day from
sediment cycles. Geochemistry, Geophysics, Geosystems 25, e2023GC011176
\<doi:10.1029/2023GC011176\>

Li, M., Kump, L.R., Hinnov, L.A. and Mann, M.E. (2018). Tracking
variable sedimentation rates and astronomical forcing in Phanerozoic
paleoclimate proxy series with evolutionary correlation coefficients and
hypothesis testing. Earth and Planetary Science Letters 501, 165-179
\<doi:10.1016/j.epsl.2018.08.041\>

Farhat, M., Auclair-Desrotour, P., Boue, G. and Laskar, J. (2022). The
resonant tidal evolution of the Earth-Moon distance. Astronomy and
Astrophysics 665, L1 \<doi:10.1051/0004-6361/202243445\>

Mann, M.E. and Lees, J.M. (1996). Robust estimation of background noise
and signal detection in climatic time series. Climatic Change 33,
409-445 \<doi:10.1007/BF00142586\>

Meyers, S.R. (2012). Seeing red in cyclic stratigraphy: Spectral noise
estimation for astrochronology. Paleoceanography 27, PA3228
\<doi:10.1029/2012PA002307\>

Vaughan, S., Bailey, R.J. and Smith, D.G. (2011). Detecting cycles in
stratigraphic data: Spectral analysis in the presence of red noise.
Paleoceanography 26, PA4211 \<doi:10.1029/2011PA002195\>

## Examples

``` r
# \donttest{
set.seed(1)
res <- astro_sedrate_test(data = mag[1:400, ], sedmin = 0.3, sedmax = 1.5, numsed = 25,
    fit = 1, n_null = 10, sensitivity = FALSE, target_sensitivity = FALSE)
res$best
#>                         score best_sedrate         value      p_fw
#> 1                r2_opt (raw)     1.311728    0.04348622        NA
#> 2  TimeOpt significance S_opt     1.500000    0.84551274 0.7272727
#> 3    tuned-band peaks S_peaks     1.003110    1.61256883 0.2727273
#> 4               joint S_joint     1.003110    0.96998545 0.6363636
#> 5               TimeOpt z_opt     1.311728    0.71897364 0.8181818
#> 6            likelihood z_lik     1.500000    1.89492606 0.6363636
#> 7                 COCO z_coco     1.500000   -0.01031668 0.8181818
#> 8             combined z_comb     1.500000    0.21500952 0.8181818
#> 9            consensus z_cons     1.500000    0.34280317 0.6363636
#> 10   likelihood B_total (raw)     1.500000 -112.30576931        NA
# }
```
