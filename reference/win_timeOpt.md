# Windowed timeOpt sedimentation rate estimation

The `win_timeOpt` function conducts a windowed timeOpt sedimentation
rate estimation. This function is based on
[`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html) and
[`eTimeOpt`](https://rdrr.io/pkg/astrochron/man/eTimeOpt.html) of the
'astrochron' package, but allows for multithreaded analysis, windows of
fixed duration, a correction of the r-squared values for their bias
towards high sedimentation rates, a mask of rates that cannot be tested
and an obliquity version of the fit.

## Usage

``` r
win_timeOpt(
  data = NULL,
  window_size = 10,
  sedmin = 0.5,
  sedmax = 2,
  numsed = 100,
  limit = FALSE,
  fit = 2,
  fitModPwr = TRUE,
  flow = NULL,
  fhigh = NULL,
  roll = 10^6,
  targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
  targetP = c(20.9, 19.9, 17.1, 17.2),
  targetO = c(53.6, 41, 39.6),
  targetOmod = 173,
  detrend = TRUE,
  normalize = TRUE,
  linLog = 1,
  window_mode = c("depth", "time"),
  window_duration = 1500,
  step = 1,
  centres = NULL,
  bias_correction = c("null", "none"),
  noise_model = c("auto", "ar1", "powerlaw"),
  n_null = 50,
  null_rho_step = 0.05,
  null_beta_step = 0.1,
  null_smooth = 2,
  alpha = 0.05,
  min_samples = 3,
  min_cycles = 2,
  null_cache = NULL,
  run_multicore = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  normalise = NULL
)
```

## Arguments

- data:

  Record to analyse (matrix or data frame with two columns: depth in
  metres in column 1, proxy value in column 2), evenly sampled (the
  depth steps must be equal within 0.1 percent; use e.g. `linterp` of
  'astrochron' to resample). Must be given (the `NULL` default fails).
  `Default=NULL`.

- window_size:

  Thickness of the moving window in metres (numeric, \> 0), only used
  with `window_mode = "depth"`. `Default=10`.

- sedmin:

  Lowest sedimentation rate of the grid (numeric, \> 0, in cm/kyr).
  `Default=0.5`.

- sedmax:

  Highest sedimentation rate of the grid (numeric, \> `sedmin`, in
  cm/kyr). `Default=2`.

- numsed:

  Number of sedimentation rates of the grid between `sedmin` and
  `sedmax` (integer, \>= 1; 1 tests `sedmin` only). `Default=100`.

- limit:

  Limit the evaluated sedimentation rates to the range in which the full
  target signal can be recovered (logical, TRUE or FALSE; passed to
  [`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html)). Only
  used with `window_mode = "depth"`; the untestable mask (`min_samples`,
  `min_cycles`) is usually the better choice. `Default=FALSE`.

- fit:

  Which amplitude modulation is tested (numeric, one of 1, 2, 3): 1
  precession modulated by eccentricity, 2 short eccentricity modulated
  by the 405-kyr cycle, 3 obliquity modulated by the obliquity
  modulation periods (`targetO`, `targetOmod`; see Details).
  `Default=2`.

- fitModPwr:

  Include the modulation periods in the spectral power fit (logical,
  TRUE or FALSE; passed to
  [`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html); not used
  for `fit = 3`). `Default=TRUE`.

- flow:

  Low-frequency cut-off of the Taner bandpass filter that isolates the
  carrier (numeric, half power point in cycles/kyr). `NULL` uses the
  default of
  [`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html) for
  `fit = 1` and `fit = 2`, and 1/46 for `fit = 3`. `Default=NULL`.

- fhigh:

  High-frequency cut-off of the Taner bandpass filter (numeric, half
  power point in cycles/kyr). `NULL` uses the default of
  [`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html) for
  `fit = 1` and `fit = 2`, and 1/36 for `fit = 3`. `Default=NULL`.

- roll:

  Roll-off rate of the Taner filter (numeric, \> 0, in dB/octave);
  larger values give steeper filter edges. `Default=10^6`.

- targetE:

  Eccentricity periods to evaluate (numeric vector, in kyr), in order of
  decreasing period with 405 kyr (or 405.7) first.
  `Default=c(405.7, 130.7, 123.8, 98.9, 94.9)`.

- targetP:

  Precession periods to evaluate (numeric vector, in kyr), in order of
  decreasing period; used for `fit = 1` and in the spectral fit.
  `Default=c(20.9, 19.9, 17.1, 17.2)`.

- targetO:

  Obliquity carrier periods (numeric vector, in kyr), only used for
  `fit = 3`. `Default=c(53.6, 41.0, 39.6)`.

- targetOmod:

  Obliquity modulation periods (numeric vector, in kyr) used in the
  envelope fit of `fit = 3`. `Default=173`.

- detrend:

  Remove a linear trend from every window before the analysis (logical,
  TRUE or FALSE); also used for the noise model fit. `Default=TRUE`.

- normalize:

  Divide the raw r-squared curve of every window by its maximum, so that
  it peaks at 1 (logical, TRUE or FALSE). The bias-corrected z-scores
  are not normalised. `Default=TRUE`.

- linLog:

  Spacing of the sedimentation rate grid (numeric, 0 or 1): 0 linear, 1
  logarithmic. `Default=1`.

- window_mode:

  How the window is defined (character, one of `"depth"`, `"time"`):
  `"depth"` a window of fixed thickness `window_size` (the behaviour of
  earlier versions); `"time"` a window of fixed duration
  `window_duration` for every trial sedimentation rate, which removes
  most of the bias towards high rates (see Details). `Default="depth"`.

- window_duration:

  Duration of the window in kyr (numeric, \> 0), only used with
  `window_mode = "time"`; its thickness at a trial rate is
  `window_duration * sedrate / 100` m. `Default=1500`.

- step:

  Spacing of the window centres in samples: every `step`-th data point
  is a window centre (integer, \>= 1); larger values are faster. Ignored
  when `centres` is given. `Default=1`.

- centres:

  Depths of the window centres to evaluate (numeric vector, in metres;
  each is moved to the nearest sample); overrides `step`. `NULL` uses
  `step`. `Default=NULL`.

- bias_correction:

  Correction of the r-squared values for their bias towards high
  sedimentation rates (character, one of `"null"`, `"none"`): `"null"`
  standardises them against noise series of a fitted noise model and
  gives z-scores and family-wise probabilities (see Details); `"none"`
  only returns the raw r-squared values (much faster). `Default="null"`.

- noise_model:

  Noise model of the correction (character, one of `"auto"`, `"ar1"`,
  `"powerlaw"`): `"ar1"` AR(1) red noise, `"powerlaw"` \\f^{-\beta}\\
  noise, `"auto"` the one of the two that fits the periodogram of the
  window best (Whittle likelihood). Only used with
  `bias_correction = "null"`. `Default="auto"`.

- n_null:

  Number of noise series per grid node of the noise model (integer, \>=
  1); more series give more stable z-scores and probabilities and take
  longer (the family-wise probabilities need at least 5). `Default=50`.

- null_rho_step:

  Step of the grid of AR(1) coefficients on which the noise level is
  computed and interpolated (numeric, between 0 and 1); smaller steps
  are more precise and slower. `Default=0.05`.

- null_beta_step:

  Step of the grid of power-law exponents on which the noise level is
  computed and interpolated (numeric, \> 0). `Default=0.1`.

- null_smooth:

  Half width (in grid steps of the sedimentation rate, integer \>= 0) of
  the moving average of the noise mean and standard deviation over
  neighbouring rates; 0 means no smoothing. `Default=2`.

- alpha:

  False-alarm level of the family-wise threshold `z_threshold` and of
  the column `significant` (numeric, between 0 and 1). `Default=0.05`.

- min_samples:

  Minimum number of samples per shortest carrier cycle (precession for
  fit 1, short eccentricity for fit 2, obliquity for fit 3) for a rate
  to be testable (numeric, \> 0; in tests timeOpt still worked with 3.4
  samples per shortest precession cycle). `Default=3`.

- min_cycles:

  Minimum number of cycles of the longest fitted period that the window
  must span for a rate to be testable (numeric, \> 0). `Default=2`.

- null_cache:

  Noise tables to reuse: the `null_cache` element of an earlier result
  of `win_timeOpt` with the same settings (list); missing grid nodes are
  added. `NULL` computes all tables. `Default=NULL`.

- run_multicore:

  Run the windows and the noise simulations in parallel on several cores
  (logical, TRUE or FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer; `NULL`: an
  automatic, safe number: physical cores minus one, at most 4, and
  limited by the free memory). Set `options(waverider.n_cores = ...)` to
  change the default. A value \> 1 also switches on parallel computing.
  Results do not depend on the number of cores. `Default=NULL`.

- verbose:

  Show progress bars of the window and noise computations (logical, TRUE
  or FALSE). `Default=FALSE`.

- normalise:

  British spelling of `normalize`; if given (not `NULL`) it overrides
  `normalize`. `Default=NULL`.

## Value

Returns a list of class `win_timeOpt.result` with the elements

- r_2_envelope, r_2_power, r_2_opt:

  matrices (sedimentation rate x window) of the raw r-squared values
  (normalised when `normalize = TRUE`)

- r_2_envelope_avg, r_2_power_avg, r_2_opt_avg:

  their averages over all windows

- sedrate:

  the sedimentation rates (cm/ka)

- x, y:

  depth and proxy value of the window centres

- linLog:

  the linLog value

- z_envelope, z_power, z_opt:

  z-scores against noise (`NULL` when `bias_correction = "none"`; `NA`
  for untestable cells)

- z_envelope_avg, z_power_avg, z_opt_avg:

  their averages over all windows

- p_fw:

  family-wise probability that noise reaches z_opt anywhere on the rate
  grid

- z_threshold:

  family-wise 1 - alpha threshold of z_opt per window

- best_sedrate, best_z, p_best, significant:

  per window: the rate with the highest z_opt among the testable rates
  (raw r2_opt when there is no correction), its z-score and family-wise
  probability, and whether `p_best < alpha`

- untestable_code, testable:

  see Details

- noise_family, noise_param:

  fitted noise model per cell (1 = AR(1) coefficient, 2 = power-law
  exponent)

- npts:

  number of points of each window

- data:

  the input data

- window_mode, window_size, window_duration, fit, bias_correction,
  noise_model, alpha:

  the settings

- null_cache:

  the noise tables, for reuse

## Details

**Bias towards high sedimentation rates.** In a window of fixed
thickness, a higher trial sedimentation rate means a shorter time span,
so fewer cycles of the target periods fit in the window. The envelope
fit and the spectral power fit then have fewer independent data points
per fitted sinusoid, and for red noise the target frequencies move to
the strong low-frequency end of the spectrum. Both raise the r-squared
values that pure noise reaches, so the raw r-squared surface has a
rising baseline, and a genuine peak at a low sedimentation rate can be
smaller than the noise level at high sedimentation rates.

**Windows of fixed duration** (`window_mode = "time"`). For every trial
sedimentation rate the window is cut to `window_duration` kyr (its
thickness is `window_duration * sedrate / 100` m), so every rate is
tested on the same number of cycles. This removes most of the bias at
its source. Rates whose window does not fit in the record around a given
depth are marked as untestable there. timeOpt is run once per rate. With
`window_mode = "depth"` (default) the window has a fixed thickness
`window_size` (the behaviour of earlier versions).

**Noise-level correction** (`bias_correction = "null"`). For each window
a noise model is fitted to the periodogram (Whittle likelihood; for time
windows on the window of the middle rate that fits in the record): an
AR(1) spectrum (flat at low frequencies, bending to f^-2) and a power
law (f^-beta, e.g. pink noise for beta = 1). With `noise_model = "auto"`
the better-fitting one is used. Noise series of the same model, length
and sample spacing are analysed with exactly the same timeOpt settings.
This gives, for every trial sedimentation rate, the mean and standard
deviation of the noise r-squared values (on a logit scale, because
r-squared is bounded by 0 and 1). These are smoothed over neighbouring
rates (`null_smooth`) and the observed values are expressed as z-scores:
\$\$z\_{env} = (logit(r^2\_{env}) - \mu\_{env}) / \sigma\_{env}\$\$ and
likewise for the spectral power fit; both are combined (Stouffer's
method) as \\z\_{opt} = (z\_{env} + z\_{pwr}) / \sqrt{2}\\. A z-score of
0 is the noise level at that sedimentation rate. To keep the computation
affordable, the noise level is computed on a grid of model parameters
(steps `null_rho_step` and `null_beta_step`; for depth windows also
three window lengths) and interpolated. The noise tables can be reused
with `null_cache`.

**Probabilities.** For each grid node the maximum z_opt over all
testable rates is recorded for every noise series. This gives a
family-wise false-alarm threshold per window (`z_threshold`, the 1 -
`alpha` quantile) and a family-wise probability per cell (`p_fw`): the
probability that noise reaches this z-score anywhere on the rate grid.

**Untestable rates** (`untestable_code`): 1 = fewer than `min_samples`
samples per shortest carrier cycle (precession for fit 1, short
eccentricity for fit 2, obliquity for fit 3); 2 = the window spans fewer
than `min_cycles` of the longest fitted period; 3 = the window does not
fit in the record. Untestable cells get no z-score.

**Obliquity version** (`fit = 3`): timeOpt fit 1 with obliquity carrier
periods (`targetO`) bandpassed between `flow` and `fhigh` (default
1/46-1/36 cycles/kyr) and the obliquity modulation periods
(`targetOmod`) in the envelope fit. A precession record read as
obliquity gives a peak at half the true rate, an obliquity record read
as precession (fit 1) a peak at twice the true rate; compare the fits
with
[`win_timeOpt_multi`](https://stratigraphy.github.io/WaverideR/reference/win_timeOpt_multi.md).

## References

Meyers, S.R. (2019). Cyclostratigraphy and the problem of
astrochronologic testing. Earth-Science Reviews 190, 190-223
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Based on the
[`timeOpt`](https://rdrr.io/pkg/astrochron/man/timeOpt.html) and
[`eTimeOpt`](https://rdrr.io/pkg/astrochron/man/eTimeOpt.html) functions
of the 'astrochron' R package.

## Examples

``` r
# \donttest{
#Conduct a windowed timeOpt on the magnetic susceptibility record
#of the Sullivan core of Pas et al., (2018).
set.seed(1)
mag_win_timeOpt <- win_timeOpt(
data = mag[1:250, ],
window_size = 6,
sedmin = 0.1,
sedmax = 1,
numsed = 5,
limit = FALSE,
fit = 2,
fitModPwr = TRUE,
flow = NULL,
fhigh = NULL,
roll = 10 ^ 6,
targetE = c(405.7, 130.7, 123.8, 98.9, 94.9),
targetP = c(20.9, 19.9, 17.1, 17.2),
detrend = TRUE,
normalize = TRUE,
linLog = 1,
bias_correction = "null",
n_null = 6,
noise_model = "ar1",
step = 10,
run_multicore = FALSE,
verbose = FALSE)
# }
```
