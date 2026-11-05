# calculate the duration of stratigraphic gaps using astronomical cycles

calculate the duration of stratigraphic gaps using the duration of
stable astronomical cycles

## Usage

``` r
dur_gaps(
  proxies = NULL,
  retracked_period_1 = NULL,
  retracked_period_2 = NULL,
  min_max = NULL,
  n_simulations = 10,
  tracked_cycle_period = NULL,
  tracked_cycle_period_unc = NULL,
  tracked_cycle_period_unc_dist = "u",
  pts = 5,
  dj = 1/200,
  lowerPeriod = 1,
  upperPeriod = 3200,
  period_max = NULL,
  period_min = NULL,
  missing_cycle_dur = NULL,
  n_cycles_missing = 1,
  missing_cycle_unc = NULL,
  missing_cycle_unc_dist = "u",
  seed_nr = 1337,
  run_multicore = FALSE,
  n_cores = NULL
)
```

## Arguments

- proxies:

  Proxy records that were used to build the astrochronological age model
  and from which the gap duration is estimated (list of two-column
  matrices or data frames: depth/height in column 1, proxy value in
  column 2). In each simulation one proxy of the list is drawn at
  random. Required (`Default=NULL`).

- retracked_period_1:

  Tracked cycle curve of the interval below/before the gap (three-column
  matrix or data frame: depth/height in column 1, period of the tracked
  cycle in column 2, its 1 standard deviation uncertainty in column 3,
  e.g. the output of
  [`retrack_wt_MC`](https://stratigraphy.github.io/WaverideR/reference/retrack_wt_MC.md)).
  In each simulation a curve is drawn from this normal uncertainty
  envelope and used to tune the proxy with
  [`curve2tune`](https://stratigraphy.github.io/WaverideR/reference/curve2tune.md).
  The gap must lie between `retracked_period_1` and
  `retracked_period_2`. Required (`Default=NULL`).

- retracked_period_2:

  Tracked cycle curve of the interval above/after the gap, in the same
  format as `retracked_period_1`. Required (`Default=NULL`).

- min_max:

  For each proxy in `proxies`, whether the gap duration is measured
  between minima (`"min"`, peaks found with
  [`min_detect`](https://stratigraphy.github.io/WaverideR/reference/min_detect.md))
  or maxima (`"max"`,
  [`max_detect`](https://stratigraphy.github.io/WaverideR/reference/max_detect.md))
  of the extracted cycle (list or character vector of the same length as
  `proxies`). Required (`Default=NULL`).

- n_simulations:

  Number of Monte Carlo simulations, i.e. number of gap durations
  calculated (integer, \> 0). `Default=10`.

- tracked_cycle_period:

  Duration in time of the tracked cycle (numeric, e.g. in kyr; e.g. 405
  for the long eccentricity cycle). The tuned records and the gap
  duration are in this time unit. Required (`Default=NULL`).

- tracked_cycle_period_unc:

  Uncertainty of `tracked_cycle_period`, in the same unit (numeric, \>=
  0): half-width of the range for `"u"`, standard deviation for `"n"`
  (see `tracked_cycle_period_unc_dist`). Required (`Default=NULL`).

- tracked_cycle_period_unc_dist:

  Distribution from which the tracked cycle period is drawn in each
  simulation (character): `"u"` uniform between
  `tracked_cycle_period -/+ tracked_cycle_period_unc` or `"n"` normal
  with mean `tracked_cycle_period` and standard deviation
  `tracked_cycle_period_unc`. `Default="u"`.

- pts:

  Number of points on each side that a point must exceed (maximum) or be
  below (minimum) to be detected as a peak in the extracted cycle;
  passed to
  [`max_detect`](https://stratigraphy.github.io/WaverideR/reference/max_detect.md)/[`min_detect`](https://stratigraphy.github.io/WaverideR/reference/min_detect.md)
  (integer, \>= 1). Larger values give more certain peaks but may miss
  minor ones. `Default=5`.

- dj:

  Spacing between successive scales of the wavelet transform of the
  tuned records, in octaves (numeric, \> 0), passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  Smaller values give a finer period resolution but a longer computation
  time. `Default=1/200`.

- lowerPeriod:

  Shortest period of the wavelet transform of the tuned records, in the
  time unit of `tracked_cycle_period` (numeric, \> 0), passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  `Default=1`.

- upperPeriod:

  Longest period of the wavelet transform of the tuned records, in the
  time unit of `tracked_cycle_period` (numeric, \> `lowerPeriod`),
  passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  `Default=3200`.

- period_max:

  Upper period boundary (time unit) of the cycle extracted from the
  tuned records with
  [`extract_signal_stable_V2`](https://stratigraphy.github.io/WaverideR/reference/extract_signal_stable_V2.md)
  (numeric, \> `period_min`). Required (`Default=NULL`).

- period_min:

  Lower period boundary (time unit) of the cycle extracted from the
  tuned records with
  [`extract_signal_stable_V2`](https://stratigraphy.github.io/WaverideR/reference/extract_signal_stable_V2.md)
  (numeric, \> 0). Required (`Default=NULL`).

- missing_cycle_dur:

  Duration of one missing cycle in the gap, in the time unit (numeric).
  If it equals `tracked_cycle_period`, the value drawn for the tracked
  cycle is also used for the missing cycle. Required (`Default=NULL`).

- n_cycles_missing:

  Number of cycles of duration `missing_cycle_dur` that are missing in
  the gap (numeric, \>= 0). `Default=1`.

- missing_cycle_unc:

  Uncertainty of `missing_cycle_dur` in the time unit (numeric, \>= 0):
  half-width for `"u"`, standard deviation for `"n"`. Only used when
  `missing_cycle_dur` differs from `tracked_cycle_period`. Required in
  that case (`Default=NULL`).

- missing_cycle_unc_dist:

  Distribution from which the missing cycle duration is drawn
  (character): `"u"` uniform or `"n"` normal (see `missing_cycle_unc`).
  `Default="u"`.

- seed_nr:

  Seed of the random number generator (integer), set with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) so that the Monte
  Carlo results can be reproduced. `Default=1337`.

- run_multicore:

  Run the simulations in parallel on several cores (logical, TRUE or
  FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore=TRUE` (integer or `NULL`).
  `NULL` picks an automatic, safe number: physical cores minus one, at
  most 4, and limited by the free memory; set
  `options(waverider.n_cores = ...)` to change that default. Results do
  not depend on the number of cores. `Default=NULL`.

## Value

a vector with all the calculated gap durations
