# Adaptive eCOCO

Adaptive COCO in physical-depth sliding windows with per-window AR(1)
Monte Carlo nulls and ridge tracking (Acycle ecoco.m adaptive branch,
ecocoAdaptiveCore.m; Mingsong Li).

## Usage

``` r
ecoco_adaptive(
  data,
  orbit9 = orbit9_default(),
  window,
  sr = seq(0.5, 20, by = 0.1),
  step = 1,
  step_depth = NULL,
  pad = 5000,
  red = 1,
  nsim = 500,
  method = c("Pearson", "Spearman"),
  max_frequency = NULL,
  amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
  seed = 1,
  center_limits = NULL,
  batch = NULL,
  verbose = TRUE,
  centre_limits = NULL
)
```

## Arguments

- data:

  Record to analyse (matrix or data frame, two columns): depth in metres
  in column 1 and the proxy in column 2. Non-finite rows are dropped,
  duplicated depths averaged, and unevenly spaced data are linearly
  interpolated to the median sample spacing (see
  [`coco_prepare_regular_data`](https://stratigraphy.github.io/WaverideR/reference/coco_prepare_regular_data.md)).
  Required, no default.

- orbit9:

  The nine astronomical target periods in kyr (numeric vector of nine
  distinct values: long eccentricity, four short eccentricity terms,
  obliquity, three precession terms), e.g. from
  [`orbit9_at_age`](https://stratigraphy.github.io/WaverideR/reference/orbit9_at_age.md)
  or
  [`calculate_orbit9`](https://stratigraphy.github.io/WaverideR/reference/calculate_orbit9.md)
  for the age of the record. The default uses the present-day periods
  (405.7, 130.7, 123.9, 98.9, 94.9, 41.0, 23.7, 22.4 and 19.0 kyr).
  `Default=orbit9_default()`.

- window:

  Length of each sliding window, in metres (numeric, \> 0, shorter than
  the record; every window must contain at least 4 data points). Each
  window is resampled to an odd number of points around its centre and
  analysed with Adaptive COCO. Required, no default.

- sr:

  Sedimentation rates to test in every window, in cm/kyr (numeric
  vector, \> 0). `Default=seq(0.5, 20, by = 0.1)`.

- step:

  Step between successive window centres, in samples of the
  (regularised) record (numeric, \> 0); only used when
  `step_depth=NULL`, giving `step_depth = step * sample interval`.
  `Default=1`.

- step_depth:

  Exact spacing of the window centres, in metres (numeric, \> 0, or
  `NULL`); overrides `step`. At most 10000 windows are allowed.
  `Default=NULL`.

- pad:

  Length (number of points) to which every window is zero-padded for the
  periodogram (integer, at least the number of points per window); a
  larger value gives a finer frequency grid. `Default=5000`.

- red:

  Background (red-noise) removal applied to every window spectrum before
  the correlation (integer): `0` no removal, `1` classical AR(1)
  background subtracted, `2` robust AR(1) background (median-smoothed
  spectrum, Mann and Lees 1996 style), `3` smoothed-window-average (SWA)
  background (needs at least 33 spectral bins). Negative power after
  subtraction is set to zero. `Default=1`.

- nsim:

  Number of AR(1) Monte Carlo simulations per window used for the global
  and local p-values (integer, \>= 0); the AR(1) coefficient and
  standard deviation are estimated in each window. `0` skips the
  simulations, so the p-values are `NA`. `Default=500`.

- method:

  Correlation coefficient between the window spectrum and the target
  spectrum (character): `"Pearson"` or `"Spearman"`.
  `Default=c("Pearson", "Spearman")`, i.e. `"Pearson"`.

- max_frequency:

  Highest frequency of the target spectrum that enters the correlation,
  in cycles/kyr (numeric, must be at least the highest orbital frequency
  `max(1/orbit9)`); `NULL` uses `1.2 * max(1/orbit9)`. Passed to
  [`coco_adaptive_evaluate`](https://stratigraphy.github.io/WaverideR/reference/coco_adaptive_evaluate.md).
  `Default=NULL`.

- amplitude_mode:

  How the amplitudes of the nine target sinusoids are set (character):
  `"four-group-area"` fits one amplitude per group (long eccentricity,
  short eccentricity, obliquity, precession) to the band areas of the
  data spectrum (Acycle default), `"per-orbit"` fits each of the nine
  amplitudes to the spectral peak at its frequency, `"fixed"` uses fixed
  weights (1.0 eccentricity, 0.8 obliquity, 0.6 precession; see
  [`coco_fixed_target_weights`](https://stratigraphy.github.io/WaverideR/reference/coco_fixed_target_weights.md)).
  The target is always synthesised as `"coherent-nine"`.
  `Default=c("four-group-area", "per-orbit", "fixed")`, i.e.
  `"four-group-area"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=1`.

- center_limits:

  Depths of the first and last window centre, in metres (numeric vector
  `c(first, last)`, or `NULL`); the windows must still fit inside the
  record. `NULL` uses
  `c(first depth + window/2, last depth - window/2)`. `Default=NULL`.

- batch:

  Number of Monte Carlo simulations processed together in one batch
  (integer, \>= 1); larger batches are faster but use more memory.
  `NULL` chooses `min(50, nsim, 3e6 / (frequencies * windows))`.
  `Default=NULL`.

- verbose:

  Print messages about the windows and the progress of the Monte Carlo
  simulations (logical, TRUE or FALSE). `Default=TRUE`.

- centre_limits:

  British spelling of `center_limits`; if given (not `NULL`) it
  overrides `center_limits`. `Default=NULL`.

## Value

list with matrices (rate-by-window): rho, p_global, p_local, norbit,
pcoco, score; \`track\` = ridge-tracked sedimentation-rate path
