# Adaptive COCO

Modern Acycle adaptive COCO: nine-term orbital target synthesised on the
finite record at every tested rate, amplitudes fitted to the spectrum,
AR(1) Monte Carlo null with global (max-over-rate) and local p-values.
Translated from Acycle corrcoefslices_rankNew.m (Mingsong Li).

## Usage

``` r
coco_adaptive(
  data,
  orbit9 = orbit9_default(),
  sr = seq(0.5, 20, by = 0.1),
  pad = 5000,
  red = 1,
  nsim = 1000,
  slices = 1,
  method = c("Pearson", "Spearman"),
  max_frequency = NULL,
  amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
  target_model = c("coherent-nine", "phase-averaged"),
  seed = 1,
  batch = NULL,
  verbose = TRUE
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

- sr:

  Sedimentation rates to test, in cm/kyr (numeric vector, \> 0); the
  correlation is computed for every value.
  `Default=seq(0.5, 20, by = 0.1)`.

- pad:

  Length (number of points) to which every slice is zero-padded for the
  periodogram (integer, at least the number of samples in the longest
  slice); a larger value gives a finer frequency grid. `Default=5000`.

- red:

  Background (red-noise) removal applied to every spectrum before the
  correlation (integer): `0` no removal, `1` classical AR(1) background
  subtracted, `2` robust AR(1) background (median-smoothed spectrum,
  Mann and Lees 1996 style), `3` smoothed-window-average (SWA)
  background (needs at least 33 spectral bins). Negative power after
  subtraction is set to zero. `Default=1`.

- nsim:

  Number of AR(1) Monte Carlo simulations used for the global and local
  p-values (integer, \>= 0); `0` skips the simulations, so no p-values
  or pCOCO are returned. `Default=1000`.

- slices:

  Number of equal-depth slices into which the record is cut; each slice
  is standardised and detrended and their spectra are averaged (integer,
  \>= 1; every slice needs at least 4 points). `Default=1`.

- method:

  Correlation coefficient between the data spectrum and the target
  spectrum (character): `"Pearson"` or `"Spearman"`.
  `Default=c("Pearson", "Spearman")`, i.e. `"Pearson"`.

- max_frequency:

  Highest frequency of the target spectrum that enters the correlation,
  in cycles/kyr (numeric, must be at least the highest orbital frequency
  `max(1/orbit9)`); `NULL` uses `1.2 * max(1/orbit9)`. `Default=NULL`.

- amplitude_mode:

  How the amplitudes of the nine target sinusoids are set (character):
  `"four-group-area"` fits one amplitude per group (long eccentricity,
  short eccentricity, obliquity, precession) to the band areas of the
  data spectrum (Acycle default; needs `target_model="coherent-nine"`),
  `"per-orbit"` fits each of the nine amplitudes to the spectral peak at
  its frequency, `"fixed"` uses fixed weights (1.0 eccentricity, 0.8
  obliquity, 0.6 precession; see
  [`coco_fixed_target_weights`](https://stratigraphy.github.io/WaverideR/reference/coco_fixed_target_weights.md)).
  `Default=c("four-group-area", "per-orbit", "fixed")`, i.e.
  `"four-group-area"`.

- target_model:

  How the target spectrum is synthesised on the finite record
  (character): `"coherent-nine"` (the spectrum of the sum of the nine
  sinusoids, including their interference) or `"phase-averaged"` (the
  sum of the individual spectra, averaged over sine and cosine phase).
  `Default=c("coherent-nine", "phase-averaged")`, i.e.
  `"coherent-nine"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=1`.

- batch:

  Number of Monte Carlo simulations processed together in one batch
  (integer, \>= 1); larger batches are faster but use more memory.
  `NULL` chooses it from `pad` (about 4 million spectral values per
  batch). `Default=NULL`.

- verbose:

  Print messages about the preprocessing and the progress of the Monte
  Carlo simulations (logical, TRUE or FALSE). `Default=TRUE`.

## Value

list with \`result\` data.frame: sr, rho, p_param, p_global, p_local,
norbit, pcoco
