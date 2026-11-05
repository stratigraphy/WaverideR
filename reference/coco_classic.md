# COCO analysis (Li et al. 2018)

Correlation coefficient between the (AR(1)-corrected) data periodogram
and an astronomical target over a grid of sedimentation rates, with the
H0 (no orbital forcing) Monte Carlo test of Li et al. (2018). Translated
from Acycle corrcoefslices.m / cyclecorr.m (Mingsong Li).

## Usage

``` r
coco_classic(
  data,
  target,
  orbit = NULL,
  sr = seq(0.2, 20, by = 0.2),
  dt = NULL,
  pad = 5000,
  red = 1,
  nsim = 0,
  slices = 1,
  adjust = 0,
  method = c("Pearson", "Spearman"),
  seed = NULL,
  keep_null = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  Record to analyse (matrix or data frame; the first two columns are
  used): depth in metres in column 1 and the proxy in column 2, evenly
  spaced (it is not interpolated; incomplete rows are dropped and the
  rows sorted by depth). Required, no default.

- target:

  Astronomical target periodogram (two-column matrix or data frame):
  frequency in cycles/kyr in column 1 and power in column 2, e.g. the
  output of
  [`coco_target_from_periods`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_periods.md),
  [`coco_target_from_series`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_series.md).
  Its peaks higher than half the mean power are rebuilt as Gaussian
  peaks with the Rayleigh width of the data at every tested rate.
  Required, no default.

- orbit:

  Astronomical periods in kyr (numeric vector), used to count how many
  astronomical parameters contribute at each rate (and, with `adjust=1`,
  to define the target bands; then seven periods are needed). `NULL`
  uses the periods stored with the target (attribute `"orbit"`, set by
  [`coco_target_from_periods`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_periods.md)),
  or else the seven periods of Li et al. (2018) (400, 120.5, 93.1, 38.6,
  23.0, 21.8 and 18.5 kyr). `Default=NULL`.

- sr:

  Sedimentation rates to test, in cm/kyr (numeric vector, \> 0).
  `Default=seq(0.2, 20, by = 0.2)`.

- dt:

  Sampling interval of `data`, in metres (numeric, \> 0); `NULL` uses
  the median spacing of the depth column. `Default=NULL`.

- pad:

  Length (number of points) to which the data (or each slice) is
  zero-padded for the periodogram (integer, at least the number of
  samples); a larger value gives a finer frequency grid. `Default=5000`.

- red:

  Red-noise treatment of the data periodogram (integer): `0` none; `1`
  the AR(1) background is subtracted (as in Li et al. 2018); `2` the
  periodogram is divided by the AR(1) background minus 1; `3` a robust
  AR(1) background (median-smoothed spectrum) is subtracted; a value
  from `50` to `99` divides by that AR(1) confidence level (in %) minus
  1, so only power above it remains. Negative values are set to zero.
  `Default=1`.

- nsim:

  Number of Monte Carlo simulations of the null hypothesis of no
  astronomical forcing (random spectra with as many peaks as the data
  spectrum), used for the p-value `p_h0` and pCOCO (integer, \>= 0; 2000
  in Li et al. 2018); `0` skips them. `Default=0`.

- slices:

  Number of equal-depth slices whose periodograms are averaged (integer,
  \>= 1); `1` uses the whole record. `Default=1`.

- adjust:

  Adjust the power of the target in four bands (long eccentricity, short
  eccentricity, obliquity, precession) to the power of the data in those
  bands at each rate (integer): `0` no, `1` yes (only for a seven-period
  target). `Default=0`.

- method:

  Correlation coefficient between the data and target spectra
  (character): `"Pearson"` (as in Li et al. 2018) or `"Spearman"`.
  `Default=c("Pearson", "Spearman")`, i.e. `"Pearson"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=NULL`.

- keep_null:

  Return the Monte Carlo null correlations (matrix, rates in rows,
  simulations in columns) as element `null` of the output (logical, TRUE
  or FALSE). `Default=FALSE`.

- verbose:

  Print progress messages of the Monte Carlo simulations (logical, TRUE
  or FALSE). `Default=TRUE`.

## Value

list with data.frame \`result\` (sr, rho, p_param, p_h0, norbit, pcoco)
