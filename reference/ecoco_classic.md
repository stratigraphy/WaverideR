# Evolutionary COCO (Li et al. 2018)

COCO in a sliding window (2019 Acycle ecoco.m driver) with best-rate
summary and ridge tracking per window.

## Usage

``` r
ecoco_classic(
  data,
  target,
  orbit = NULL,
  window,
  step = 1,
  sr = seq(0.2, 20, by = 0.2),
  dt = NULL,
  pad = 5000,
  red = 1,
  nsim = 0,
  slices = 1,
  adjust = 0,
  delinear = 1,
  method = c("Pearson", "Spearman"),
  mc = c("first-window", "each-window"),
  seed = NULL,
  track_ridge = TRUE,
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
  [`coco_target_from_periods`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_periods.md)
  or
  [`coco_target_from_series`](https://stratigraphy.github.io/WaverideR/reference/coco_target_from_series.md);
  passed to
  [`coco_classic`](https://stratigraphy.github.io/WaverideR/reference/coco_classic.md)
  for every window. Required, no default.

- orbit:

  Astronomical periods in kyr (numeric vector), used to count how many
  astronomical parameters contribute at each rate (and, with `adjust=1`,
  to define the target bands; then seven periods are needed). `NULL`
  uses the periods stored with the target (attribute `"orbit"`), or else
  the seven periods of Li et al. (2018) (400, 120.5, 93.1, 38.6, 23.0,
  21.8 and 18.5 kyr). `Default=NULL`.

- window:

  Length of the sliding window, in metres (numeric, \> 0; it must
  contain at least 3 samples and at most the whole record). Required, no
  default.

- step:

  Step between successive windows, in number of samples (integer, \>= 1
  and smaller than half the number of samples). `Default=1`.

- sr:

  Sedimentation rates to test in every window, in cm/kyr (numeric
  vector, \> 0). `Default=seq(0.2, 20, by = 0.2)`.

- dt:

  Sampling interval of `data`, in metres (numeric, \> 0); `NULL` uses
  the median spacing of the depth column. It also sets the number of
  samples per window, `trunc(window / dt)`. `Default=NULL`.

- pad:

  Length (number of points) to which each window is zero-padded for the
  periodogram (integer, at least the number of samples per window); a
  larger value gives a finer frequency grid. `Default=5000`.

- red:

  Red-noise treatment of the window periodograms, passed to
  [`coco_classic`](https://stratigraphy.github.io/WaverideR/reference/coco_classic.md)
  (integer): `0` none; `1` the AR(1) background is subtracted (as in Li
  et al. 2018); `2` the periodogram is divided by the AR(1) background
  minus 1; `3` a robust AR(1) background (median-smoothed spectrum) is
  subtracted; a value from `50` to `99` divides by that AR(1) confidence
  level (in %) minus 1. Negative values are set to zero. `Default=1`.

- nsim:

  Number of Monte Carlo simulations of the null hypothesis of no
  astronomical forcing, used for the p-value `p_h0` and pCOCO (integer,
  \>= 0; 2000 in Li et al. 2018); `0` skips them and the score becomes
  `norbit * rho`. See `mc` for where the null is simulated. `Default=0`.

- slices:

  Number of equal-depth slices of each window whose periodograms are
  averaged (integer, \>= 1); `1` uses the whole window. `Default=1`.

- adjust:

  Adjust the power of the target in four bands (long eccentricity, short
  eccentricity, obliquity, precession) to the power of the data in those
  bands at each rate (integer): `0` no, `1` yes (only for a seven-period
  target). `Default=0`.

- delinear:

  Remove the mean of each window before the analysis, as in Acycle
  ecoco.m (numeric): `1` yes, any other value no. `Default=1`.

- method:

  Correlation coefficient between the data and target spectra
  (character): `"Pearson"` (as in Li et al. 2018) or `"Spearman"`.
  `Default=c("Pearson", "Spearman")`, i.e. `"Pearson"`.

- mc:

  Where the Monte Carlo null distribution is simulated (character):
  `"first-window"` reproduces Acycle 2018-2019, the null is simulated
  once in the first window and reused for all windows; `"each-window"`
  simulates a separate null for every window (slower). Only used when
  `nsim > 0`. `Default=c("first-window", "each-window")`, i.e.
  `"first-window"`.

- seed:

  Seed of the random number generator (integer), so that the Monte Carlo
  results can be reproduced; `NULL` leaves the current random state
  unchanged. `Default=NULL`.

- track_ridge:

  Track a smooth sedimentation-rate ridge through the score map with
  [`ecoco_track_ridge`](https://stratigraphy.github.io/WaverideR/reference/ecoco_track_ridge.md)
  (logical, TRUE or FALSE); `FALSE` returns the best rate per window
  ([`ecoco_best_per_window`](https://stratigraphy.github.io/WaverideR/reference/ecoco_best_per_window.md))
  as `track`. `Default=TRUE`.

- verbose:

  Print progress messages for every window (logical, TRUE or FALSE).
  `Default=TRUE`.

## Value

list of matrices (rows = sr, columns = windows): rho, p_h0, norbit,
pcoco, score (= norbit/norb \* pcoco) plus \`track\` (best rate per
window)
