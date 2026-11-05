# Interleaved eCOCO (odd/even cross-validated evolutionary COCO)

Interleaved eCOCO evaluates complete sliding windows of the cleaned
(finite, sorted, duplicate-averaged) record. Within each window the
globally odd observations train a frozen four-group coherent-nine
astronomical target that is validated on the globally even observations,
and vice versa (Interleaved cvCOCO applied per window). Each fold is
linearly regularised at its own median spacing after the split. The
consensus correlation is the same-rate minimum of the two directional
curves; its global p-value is a plus-one Monte Carlo p against the null
maximum over the complete rate grid. Each null realisation is one
stationary Gaussian AR(1) process on the complete window order, split
and interpolated like the data. Windows that cannot be resolved for the
rate grid are returned empty. pCOCO = rho \* \|log10(p_global)\| is the
ridge score. Translated from Acycle v3.1 ecocoInterleavedCore.m,
interleavedcvcoco.m and cvcoco.m (Mingsong Li). Requires the cocoR
helpers (part of this package).

## Usage

``` r
ecoco_interleaved(
  data = NULL,
  orbit9 = c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.682, 22.3758,
    18.9519),
  window = NULL,
  sr = NULL,
  step_depth = NULL,
  step = 1,
  pad = NULL,
  red = 0,
  n_sim = 2000,
  method = "Pearson",
  max_frequency = NULL,
  seed = 1,
  window_mode = "physical-depth",
  dt = NULL,
  batch_size = 20,
  genplot = FALSE,
  verbose = FALSE
)
```

## Arguments

- data:

  Record to analyse (matrix or data frame, two columns): depth in metres
  in column 1 and the proxy in column 2 (at least 16 unique depths).
  Give the raw, possibly irregular data; do not interpolate beforehand,
  because the odd and even observations are split before each fold is
  regularised. Non-finite rows are dropped, rows sorted and duplicated
  depths averaged. Required (`Default=NULL`).

- orbit9:

  The nine astronomical target periods in kyr (numeric vector of nine
  distinct positive values: long eccentricity, four short eccentricity
  terms, obliquity, three precession terms), e.g. from
  [`orbit9_at_age`](https://stratigraphy.github.io/WaverideR/reference/orbit9_at_age.md)
  for the age of the record.
  `Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)`.

- window:

  Length of the sliding window, in metres (numeric, \> 0, at most the
  record span; windows with fewer than 16 observations are left empty).
  `NULL` uses 25% of the record length. `Default=NULL`.

- sr:

  Sedimentation rates to test, in cm/kyr (numeric vector, \> 0, strictly
  increasing with a constant step). `NULL` uses the Acycle GUI default:
  from `2 * dt * 0.065 * 100` up to 100 cm/kyr (reduced so that the
  longest period fits in the record), in steps of 0.1 cm/kyr or at most
  300 steps. `Default=NULL`.

- step_depth:

  Spacing of the window centres, in metres (numeric, \> 0); only used
  with `window_mode="physical-depth"`. `NULL` uses `dt`, or
  `dt * ceiling(n / 300)` for more than 300 samples. `Default=NULL`.

- step:

  Step between successive windows, in samples (integer, \>= 1); only
  used with `window_mode="legacy-count"`. `Default=1`.

- pad:

  Length (number of points) to which the periodogram of each fold is
  zero-padded; each fold uses `max(pad, fold length)` (integer, \> 0).
  `NULL` uses the Acycle default: 5000 for up to 2500 samples, 10000 for
  up to 5000 samples, and `floor(n / 5000) * 5000 + 5000` above that.
  `Default=NULL`.

- red:

  Background (red-noise) removal applied to every spectrum (integer):
  `0` none, `1` classical AR(1) background, `2` robust AR(1) background,
  `3` smoothed-window-average (SWA) background. `Default=0`.

- n_sim:

  Number of AR(1) Monte Carlo simulations per window used for the global
  and local p-values (integer, \>= 0). `Default=2000`.

- method:

  Correlation coefficient between the data and target spectra
  (character): `"Pearson"` or `"Spearman"`. `Default="Pearson"`.

- max_frequency:

  Highest frequency of the target spectrum that enters the correlation,
  in cycles/kyr (numeric, at least `max(1/orbit9)`); `NULL` uses
  `1.2 * max(1/orbit9)`. `Default=NULL`.

- seed:

  Parent seed of the random number generator (integer). Each window uses
  a seed derived from it and from the window position (Acycle
  deriveWindowSeed), so results are reproducible and independent of the
  other windows; `NULL` leaves the current random state unchanged.
  `Default=1`.

- window_mode:

  How the windows are defined (character): `"physical-depth"` (windows
  of fixed length in metres centred every `step_depth`) or
  `"legacy-count"` (windows of a fixed, even number of points,
  `window / dt + 1` rounded, at least 16, moved by `step` samples).
  `Default="physical-depth"`.

- dt:

  Nominal sampling interval, in metres (numeric, \> 0), used for the
  number of points per window in `"legacy-count"` mode and for the
  defaults of `step_depth` and `sr`. `NULL` uses the median spacing of
  the cleaned data. `Default=NULL`.

- batch_size:

  Number of Monte Carlo simulations processed together (integer, \>= 1);
  larger batches are faster but use more memory. `Default=20`.

- genplot:

  Plot the result with
  [`plot_ecoco_cv`](https://stratigraphy.github.io/WaverideR/reference/plot_ecoco_cv.md)
  (logical, TRUE or FALSE). `Default=FALSE`.

- verbose:

  Print messages and a progress bar over the windows (logical, TRUE or
  FALSE). `Default=FALSE`.

## Value

A list (class "ecoco_cv") with `method`, `sr`, `depth`, rate-by-window
matrices `rho`, `p_global`, `p_local`, `norbit`, `pcoco`, `score`, lists
`forward` (odd -\> even), `backward` (even -\> odd), `consensus` and
`strict` (consensus curve with the larger directional p), `windows`
(data frame: bounds, point counts, parity phase, derived seed, success,
window-level p-values, fold spacings, all-nine rate ranges, AR(1) rho),
`best` and `track`.

## References

Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking variable
sedimentation rates and astronomical forcing in Phanerozoic paleoclimate
proxy series with evolutionary correlation coefficients and hypothesis
testing. Earth and Planetary Science Letters 501, 165-179.
[doi:10.1016/j.epsl.2018.08.041](https://doi.org/10.1016/j.epsl.2018.08.041)

## Author

Michiel Arts (R translation); original MATLAB code: Mingsong Li

## Examples

``` r
# \donttest{
set.seed(2)
dt <- 0.05
depth <- seq(0, 30, by = dt)
tt <- depth * 100 / 5
x <- sin(2 * pi * tt / 405.7) + 0.8 * sin(2 * pi * tt / 98.9) +
  0.5 * sin(2 * pi * tt / 40.99) + 0.4 * sin(2 * pi * tt / 23.68) +
  0.3 * rnorm(length(tt))
res <- ecoco_interleaved(cbind(depth, x), window = 15,
                         sr = seq(3, 8, 0.5), step_depth = 2, pad = 400,
                         n_sim = 20)
res$track[, c("depth", "sr")]
#>   depth  sr
#> 1   7.5 5.5
#> 2   9.5 5.5
#> 3  11.5 5.5
#> 4  13.5 3.5
#> 5  15.5 5.5
#> 6  17.5 5.5
#> 7  19.5 3.5
#> 8  21.5 3.5
# }
```
