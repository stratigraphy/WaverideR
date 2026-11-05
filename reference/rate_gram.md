# Rate-gram: the astronomical template as a spectrum over sedimentation rate, per window

A purely spectral estimate of the sedimentation rate along a record,
without any tracking: for every window of
[`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md)
(`criterion = "template"`) the template score (the summed log
line-to-continuum ratio at the periods of all astronomical cycles, the
cyclostratigraphic counterpart of the harmonic product spectrum used for
pitch detection in speech and music) is taken at its best gradient and
shown against the rate; its maximum per window is the spectral rate
estimate, optionally median-filtered over a few windows. On CIP Case 3
the maximum per window correlates 0.56 with the true rate, the median
over 11 windows 0.73 (the tracked rate of
[`track_chirp_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_chirp_hmm.md),
which adds continuity, 0.97).

## Usage

``` r
rate_gram(
  ce,
  k = 11,
  plot = TRUE,
  rate_unit = c(1, "depth units / kyr"),
  reference = NULL,
  main = NULL
)
```

## Arguments

- ce:

  Chirp-matched EHA result: the list returned by
  [`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md)
  run with `criterion = "template"` (it must contain
  `template_score_all` and `sedrate_grid`). Required, no default.

- k:

  Width, in windows, of the running median applied to the best rate per
  window (integer, \>= 1; even values are rounded up to the next odd
  number, 1 means no smoothing). `Default=11`.

- plot:

  Draw the rate-gram (map of the template score against sedimentation
  rate and depth, with the best rate per window as white points and its
  running median as a white line) (logical, TRUE or FALSE); with `FALSE`
  only the list is returned. `Default=TRUE`.

- rate_unit:

  Conversion factor and label of the rate axis: a vector
  `c(factor, "label")` whose first element multiplies the rates (in
  depth units per kyr) for the axis labels, e.g. `c(100, "cm/kyr")` for
  a depth axis in metres. Only changes the labelling.
  `Default=c(1, "depth units / kyr")`.

- reference:

  Reference sedimentation-rate curve drawn as a red dashed line for
  comparison: a two-column matrix or data frame with depth in column 1
  and rate (depth units per kyr) in column 2; `NULL` draws none.
  `Default=NULL`.

- main:

  Title of the plot (character); `NULL` gives a default title.
  `Default=NULL`.

## Value

A list with `depth`, `sedrate` (rate grid), `map` (rates x windows),
`best` (the maximum per window) and `smooth` (its running median).
