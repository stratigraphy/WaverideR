# Build a track (depth, sedimentation rate, time) from a tracked period curve, a rate curve or an age model

The tracking functions return `wr_track` objects; this builds one from
the curves older workflows produce, so that they can be given to
[`analyse_dechirped`](https://stratigraphy.github.io/WaverideR/reference/analyse_dechirped.md)
and
[`check_agemodel`](https://stratigraphy.github.io/WaverideR/reference/check_agemodel.md):
a tracked period curve
([`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md):
depth and the period of one astronomical cycle in depth units), a
sedimentation-rate curve (e.g. from an EHA/sedrate workflow) or an
age-depth model.

## Usage

``` r
as_wr_track(
  x,
  type = c("period", "sedrate", "age"),
  cycle = 405,
  rate_factor = 1,
  step = NULL
)
```

## Arguments

- x:

  Curve to convert (data frame or matrix, two columns; further columns
  are ignored): depth in column 1 and, in column 2, the tracked period
  (in depth units) for `type="period"`, the sedimentation rate (depth
  units per kyr after multiplying by `rate_factor`) for
  `type="sedrate"`, or the age (kyr after multiplying by `rate_factor`)
  for `type="age"`. Rows with non-finite values are dropped, duplicated
  depths are averaged, and at least two points are needed. Required, no
  default.

- type:

  What column 2 of `x` holds (character): `"period"` (the period of one
  astronomical cycle in depth units, e.g. from
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md);
  the rate is period / `cycle`), `"sedrate"` (a sedimentation rate; must
  be positive) or `"age"` (an age-depth model; must be strictly
  monotonic, and time is counted from the top of the curve).
  `Default=c("period", "sedrate", "age")`, i.e. `"period"`.

- cycle:

  Duration of the tracked astronomical cycle in kyr (numeric, \> 0),
  used to turn the tracked period into a sedimentation rate when
  `type="period"`; ignored otherwise. `Default=405`.

- rate_factor:

  Multiplier applied to column 2 of `x` to bring it to depth units per
  kyr (`type="sedrate"`) or to kyr (`type="age"`) (numeric); e.g. `0.01`
  for a rate in cm/kyr with depth in metres, `1000` for an age in Ma.
  Ignored when `type="period"`. `Default=1`.

- step:

  Depth step of the output track, in depth units (numeric, \> 0); the
  rate is interpolated (log-linearly) onto a regular grid with this step
  and the time is recalculated. `NULL` keeps the depth positions of `x`.
  `Default=NULL`.

## Value

A `wr_track` data frame with `depth`, `sedrate` and `time` (kyr from the
top of the curve, increasing with depth).

## Examples

``` r
curve <- data.frame(depth = 0:100, period = 12 + 3 * sin(0:100 / 10))
tr <- as_wr_track(curve, "period", cycle = 405)
head(tr)
#>   depth    sedrate      time
#> 1     0 0.02962963   0.00000
#> 2     1 0.03036914  33.33908
#> 3     2 0.03110125  65.87969
#> 4     3 0.03181867  97.67026
#> 5     4 0.03251421 128.76219
#> 6     5 0.03318093 159.20898
```
