# Plot an age-model check: the de-chirped maps, the rate corrections the cycles ask for and the flagged intervals

Plot an age-model check: the de-chirped maps, the rate corrections the
cycles ask for and the flagged intervals

## Usage

``` r
plot_agemodel_check(
  x,
  transforms = names(x$dechirped),
  view = c("tuned", "depth"),
  rate_unit = c(1, "depth units / kyr"),
  main = NULL,
  ...
)
```

## Arguments

- x:

  Age-model check to plot: the output of
  [`check_agemodel`](https://stratigraphy.github.io/WaverideR/reference/check_agemodel.md)
  (an object of class `wr_agemodel_check`). Required, no default.

- transforms:

  De-chirped maps to show, one panel each (character vector, a subset of
  `names(x$dechirped)`, i.e. of the transforms used in
  [`check_agemodel`](https://stratigraphy.github.io/WaverideR/reference/check_agemodel.md)).
  `Default=names(x$dechirped)` (all computed maps).

- view:

  Horizontal axis of the plots (character): `"tuned"` shows the maps in
  the time of the track (kyr), where the astronomical cycles should be
  horizontal lines; `"depth"` shows them mapped back to depth.
  `Default=c("tuned", "depth")`, i.e. `"tuned"`.

- rate_unit:

  Factor and label of the rate axis, as a vector `c(factor, label)`,
  e.g. `c(100, "cm/kyr")`. In the current version it has no effect on
  the plot (the correction panel shows a dimensionless correction
  factor). `Default=c(1, "depth units / kyr")`.

- main:

  Title of the top (proxy) panel (character); `NULL` uses a default
  title. `Default=NULL`.

- ...:

  Further arguments passed to
  [`plot_dechirped`](https://stratigraphy.github.io/WaverideR/reference/plot_dechirped.md)
  for the map panels (e.g. `palette`, `level`, `scale`, `outline_lwd`,
  `ylab`); `style`, `show` and `main` are set by this function.

## Value

Invisibly `x`.
