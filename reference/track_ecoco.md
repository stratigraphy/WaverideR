# Track the sedimentation rate through an eCOCO map

Interactively select points in the map of an eCOCO result (sedimentation
rate against depth) to trace how the sedimentation rate changes along
the record. This is the eCOCO counterpart of
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md):
the map of one of the eCOCO statistics (by default the combined
correlation and significance score `pcoco`) is plotted, the local maxima
(ridges) are marked, and the sedimentation rate is tracked by clicking
ridge points.

Selected points are shown as white dots. Clicking a selected point again
deselects it (red dot). Press the right mouse button or Esc (depending
on the graphics device) to finish. An interactive device that supports
[`identify`](https://rdrr.io/r/graphics/identify.html) is needed (on
Windows this works out of the box, on macOS XQuartz or the Quartz
device, on Linux an X11 device).

The tracking can also be done without a mouse: pass the points as
`points`. They are snapped to the nearest ridge point, which makes a
tracking reproducible and scriptable.

## Usage

``` r
track_ecoco(
  ecoco = NULL,
  map = c("pcoco", "rho", "score", "p_global", "norbit"),
  track_peaks = TRUE,
  n_peaks = 3,
  rho_min = NULL,
  p_max = NULL,
  norbit_min = NULL,
  points = NULL,
  interpolate = FALSE,
  n.levels = 100,
  plot_horizontal = TRUE,
  plot_dir = TRUE,
  x_lab = "depth (metres)",
  y_lab = "sedimentation rate (cm/kyr)",
  palette_name = "viridis",
  genplot = FALSE
)
```

## Arguments

- ecoco:

  eCOCO result to track (required): the list returned by
  [`ecoco_adaptive`](https://stratigraphy.github.io/WaverideR/reference/ecoco_adaptive.md),
  [`ecoco_blocked`](https://stratigraphy.github.io/WaverideR/reference/ecoco_blocked.md)
  or
  [`ecoco_interleaved`](https://stratigraphy.github.io/WaverideR/reference/ecoco_interleaved.md),
  with the vectors `sr` (sedimentation rates, cm/kyr) and `depth`
  (window centres) and rate-by-window matrices (one row per rate, one
  column per window). `Default=NULL`.

- map:

  The eCOCO statistic used for tracking and plotting (character):
  `"pcoco"` (combined correlation and significance score), `"rho"`
  (correlation coefficient), `"score"`, `"p_global"` (global p-value,
  shown and tracked as minus log10 p) or `"norbit"` (number of
  contributing orbital periods); it must be present in `ecoco`.
  `Default=c("pcoco", "rho", "score", "p_global", "norbit")`, i.e.
  `"pcoco"`.

- track_peaks:

  Only allow the local maxima of the map along the sedimentation rate
  axis to be selected (logical); `FALSE` allows every grid cell.
  `Default=TRUE`.

- n_peaks:

  Maximum number of peaks per window that can be selected, the highest
  ones (integer, \> 0; `Inf` keeps all); only used with
  `track_peaks = TRUE`. `Default=3`.

- rho_min:

  Minimum correlation coefficient (`rho` of `ecoco`) a point must have
  to be selectable (numeric, or `NULL` for no criterion).
  `Default=NULL`.

- p_max:

  Maximum global p-value (`p_global` of `ecoco`) a point may have to be
  selectable (numeric, between 0 and 1, or `NULL` for no criterion).
  `Default=NULL`.

- norbit_min:

  Minimum number of contributing orbital periods (`norbit` of `ecoco`) a
  point must have to be selectable (integer, or `NULL` for no
  criterion). `Default=NULL`.

- points:

  Points to select without mouse input (data frame or matrix with depth
  in column 1 and sedimentation rate in cm/kyr in column 2, or `NULL`
  for interactive selection). Each point is snapped to the nearest
  selectable point (distance scaled by the ranges of depth and rate).
  `Default=NULL`.

- interpolate:

  Linearly interpolate the tracked sedimentation rate to all windows of
  the eCOCO result, constant outside the tracked interval (logical); the
  result then only has the columns `depth` and `sedrate`.
  `Default=FALSE`.

- n.levels:

  Number of colour levels of the map (integer, \> 1); the breaks are
  evenly spaced. `Default=100`.

- plot_horizontal:

  Orientation of the map (logical): depth on the x-axis (`TRUE`) or on
  the y-axis (`FALSE`). `Default=TRUE`.

- plot_dir:

  Direction of the depth axis (logical): `TRUE` draws depth increasing
  along the axis, `FALSE` reverses it (for example for borehole data
  with depth increasing downwards). `Default=TRUE`.

- x_lab:

  Label of the depth axis (character). `Default="depth (metres)"`.

- y_lab:

  Label of the sedimentation rate axis (character).
  `Default="sedimentation rate (cm/kyr)"`.

- palette_name:

  Name of a [`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html)
  palette for the colour scale, a character string, for example
  `"viridis"`, `"Plasma"`, `"Inferno"`, `"YlOrRd"`, `"Spectral"`,
  `"Blue-Red"`, `"Zissou 1"`; run
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  to list all 115 names. Low values get the first colour of the palette.
  Names of the viridis family (`"viridis"`, `"magma"`, `"inferno"`,
  `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`) give the
  colours of the viridis package, as in the wavelet plots.
  `Default="viridis"`.

- genplot:

  Plot the map with the selected points and the tracked curve when
  `points` is used (logical); the interactive mode always plots.
  `Default=FALSE`.

## Value

A data frame with the tracked windows: `depth`, `sedrate` (cm/kyr) and,
when available in `ecoco`, `rho`, `p_global`, `norbit` and the value of
the map. The result can be used like the output of
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)
(for example to build an age model).

## See also

[`ecoco_track`](https://stratigraphy.github.io/WaverideR/reference/ecoco_track.md),
[`ecoco_track_ridge`](https://stratigraphy.github.io/WaverideR/reference/ecoco_track_ridge.md),
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)

## Author

Michiel Arts. Inspired by
[`track_period`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)
and the tracking routines of Acycle (Mingsong Li); see
[`ecoco_track`](https://stratigraphy.github.io/WaverideR/reference/ecoco_track.md)
for the automatic version.

## Examples

``` r
# synthetic eCOCO-like map: sedimentation rate increases with depth
sr <- seq(1, 10, by = 0.25)
depth <- seq(0, 20, by = 0.5)
true_sr <- 4 + 2 * depth / 20
rho <- 0.8 * exp(-outer(sr, true_sr, "-")^2 / 0.5)
ecoco <- list(sr = sr, depth = depth, rho = rho, pcoco = 2 * rho,
              p_global = matrix(0.01, length(sr), length(depth)),
              norbit = matrix(9, length(sr), length(depth)))
# scripted tracking: points are snapped to the nearest ridge point
trk <- track_ecoco(ecoco, points = data.frame(depth = c(2, 10, 18), sr = c(4.2, 5, 5.8)))
trk
#>   depth sedrate     rho p_global norbit   pcoco
#> 1     2    4.25 0.79601     0.01      9 1.59202
#> 2    10    5.00 0.80000     0.01      9 1.60000
#> 3    18    5.75 0.79601     0.01      9 1.59202
# \donttest{
# interactive tracking (needs a graphics device with mouse input)
if (interactive()) {
  trk <- track_ecoco(ecoco)
}
# }
```
