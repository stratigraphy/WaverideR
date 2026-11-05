# Plot an evolutive MESA

Plots the spectral image of an evolutive maximum entropy spectrum
([`mesa_evolutive`](https://stratigraphy.github.io/WaverideR/reference/mesa_evolutive.md)):
one panel for every window length, side by side, with the position
(depth or time) on one axis and the period or frequency on the other, as
in the spectral image columns of log-analysis software. The spectrum of
every window position is divided by its maximum, so that the ridges of
the spectral peaks show along the whole record. The peaks of every
position, the tracked ridge
([`track_mesa`](https://stratigraphy.github.io/WaverideR/reference/track_mesa.md))
and the average amplitude spectrum
([`mesa_average`](https://stratigraphy.github.io/WaverideR/reference/mesa_average.md))
can be added.

## Usage

``` r
plot_mesa_evolutive(
  x,
  vertical = TRUE,
  xaxis = c("period", "frequency"),
  period_lim = NULL,
  frequency_lim = NULL,
  log_scale = FALSE,
  peaks = FALSE,
  peak_min = 0.3,
  track = NULL,
  average = FALSE,
  colors = NULL,
  position_lim = NULL,
  main = NULL,
  log_axis = TRUE,
  scale = "auto",
  palette = "viridis",
  normalise = c("window", "period", "none"),
  ...
)
```

## Arguments

- x:

  Evolutive maximum entropy spectrum to plot: the list returned by
  [`mesa_evolutive`](https://stratigraphy.github.io/WaverideR/reference/mesa_evolutive.md)
  (one or more window lengths). Required, no default.

- vertical:

  Orientation of the plot (logical). `TRUE`: the position (depth/time)
  is the vertical axis, increasing downwards as for a depth log, and the
  window lengths are panels side by side; `FALSE`: the position is the
  horizontal axis and the window lengths are stacked panels.
  `Default=TRUE`.

- xaxis:

  Spectral axis (character): `"period"` or `"frequency"` (cycles per
  unit of position). `Default="period"` (the first element of
  `c("period", "frequency")`).

- period_lim:

  Range of the period axis, in the unit of the position axis (numeric
  vector of length 2, `c(min, max)`); used when `xaxis = "period"`.
  `NULL` runs from twice the sampling interval to the longest window
  length. `Default=NULL`.

- frequency_lim:

  Range of the frequency axis, in cycles per unit of position (numeric
  vector of length 2); used when `xaxis = "frequency"`. `NULL` uses the
  full frequency range. `Default=NULL`.

- log_scale:

  Kept for compatibility (logical): `TRUE` is the same as
  `scale = "log"`; `FALSE` uses `scale`. `Default=FALSE`.

- peaks:

  Draw the spectral peaks of every window position as white points
  (logical, TRUE or FALSE). `Default=FALSE`.

- peak_min:

  Smallest relative power of a peak that is drawn with `peaks = TRUE`
  (numeric, 0 to 1, power relative to the maximum of the spectrum at the
  same position). `Default=0.3`.

- track:

  Ridge track(s) drawn as red lines on the panel of the matching window
  length: a `mesa_track` object returned by
  [`track_mesa`](https://stratigraphy.github.io/WaverideR/reference/track_mesa.md)
  or a list of them; `NULL` draws none. `Default=NULL`.

- average:

  Draw the average amplitude spectrum of every window length
  ([`mesa_average`](https://stratigraphy.github.io/WaverideR/reference/mesa_average.md))
  beside or below its image, with the periods (or frequencies) of its
  three strongest peaks labelled (logical). `Default=FALSE`.

- colors:

  Colours to interpolate for the colour scale (character vector of
  colours, e.g. `c("white", "orange", "red")`); if given it overrides
  `palette`. `NULL` uses `palette`. `Default=NULL`.

- position_lim:

  Range of the position (depth/time) axis (numeric vector of length 2);
  `NULL` uses the full range of window positions. `Default=NULL`.

- main:

  Overall title, drawn in the outer top margin (character); `NULL` draws
  none. `Default=NULL`.

- log_axis:

  Draw the period or frequency axis on a log2 scale (logical): the
  spectrum of a cyclostratigraphic record spans several octaves, and on
  a linear axis the long periods take most of the plot. `Default=TRUE`.

- scale:

  Spacing of the colour breaks of the colour scale, as in all image
  plots of the package (character): `"quantile"` (histogram-equalised),
  `"log"`, `"sqrt"`, `"linear"` or `"auto"` (logarithmic for mostly
  empty maps or values spanning more than three orders of magnitude,
  quantiles otherwise; see
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)).
  A maximum-entropy spectrum spans many decades, so `"auto"` is usually
  logarithmic for it. `Default="auto"`.

- palette:

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

- normalise:

  How the spectrum is scaled before plotting (character), one of:
  `"window"` (the spectrum of each window position divided by its
  maximum, so the ridges show along the whole record), `"period"` (each
  period divided by its mean along the record, which removes the
  red-noise slope and shows where a cycle is stronger or weaker) or
  `"none"` (raw power). `Default="window"` (the first element of
  `c("window", "period", "none")`).

- ...:

  Not used; accepted (and ignored) so that extra arguments do not give
  an error.

## Value

Called for its plot; invisibly returns `x`.

## Examples

``` r
set.seed(2)
t <- seq(0, 800, by = 2)
y <- sin(2 * pi * cumsum(1 / ifelse(t < 400, 40, 80)) * 2) + rnorm(length(t), sd = 0.4)
ev <- mesa_evolutive(data.frame(t, y), window = c(150, 250), order = 20)
plot_mesa_evolutive(ev, average = TRUE)
```
