# Plots a cross superlet scalogram

Plot cross superlet scalogram using the outcome of the
[`analyse_Xsuperlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md)
function.

## Usage

``` r
plot_Xsuperlet(
  Xsuperlet = NULL,
  lowerPeriod = NULL,
  upperPeriod = NULL,
  n.levels = 100,
  palette_name = "viridis",
  color_brewer = "viridis",
  useRaster = TRUE,
  periodlab = "Period (metres)",
  x_lab = "depth (metres)",
  keep_editable = FALSE,
  dev_new = wr_dev_new_default(),
  plot_dir = TRUE,
  plot_arrows = TRUE,
  pwr_quant = 0.85,
  n_arrows = 75,
  T_unit = NULL,
  P_unit = NULL,
  arrow_length = 0.075,
  arrow_lwd = 1.9,
  arrow_angle = 25,
  add_lines = NULL,
  add_points = NULL,
  add_abline_h = NULL,
  add_abline_v = NULL,
  add_data = TRUE,
  add_avg = FALSE,
  plot_horizontal = TRUE,
  colour_brewer = NULL
)
```

## Arguments

- Xsuperlet:

  Cross superlet result to plot: the list returned by
  [`analyse_Xsuperlet`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md)
  (computed from the two records `data_1` and `data_2`). Required
  (`Default=NULL` gives an error).

- lowerPeriod:

  Shortest period shown on the period axis, in the depth/time unit of
  the input records (numeric, \> 0). `NULL` uses the shortest period of
  `Xsuperlet`. `Default=NULL`.

- upperPeriod:

  Longest period shown on the period axis, in the depth/time unit of the
  input records (numeric, \> `lowerPeriod`). `NULL` uses the longest
  period of `Xsuperlet`. `Default=NULL`.

- n.levels:

  Number of colour levels of the power scale (integer, \>= 2); the level
  breaks are the quantiles of the cross power, so each colour covers the
  same share of the scalogram. `Default=100`.

- palette_name:

  Name of the colour palette of the power (colour) scale, a character
  string. Which names are valid depends on `color_brewer`:

  - `color_brewer = "grDevices"`: `"rainbow"`, `"heat.colors"`,
    `"terrain.colors"`, `"topo.colors"`, `"cm.colors"`, or any of the
    115 [`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html)
    palettes (list them with
    [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)),
    for example `"viridis"`, `"Plasma"`, `"Inferno"`, `"YlOrRd"`,
    `"Spectral"`, `"Blue-Red"`, `"Zissou 1"`.

  - `color_brewer = "viridis"`: `"viridis"`, `"magma"`, `"inferno"`,
    `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`.

  - `color_brewer = "RColorBrewer"`: `"BrBG"`, `"PiYG"`, `"PRGn"`,
    `"PuOr"`, `"RdBu"`, `"RdGy"`, `"RdYlBu"`, `"RdYlGn"`, `"Spectral"`,
    `"Accent"`, `"Dark2"`, `"Paired"`, `"Pastel1"`, `"Pastel2"`,
    `"Set1"`, `"Set2"`, `"Set3"`, `"Blues"`, `"BuGn"`, `"BuPu"`,
    `"GnBu"`, `"Greens"`, `"Greys"`, `"Oranges"`, `"OrRd"`, `"PuBu"`,
    `"PuBuGn"`, `"PuRd"`, `"Purples"`, `"RdPu"`, `"Reds"`, `"YlGn"`,
    `"YlGnBu"`, `"YlOrBr"`, `"YlOrRd"` (also listed by
    [`wr_brewer_palettes()`](https://stratigraphy.github.io/WaverideR/reference/wr_coloured_noise.md)).

  - `color_brewer = "colorRamps"`: `"blue2green"`, `"blue2green2red"`,
    `"blue2red"`, `"blue2yellow"`, `"cyan2yellow"`, `"green2red"`,
    `"magenta2green"`, `"matlab.like"`, `"matlab.like2"`, `"ygobb"`.

  - `color_brewer = "scico"` (all except `"batlowK"` need the scico
    package): `"acton"`, `"bam"`, `"bamako"`, `"bamO"`, `"batlow"`,
    `"batlowK"`, `"batlowW"`, `"berlin"`, `"bilbao"`, `"broc"`,
    `"brocO"`, `"buda"`, `"bukavu"`, `"cork"`, `"corkO"`, `"davos"`,
    `"devon"`, `"fes"`, `"glasgow"`, `"grayC"`, `"hawaii"`, `"imola"`,
    `"lajolla"`, `"lapaz"`, `"lipari"`, `"lisbon"`, `"managua"`,
    `"navia"`, `"nuuk"`, `"oleron"`, `"oslo"`, `"roma"`, `"romaO"`,
    `"tofino"`, `"tokyo"`, `"turku"`, `"vanimo"`, `"vik"`, `"vikO"`.

  Low power is drawn with the first colour of the palette and high power
  with the last. `Default="viridis"`.

- color_brewer:

  Palette family from which `palette_name` is taken, a character string:
  `"grDevices"`, `"viridis"`, `"RColorBrewer"`, `"colorRamps"` or
  `"scico"` (the last one, except for `"batlowK"`, needs the scico
  package). See `palette_name` for the palette names of each family.
  `Default="viridis"`.

- useRaster:

  Draw the scalogram as a raster image (`TRUE`) or as vector polygons
  (`FALSE`) (logical). WARNING: plotting as a vector image is
  computationally intensive and gives large files. `Default=TRUE`.

- periodlab:

  Label of the period axis (character). `Default="Period (metres)"`.

- x_lab:

  Label of the depth/time axis (character). `Default="depth (metres)"`.

- keep_editable:

  Keep the plot layout and graphical parameters (`par`) after plotting,
  so that extra elements can be added to the scalogram panel with
  `lines`, `points` etc. (logical). With `FALSE` the previous `par`
  settings are restored on exit. `Default=FALSE`.

- dev_new:

  Open a new plotting window for the plot (logical, TRUE or FALSE). The
  wavelet plots contain a raster image that needs more room than the
  small RStudio plot pane, which gives "figure margins too large" errors
  or clipped text, so a separate window is opened by default in the
  RStudio GUI. In every other front end (plain R, ESS, Positron, knitr,
  Rscript) the plot is drawn on the current device and no window is
  opened. Set TRUE or FALSE to force either behaviour, or set
  `options(waverider.dev_new = TRUE/FALSE)`.
  `Default=wr_dev_new_default()`: TRUE in the RStudio GUI, FALSE
  elsewhere.

- plot_dir:

  Direction of the depth/time axis (logical). `TRUE` plots increasing
  depth/time values from left to right (or bottom to top when
  `plot_horizontal=FALSE`); `FALSE` reverses the axis. `Default=TRUE`.

- plot_arrows:

  Draw the phase arrows of the cross superlet transform on the scalogram
  (logical). An arrow is drawn at each local maximum of cross power
  along the period axis whose power exceeds the `pwr_quant` quantile;
  its direction gives the phase difference between the two records
  (radians, measured anticlockwise from the depth/time axis).
  `Default=TRUE`.

- pwr_quant:

  Quantile of the cross power (numeric, between 0 and 1) above which
  phase arrows are drawn; higher values keep only the arrows on the
  strongest power ridges. `Default=0.85`.

- n_arrows:

  Approximate number of arrow positions along the depth/time axis
  (integer, \>= 1); every `round(n/n_arrows)`-th depth/time step of the
  `n` steps gets arrows. `Default=75`.

- T_unit:

  Length of the arrows along the depth/time axis, in depth/time units
  (numeric, \> 0). `NULL` (or `P_unit=NULL`) gives arrows of a fixed
  visual length of 0.3 inch whatever the record length. `Default=NULL`.

- P_unit:

  Length of the arrows along the period axis, in log2(period) units
  (numeric, \> 0). `NULL` (or `T_unit=NULL`) uses the automatic arrow
  length (see `T_unit`). `Default=NULL`.

- arrow_length:

  Length of the arrow head edges in inches (numeric, \> 0; `length`
  argument of [`arrows`](https://rdrr.io/r/graphics/arrows.html)).
  `Default=0.075`.

- arrow_lwd:

  Line width of the arrows (numeric, \> 0). `Default=1.9`.

- arrow_angle:

  Angle between the arrow shaft and the edges of the arrow head, in
  degrees (numeric, 0-90). `Default=25`.

- add_lines:

  Lines to draw on the scalogram: a matrix or data frame with depth/time
  in column 1 and, in each further column, period values (same unit as
  the period axis) of one line; `NULL` draws none. `Default=NULL`.

- add_points:

  Points to draw on the scalogram: a matrix or data frame with
  depth/time in column 1 and, in each further column, period values of
  one series of points; `NULL` draws none. `Default=NULL`.

- add_abline_h:

  Numeric vector of positions of horizontal lines, e.g. `c(2, 3, 5, 6)`:
  periods when `plot_horizontal=TRUE`, depth/time values when
  `plot_horizontal=FALSE`; `NULL` draws none. `Default=NULL`.

- add_abline_v:

  Numeric vector of positions of vertical lines, e.g. `c(2, 3, 5, 6)`:
  depth/time values when `plot_horizontal=TRUE` (also drawn on the data
  panel), periods when `plot_horizontal=FALSE`; `NULL` draws none.
  `Default=NULL`.

- add_data:

  Plot the two input records (record 1 in black, record 2 in red) in a
  panel beside the scalogram (logical). `Default=TRUE`.

- add_avg:

  Plot the time-averaged cross power spectrum in a panel beside the
  scalogram (logical). `Default=FALSE`.

- plot_horizontal:

  Orientation of the plot (logical): `TRUE` puts depth/time on the x
  axis and period on the y axis; `FALSE` puts depth/time on the y axis
  and period on the x axis. `Default=TRUE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

The output is a plot of a cross superlet scalogram.

## References

Roesch, A., & Schmidbauer, H. (2018). WaveletComp: Computational Wavelet
Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
Time-frequency super-resolution with superlets. Nature Communications,
12(1), 337.
[doi:10.1038/s41467-020-20539-9](https://doi.org/10.1038/s41467-020-20539-9)

## Author

plotting code based on the "wt.image" function of the 'WaveletComp' R
package Whereas the "analyse_Xsuperlet" that generates the input for the
plotting function is based on the matlab code in Moca et al. (2021) and
the the "analyze.coherency" function of the 'WaveletComp' R package

## Examples

``` r
# \donttest{
#Example 1. A cross superlet plot of two etp solutions with noise overprint


etp_1 <- astrochron::etp(
 tmin = 0,
 tmax = 500,
 dt = 1,
 eWt = 1.5,
 oWt = 0.75,
 pWt = 1,
 esinw = TRUE,
 standardize = TRUE,
 genplot = FALSE,
 verbose = FALSE
)

etp_2 <- astrochron::etp(
 tmin = 0,
 tmax = 500,
 dt = 1,
 eWt = 1,
 oWt = 0.5,
 pWt = 1.5,
 esinw = TRUE,
 standardize = TRUE,
 genplot = FALSE,
verbose = FALSE
)

etp_1[, 2] <- etp_1[, 2] + wr_coloured_noise(
 nrow(etp_1),
 sd = sd(etp_1[, 2]) / 1.5,
 mean = mean(etp_1[, 2]),
 phi = 0.9
)
etp_2[, 2] <- etp_2[, 2] + wr_coloured_noise(
 nrow(etp_2),
 sd = sd(etp_2[, 2]) / 1.5,
 mean = mean(etp_2[, 2]),
 phi = 0.9
)

X_super_etp <- analyse_Xsuperlet(
 data_1 = etp_1,
 data_2  = etp_2,
 upperPeriod = 1024,
 lowerPeriod = 2,
 Nf = 32,
 c1 = 3,
 o = c(1, 10),
 mult = TRUE,
 verbose = FALSE
)
plot_Xsuperlet(Xsuperlet = X_super_etp, lowerPeriod = 2, upperPeriod = 1024,
              n.levels = 100, palette_name = "rainbow", color_brewer = "grDevices",
              useRaster = TRUE, periodlab = "Period (metres)", x_lab = "depth (metres)",
              keep_editable = FALSE, dev_new = TRUE, plot_dir = TRUE,
              plot_arrows = TRUE,
              pwr_quant = 0.85,
              n_arrows = 50,
              T_unit = NULL,
              P_unit  = NULL,
              arrow_length = 0.1,
              arrow_lwd = 2,
              arrow_angle = 25,add_lines = NULL, add_points = NULL, add_abline_h = NULL,
              add_abline_v = NULL, add_data = TRUE, add_avg = FALSE, plot_horizontal = TRUE)

# }
```
