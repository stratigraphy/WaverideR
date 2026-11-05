# Plots a cross wavelet scalogram

Plot cross wavelet scalogram using the outcome of the
[`analyse_wavelet_coherence`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_coherence.md)
function.

## Usage

``` r
plot_wavelet_coherence(
  wavelet_coh = NULL,
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

- wavelet_coh:

  Wavelet coherence object to plot (required): the list returned by
  [`analyse_wavelet_coherence`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_coherence.md).
  Its `Coherence`, `Phase`, `Period`, `x1`/`y1` and `x2`/`y2` elements
  are used. `Default=NULL`.

- lowerPeriod:

  Shortest period shown on the period axis, in the depth/time unit of
  the records (numeric, \> 0). `NULL` uses the shortest period of
  `wavelet_coh`. `Default=NULL`.

- upperPeriod:

  Longest period shown on the period axis, in the depth/time unit of the
  records (numeric, \> `lowerPeriod`). `NULL` uses the longest period of
  `wavelet_coh`. `Default=NULL`.

- n.levels:

  Number of colour levels of the coherence scale (integer, \> 1); the
  breaks are quantiles of the coherence. `Default=100`.

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

  Draw the scalogram as a raster image (`TRUE`) or as a vector image of
  individual cells (`FALSE`) (logical). Vector images are slow to draw
  and give large files. `Default=TRUE`.

- periodlab:

  Label of the period axis (character). `Default="Period (metres)"`.

- x_lab:

  Label of the depth/time axis (character). `Default="depth (metres)"`.

- keep_editable:

  Keep the graphical parameters (`par`, layout) of the plot after the
  function returns, so that extra elements can be added to the plot
  afterwards (logical). With `FALSE` the previous graphical parameters
  are restored on exit. `Default=FALSE`.

- dev_new:

  Open a new plotting window for the plot. The wavelet plots contain a
  raster image that needs more room than the small RStudio plot pane,
  which gives "figure margins too large" errors or clipped text, so a
  separate window is opened by default in the RStudio GUI. In every
  other front end (plain R, ESS, Positron, knitr, Rscript) the plot is
  drawn on the current device and no window is opened. Set TRUE or FALSE
  to force either behaviour, or set
  `options(waverider.dev_new = TRUE/FALSE)`. Logical.
  `Default=wr_dev_new_default()`: TRUE in the RStudio GUI, FALSE
  elsewhere.

- plot_dir:

  Direction of the depth/time axis (logical): `TRUE` draws it with
  values increasing from left to right (or top to bottom in a vertical
  plot), `FALSE` reverses the axis. `Default=TRUE`.

- plot_arrows:

  Draw phase arrows on the scalogram (logical). Arrows are drawn at the
  local coherence maxima along the period axis whose coherence is above
  the `pwr_quant` quantile; their direction gives the phase difference
  between the two records (pointing right: in phase, left: anti-phase).
  `Default=TRUE`.

- pwr_quant:

  Quantile of the coherence (numeric, between 0 and 1) above which phase
  arrows are drawn; higher values give fewer arrows. `Default=0.85`.

- n_arrows:

  Approximate number of arrow positions along the depth/time axis
  (integer, \> 0); the columns of the scalogram are thinned to about
  this number. `Default=75`.

- T_unit:

  Horizontal size of the arrows in depth/time units (numeric). `NULL`
  (or `P_unit = NULL`) uses an automatic, fixed visual length of 0.3
  inch whatever the record length. `Default=NULL`.

- P_unit:

  Vertical size of the arrows in log2 period units (numeric). `NULL` (or
  `T_unit = NULL`) uses the automatic size (see `T_unit`).
  `Default=NULL`.

- arrow_length:

  Length of the edges of the arrow head, in inches (numeric; `length` of
  [`arrows`](https://rdrr.io/r/graphics/arrows.html)). `Default=0.075`.

- arrow_lwd:

  Line width of the arrows (numeric). `Default=1.9`.

- arrow_angle:

  Angle between the shaft and the edges of the arrow head, in degrees
  (numeric; `angle` of
  [`arrows`](https://rdrr.io/r/graphics/arrows.html)). `Default=25`.

- add_lines:

  Lines to draw on the scalogram, e.g. tracked periods (matrix or data
  frame, or `NULL` for none): depth/time in column 1 and one or more
  columns of period values after it; each period column is drawn as a
  separate line. `Default=NULL`.

- add_points:

  Points to draw on the scalogram (matrix or data frame, or `NULL` for
  none): depth/time in column 1 and one or more columns of period values
  after it. `Default=NULL`.

- add_abline_h:

  Positions of extra horizontal lines (numeric vector, e.g.
  `c(2, 3, 5)`, or `NULL` for none). In a horizontal plot these are
  periods; in a vertical plot (`plot_horizontal = FALSE`) they are
  depth/time values. `Default=NULL`.

- add_abline_v:

  Positions of extra vertical lines (numeric vector or `NULL` for none).
  In a horizontal plot these are depth/time values; in a vertical plot
  they are periods. `Default=NULL`.

- add_data:

  Draw the two data series in a panel next to the scalogram, the second
  one in red (logical). `Default=TRUE`.

- add_avg:

  Draw the average coherence per period (`Coh.avg`) in a panel to the
  side of the scalogram (logical). `Default=FALSE`.

- plot_horizontal:

  Orientation of the plot (logical): `TRUE` puts depth/time on the
  x-axis and period on the y-axis; `FALSE` draws a vertical plot with
  depth/time on the y-axis. `Default=TRUE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

The output is a plot of a cross wavelet scalogram.

## References

Roesch, A., & Schmidbauer, H. (2018). WaveletComp: Computational Wavelet
Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A., & Mureșan, R. C. (2021).
Time-frequency super-resolution with superlets. Nature Communications,
12(1), 337.
[doi:10.1038/s41467-020-20539-9](https://doi.org/10.1038/s41467-020-20539-9)

## Author

plotting code based on the "wt.image" and "analyze.coherency" functions
of the 'WaveletComp' R package

## Examples

``` r
# \donttest{
#Example 1. A cross superlet plot of two etp solutions with noise overprint


etp_1 <- astrochron::etp(
 tmin = 0,
 tmax = 1500,
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
 tmax = 1500,
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

coh_etp <- analyse_wavelet_coherence(
 data_1 = etp_1,
 data_2  = etp_2,
 upperPeriod = 1024,
 lowerPeriod = 2
)
plot_wavelet_coherence(wavelet_coh = coh_etp, lowerPeriod = 2, upperPeriod = 1024,
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
