# Plots an log2 spaced EHA spectra

Plot wavelet scalogram using the outcome of the
[`eha_log2`](https://stratigraphy.github.io/WaverideR/reference/eha_log2.md)
function.

## Usage

``` r
plot_eha_log2(
  eha_log2 = NULL,
  plot_opt = "Amplitude",
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
  add_lines = NULL,
  add_points = NULL,
  add_abline_h = NULL,
  add_abline_v = NULL,
  add_MTM_peaks = FALSE,
  add_data = TRUE,
  add_avg = FALSE,
  pval_abline = c(0.1, 0.05),
  pval_cutoff = c(0.1),
  add_MTM = FALSE,
  mtm_siglvl = 0.95,
  demean_mtm = TRUE,
  detrend_mtm = TRUE,
  padfac_mtm = 5,
  tbw_mtm = 3,
  plot_horizontal = TRUE,
  colour_brewer = NULL,
  scale = "auto"
)
```

## Arguments

- eha_log2:

  Evolutive harmonic analysis result to plot: the list returned by
  [`eha_log2`](https://stratigraphy.github.io/WaverideR/reference/eha_log2.md).
  Required (`Default=NULL` gives an error).

- plot_opt:

  Quantity shown in the spectrogram and in the average panel
  (character), one of `"Amplitude"`, `"Power"`, `"Probability"`
  (harmonic F-test probability) or `"F_test"` (F-test value); any other
  value plots the amplitude. `Default="Amplitude"`.

- lowerPeriod:

  Shortest period shown on the period axis, in the depth/time unit of
  the input record (numeric, \> 0). `NULL` uses the shortest period of
  `eha_log2`. `Default=NULL`.

- upperPeriod:

  Longest period shown on the period axis, in the depth/time unit of the
  input record (numeric, \> `lowerPeriod`). `NULL` uses the longest
  period of `eha_log2`. `Default=NULL`.

- n.levels:

  Number of colour levels of the colour scale (integer, \>= 2).
  `Default=100`.

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

  Draw the spectrogram as a raster image (`TRUE`) or as vector polygons
  (`FALSE`) (logical). WARNING: plotting as a vector image is
  computationally intensive and gives large files. `Default=TRUE`.

- periodlab:

  Label of the period axis (character). `Default="Period (metres)"`.

- x_lab:

  Label of the depth/time axis (character). `Default="depth (metres)"`.

- keep_editable:

  Keep the plot layout and graphical parameters (`par`) after plotting,
  so that extra elements can be added to the spectrogram panel with
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

- add_lines:

  Lines to draw on the spectrogram: a matrix or data frame with
  depth/time in column 1 and, in each further column, period values
  (same unit as the period axis) of one line; `NULL` draws none.
  `Default=NULL`.

- add_points:

  Points to draw on the spectrogram: a matrix or data frame with
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

- add_MTM_peaks:

  Run a multitaper (MTM) analysis of the whole record
  ([`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)) and draw the
  periods of its significant peaks (at `mtm_siglvl`) as dotted lines on
  the spectrogram (and in red on the MTM panels); the MTM result with a
  `Peak` column (1 = significant peak) is then returned invisibly
  (logical). `Default=FALSE`.

- add_data:

  Plot the input record in a panel beside the spectrogram (logical).
  `Default=TRUE`.

- add_avg:

  Plot the average of the plotted quantity (`plot_opt`) over the record
  in a panel beside the spectrogram (logical). `Default=FALSE`.

- pval_abline:

  Not used; kept for compatibility with older versions (numeric vector
  of p-values). `Default=c(0.1, 0.05)`.

- pval_cutoff:

  Not used; kept for compatibility with older versions (numeric
  p-value). `Default=c(0.1)`.

- add_MTM:

  Add three panels beside the spectrogram with the MTM power spectrum of
  the whole record (with its AR(1) confidence levels), the AR(1)
  confidence level and the harmonic F-test confidence level (in %),
  computed with [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)
  (logical). `Default=FALSE`.

- mtm_siglvl:

  Significance level used to select the MTM peaks (numeric, between 0
  and 1; `siglevel` of
  [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)). `Default=0.95`.

- demean_mtm:

  Remove the mean from the record before the MTM analysis (logical;
  `demean` of [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)).
  `Default=TRUE`.

- detrend_mtm:

  Remove a linear trend from the record before the MTM analysis
  (logical; `detrend` of
  [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)). `Default=TRUE`.

- padfac_mtm:

  Zero-padding factor of the MTM analysis (numeric, \>= 1; `padfac` of
  [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)). `Default=5`.

- tbw_mtm:

  Time-bandwidth product of the MTM analysis (numeric, e.g. 2, 3 or 4;
  `tbw` of [`mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)); higher
  values give smoother spectra with lower frequency resolution.
  `Default=3`.

- plot_horizontal:

  Orientation of the plot (logical): `TRUE` puts depth/time on the x
  axis and period on the y axis; `FALSE` puts depth/time on the y axis
  and period on the x axis. `Default=TRUE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

- scale:

  Spacing of the colour breaks of the colour scale, as in all image
  plots of the package (character): `"quantile"` (histogram-equalised),
  `"log"`, `"sqrt"`, `"linear"` or `"auto"` (logarithmic for mostly
  empty maps or values spanning more than three orders of magnitude,
  quantiles otherwise); see
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md).
  `Default="auto"`.

## Value

The output is a plot of a EHA spectra. if add_MTM_peaks = TRUE then the
output of the MTM analysis will given as matrix

## References

Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
Wavelet Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

S.R. Meyers, 2012, Seeing Red in Cyclic Stratigraphy: Spectral Noise
Estimation for Astrochronology: Paleoceanography, 27, PA3228,
\<doi:10.1029/2012PA002307\>

S.R. Meyers, 2019 Cyclostratigraphy and the problem of astrochronologic
testing, Earth-Science Reviews, Volume 190,
\<doi.org/10.1016/j.earscirev.2018.11.015.\>

## Author

Code based on the wt.image" functions of the 'WaveletComp' R package The
EHA and MTM analysis parts are from the astrochron R package of Meyers
et al., (2012)

## Examples

``` r
# \donttest{
#Example 1. A plot of a wavelet spectra using the Total Solar Irradiance
# data set of Steinhilber et al., (2012)

TSI_eha_log2 <- eha_log2(data = TSI,
win =  8192,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 8192,
lowerPeriod = 50,
pad = NULL,
padding = "noise")

plot_eha_log2(
eha_log2 = TSI_eha_log2,
plot_opt = "Amplitude",
lowerPeriod = 50,
upperPeriod = 8192,
n.levels = 100,
palette_name = "rainbow",
color_brewer = "grDevices",
useRaster = TRUE,
periodlab = "Period (metres)",
x_lab = "depth (metres)",
keep_editable = FALSE,
dev_new = TRUE,
plot_dir = TRUE,
add_lines = NULL,
add_points = NULL,
add_abline_h = NULL,
add_abline_v = NULL,
add_MTM_peaks = FALSE,
add_data = TRUE,
add_avg = FALSE,
add_MTM = FALSE,
mtm_siglvl = 0.95,
demean_mtm = TRUE,
detrend_mtm = TRUE,
padfac_mtm = 5,
tbw_mtm = 3,
plot_horizontal = TRUE)


#Example 2. A plot of a wavelet spectra using the magnetic susceptibility
#data set of Pas et al., (2018)
mag_eha_log2 <- eha_log2(data = mag,
win =  50,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 50,
lowerPeriod = 1,
pad = NULL,
padding = "noise")

plot_eha_log2(
eha_log2 = mag_eha_log2,
plot_opt = "Amplitude",
lowerPeriod = 1,
upperPeriod = 50,
n.levels = 100,
palette_name = "rainbow",
color_brewer = "grDevices",
useRaster = TRUE,
periodlab = "Period (metres)",
x_lab = "depth (metres)",
keep_editable = FALSE,
dev_new = TRUE,
plot_dir = TRUE,
add_lines = NULL,
add_points = NULL,
add_abline_h = NULL,
add_abline_v = NULL,
add_MTM_peaks = FALSE,
add_data = TRUE,
add_avg = FALSE,
add_MTM = FALSE,
mtm_siglvl = 0.95,
demean_mtm = TRUE,
detrend_mtm = TRUE,
padfac_mtm = 5,
tbw_mtm = 3,
plot_horizontal = TRUE)



#Example 3. A plot of a wavelet spectra using the greyscale
# data set of Zeeden et al., (2013)
grey_eha_log2 <- eha_log2(data = grey,
win =  20,
tbw = 4,
demean = TRUE,
detrend = TRUE,
upperPeriod = 20,
lowerPeriod = 1,
pad = NULL,
padding = "noise")


plot_eha_log2(
eha_log2 = grey_eha_log2,
plot_opt = "Amplitude",
lowerPeriod = 1,
upperPeriod = 20,
n.levels = 100,
palette_name = "rainbow",
color_brewer = "grDevices",
useRaster = TRUE,
periodlab = "Period (metres)",
x_lab = "depth (metres)",
keep_editable = FALSE,
dev_new = TRUE,
plot_dir = TRUE,
add_lines = NULL,
add_points = NULL,
add_abline_h = NULL,
add_abline_v = NULL,
add_MTM_peaks = FALSE,
add_data = TRUE,
add_avg = FALSE,
add_MTM = FALSE,
mtm_siglvl = 0.95,
demean_mtm = TRUE,
detrend_mtm = TRUE,
padfac_mtm = 5,
tbw_mtm = 3,
plot_horizontal = TRUE)

# }
```
