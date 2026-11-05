# Track the period of a cycle in a wavelet spectra

Interactively select points in a wavelet spectra to trace a period in a
wavelet spectra. The `track_period_wavelet` function plots a wavelet
spectra in which spectral peaks can selected allowing one to track a
ridge hence one can track the a cycle with a changing period. Tracking
points can be selected in the Interactive interface and will be shown as
white dots when one wants to deselect a point the white dots can be
re-clicked/re-selected and will turn red which indicates that the
previously selected point is deselected. Deselecting points can be quite
tricky due to the close spacing of points and such the
[`delpts_tracked_period_wt`](https://stratigraphy.github.io/WaverideR/reference/delpts_tracked_period_wt.md)
can be used to delete points were previously selected using the
`track_period_wavelet` function.

## Usage

``` r
track_period_wavelet(
  wavelet = NULL,
  astro_cycle = 405,
  n.levels = 100,
  track_peaks = TRUE,
  periodlab = "Period (metres)",
  x_lab = "depth (metres)",
  palette_name = "viridis",
  color_brewer = "viridis",
  plot_horizontal = TRUE,
  plot_dir = TRUE,
  lowerPeriod = NULL,
  upperPeriod = NULL,
  add_lines = NULL,
  add_points = NULL,
  add_abline_h = NULL,
  add_abline_v = NULL,
  colour_brewer = NULL
)
```

## Arguments

- wavelet:

  Wavelet object (class `analyze.wavelet`) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  in which the cycle is tracked; the elements `x`, `y`, `Power`,
  `Period`, `axis.2`, `coi.1`, `coi.2` and `Power.avg` are used. Must be
  given (the `NULL` default fails). `Default=NULL`.

- astro_cycle:

  Duration of the astronomical cycle that is tracked (numeric, \> 0, in
  kyr, e.g. 405 for the long eccentricity cycle). It is used to convert
  the tracked period into the sedimentation rate of the output (period /
  `astro_cycle` \* 100, i.e. cm/kyr for a record in metres).
  `Default=405`.

- n.levels:

  Number of colour levels of the power scale (integer, \> 1); the breaks
  are quantiles of the power. `Default=100`.

- track_peaks:

  Restrict the selection to spectral peaks (logical, TRUE or FALSE).
  With `TRUE` the local power maxima along the period axis are drawn as
  black circles and a click selects the nearest one; with `FALSE` any
  point of the wavelet spectrum can be clicked (the click is snapped to
  the nearest depth of the wavelet and its period is used as clicked).
  In both modes clicking a selected point again removes it; finish with
  right-click/Esc (or "Finish" in RStudio). `Default=TRUE`.

- periodlab:

  Label of the period axis (character). `Default="Period (metres)"`.

- x_lab:

  Label of the depth/time axis (character). `Default="depth (metres)"`.

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

- plot_horizontal:

  Layout of the plot (logical, TRUE or FALSE): `TRUE` depth/time on the
  x-axis and period on the y-axis; `FALSE` period on the x-axis and
  depth/time on the y-axis. The output is the same in both layouts.
  `Default=TRUE`.

- plot_dir:

  Direction of the depth/time axis (logical, TRUE or FALSE): `TRUE`
  plots it from low to high values (e.g. borehole data that get older
  with increasing depth); `FALSE` reverses the axis (e.g. outcrop logs
  in which 0 m is the base of the section). `Default=TRUE`.

- lowerPeriod:

  Shortest period shown on the period axis (numeric, in the depth/time
  unit of the record); `NULL` uses the shortest period of `wavelet`.
  `Default=NULL`.

- upperPeriod:

  Longest period shown on the period axis (numeric, in the depth/time
  unit of the record); `NULL` uses the longest period of `wavelet`.
  `Default=NULL`.

- add_lines:

  Optional lines to draw on the wavelet plot (matrix or data frame:
  depth/time in column 1, one or more columns of periods after it; one
  line per period column). `Default=NULL`.

- add_points:

  Optional points to draw on the wavelet plot (matrix or data frame:
  depth/time in column 1, one or more columns of periods after it).
  `Default=NULL`.

- add_abline_h:

  Optional horizontal reference lines (numeric vector, e.g.
  `c(2, 3, 5)`), drawn at these depths/times; only used in the vertical
  layout (`plot_horizontal = FALSE`), the horizontal layout ignores it.
  `Default=NULL`.

- add_abline_v:

  Optional vertical reference lines (numeric vector), drawn at these
  periods; only used in the vertical layout (`plot_horizontal = FALSE`),
  the horizontal layout ignores it. `Default=NULL`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

Results of the tracking of a cycle in the wavelet spectra is a matrix
with 3 columns. The first column is depth/time The second column is the
period of the tracked cycle The third column is the sedimentation rate
based on the duration (in time) of the tracked cycle

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

The function is based/inspired on the
[traceFreq](https://rdrr.io/pkg/astrochron/man/traceFreq.html) function
of the 'astrochron' R package

## Examples

``` r
# \donttest{
if (interactive()) { # needs mouse input (identify())
#Track the 405kyr eccentricity cycle in the magnetic susceptibility record
# of the Sullivan core of Pas et al., (2018)

mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)

mag_track <- track_period_wavelet(wavelet = mag_wt,
astro_cycle = 405,
n.levels = 100,
track_peaks = TRUE,
periodlab =  "Period (metres)",
x_lab = "depth (metres)",
palette_name = "rainbow",
color_brewer = "grDevices",
plot_horizontal = TRUE,
plot_dir = TRUE,
lowerPeriod = NULL,
upperPeriod = NULL,
add_lines = NULL,
add_points = NULL,
add_abline_h = NULL,
add_abline_v = NULL)
}# }
```
