# Extract a signal in between tracked boundaries in a wavelet scalogram

Interactively select points in a wavelet scalogram to trace the upper
and lower period of an cycle. The `dynamic_extraction` function plots a
wavelet scalogram in which points peaks can selected allowing one to
track the lower and upper period of a cycle. First track the upper or
lower period of the to be extracted cycle and then track the other
boundary. Tracking points can be selected in the Interactive interface
and will be shown as white dots connected by a black line. When one
wants to deselect a point the white dots can be re-clicked/re-selected
and will turn red which indicates that the previously selected point is
deselected. Deselecting points can be quite tricky. After tracking the
lower and upper boundaries of the cycle the `dynamic_extraction`
function will extract the signal in between the boundaries. the output
can then used as input for the
[`minimal_tuning`](https://stratigraphy.github.io/WaverideR/reference/minimal_tuning.md)
function to create an age model.

## Usage

``` r
dynamic_extraction(
  wavelet = NULL,
  n.levels = 100,
  add_peaks = FALSE,
  periodlab = "Period (metres)",
  x_lab = "depth (metres)",
  palette_name = "viridis",
  color_brewer = "viridis",
  plot_horizontal = TRUE,
  smooth = FALSE,
  add_mean = TRUE,
  colour_brewer = NULL
)
```

## Arguments

- wavelet:

  Wavelet object (list) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  its scalogram is plotted for the tracking and its reconstructed
  components are summed between the tracked boundaries. Required
  (`Default=NULL`).

- n.levels:

  Number of colour levels of the power scale (integer, \> 1); the levels
  are quantiles of the power. `Default=100`.

- add_peaks:

  Add the spectral peaks (local power maxima along the period axis of
  each column of the scalogram) as points to the tracking plot, to help
  the tracking (logical, TRUE or FALSE). `Default=FALSE`.

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

  Orientation of the tracking plot (logical): `TRUE` draws depth/time on
  the x-axis and period on the y-axis, `FALSE` draws period on the
  x-axis and depth/time on the y-axis. `Default=TRUE`.

- smooth:

  Smooth each of the two tracked boundaries with
  [`loess_auto`](https://stratigraphy.github.io/WaverideR/reference/loess_auto.md)
  after interpolating them to the depth/time grid of the wavelet
  (logical, TRUE or FALSE). `Default=FALSE`.

- add_mean:

  Add the mean of the original data series to the extracted signal
  (logical, TRUE or FALSE); `FALSE` returns the signal centred on zero.
  `Default=TRUE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

Results of the tracking of a cycle in the wavelet spectra is a matrix
with 3 columns. The first column is depth/time The second column is the
extracted tracked cycle The third column is upper tracked period The
fourth column is lower tracked period

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
if (FALSE) { # \dontrun{
#Track the 405kyr upper and lower periods of the eccentricity cycle in the
#magnetic susceptibility record of the Sullivan core of Pas et al., (2018)

mag_wt <- analyse_wavelet(
 data = mag,
 dj = 1 / 100,
 lowerPeriod = 0.1,
 upperPeriod = 254,
 verbose = FALSE,
 omega_nr = 10
)

mag_ext <- dynamic_extraction(
 wavelet = mag_wt,
 n.levels = 100,
 add_peaks = FALSE,
 periodlab = "Period (metres)",
 x_lab = "depth (metres)",
 palette_name = "rainbow",
 color_brewer = "grDevices",
 plot_horizontal = TRUE,
 smooth = TRUE,
 add_mean = TRUE
)
} # }
```
