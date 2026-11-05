# Remove tracking points which were tracked in a wavelet spectra

Interactively select points for deletion With the
[`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
function it is possible to track points in a wavelet spectra, however
errors can be made and as such it is possible to delete these points
with the `delpts_tracked_period_wt` function. This function allows one
to select points for deletion. \#'

## Usage

``` r
delpts_tracked_period_wt(
  tracking_pts = NULL,
  wavelet = NULL,
  n.levels = 100,
  periodlab = "Period (metres)",
  x_lab = "depth (metres)",
  palette_name = "viridis",
  color_brewer = "viridis",
  colour_brewer = NULL
)
```

## Arguments

- tracking_pts:

  Tracked points to clean (matrix or data frame, at least two columns):
  depth in column 1 and the tracked period in column 2 (same units as
  the wavelet), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md).
  The points are drawn on the scalogram and the ones clicked are
  removed. Required (`Default=NULL` gives an error).

- wavelet:

  Scalogram on which the points are shown: the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (class `analyze.wavelet`). Required (`Default=NULL` gives an error).

- n.levels:

  Number of colour levels of the power scale (integer, \> 1).
  `Default=100`.

- periodlab:

  Label of the y-axis (character). `Default="Period (metres)"`.

- x_lab:

  Label of the x-axis (character). `Default="depth (metres)"`.

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

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

The results of the deletion of the tracking points is a matrix with 3
columns. The first column is depth/time The second column is the period
of the tracked cycle The third column is the sedimentation rate based on
the duration (in time) of the tracked cycle

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

#mag_track <- track_period_wavelet(astro_cycle = 405,
#                                   wavelet=mag_wt,
#                                   n.levels = 100,
#                                   periodlab = "Period (metres)",
#                                   x_lab = "depth (metres)"
#                                    palette_name ="rainbow",
#                                    color_brewer ="grDevices)

#load the mag_track_solution data set to get an example data set from which
#data points can be deleted


mag_track_corr <- delpts_tracked_period_wt(tracking_pts = mag_track_solution,
                                    wavelet = mag_wt,
                                    n.levels = 100,
                                    periodlab = "Period (metres)",
                                    x_lab = "depth (metres)",
                                   palette_name ="rainbow",
                                   color_brewer ="grDevices")
}# }
```
