# Windowed fft based spectral analysis

The `win_fft` function conducts a windowed spectral analysis based using
the fft

## Usage

``` r
win_fft(
  data = NULL,
  padfac = 5,
  window_size = NULL,
  run_multicore = FALSE,
  n_cores = NULL,
  genplot = FALSE,
  x_lab = c("depth (m)"),
  y_lab = c("frequency cycle/metre"),
  plot_res = 1,
  perc_vis = 0,
  freq_max = NULL,
  freq_min = NULL,
  palette_name = "viridis",
  color_brewer = "viridis",
  keep_editable = FALSE,
  verbose = FALSE,
  dev_new = FALSE,
  colour_brewer = NULL
)
```

## Arguments

- data:

  Record to analyse (matrix or data frame with two columns: depth/time
  in column 1, proxy value in column 2), evenly sampled (the sampling
  interval is taken from the first two rows; rows with missing values
  are removed). Must be given (the `NULL` default fails).
  `Default=NULL`.

- padfac:

  Zero-padding factor of each window (numeric, \>= 1): every window is
  padded with zeros to `padfac` times its length before the FFT, which
  gives a finer (interpolated) frequency grid and a smoother spectrum; 1
  or less means no padding. `Default=5`.

- window_size:

  Length of the running window (numeric, \> 0, in the depth/time unit of
  `data`, several sampling intervals long). A spectrum is computed for a
  window centred on every sample; at both ends the record is extended by
  half a window of AR(1) noise with the mean, standard deviation and
  lag-1 autocorrelation of the first and last half window. Must be
  given. `Default=NULL`.

- run_multicore:

  Compute the windows in parallel on several cores (logical, TRUE or
  FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer; `NULL`: an
  automatic, safe number: physical cores minus one, at most 4, and
  limited by the free memory). Set `options(waverider.n_cores = ...)` to
  change the default. A value \> 1 also switches on parallel computing.
  Results do not depend on the number of cores. `Default=NULL`.

- genplot:

  Plot the matrix chosen with `plot_res` against depth and frequency,
  with the mean spectrum and a colour key (logical, TRUE or FALSE);
  needs `freq_min` and `freq_max`. `Default=FALSE`.

- x_lab:

  Label of the x-axis (depth/time) of the plot (character).
  `Default="depth (m)"`.

- y_lab:

  Label of the y-axis (frequency) of the plot (character).
  `Default="frequency cycle/metre"`.

- plot_res:

  Which matrix of the result is plotted (integer, 1 to 9): 1 amplitude,
  2 power, 3 phase, 4 frequency uncertainty, 5 AR(1) confidence level
  (percent), 6 AR(1) fit, 7 AR(1) 90% power level, 8 AR(1) 95% power
  level, 9 AR(1) 99% power level. `Default=1`.

- perc_vis:

  Lower cut-off of the colour scale as a quantile of the plotted values
  (numeric, between 0 and 1): values below this quantile get the lowest
  colour, e.g. 0.5 shows only the upper half of the values in colour.
  `Default=0`.

- freq_max:

  Highest frequency shown in the plot (numeric, in cycles per depth/time
  unit, at most the Nyquist frequency); required when `genplot = TRUE`.
  `Default=NULL`.

- freq_min:

  Lowest frequency shown in the plot (numeric, in cycles per depth/time
  unit, \> 0); required when `genplot = TRUE`. `Default=NULL`.

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

- keep_editable:

  Keep the graphical parameters and layout after plotting, so that extra
  features can be added (logical, TRUE or FALSE); `FALSE` restores the
  previous settings. `Default=FALSE`.

- verbose:

  Show a progress bar of the window computations (logical, TRUE or
  FALSE). `Default=FALSE`.

- dev_new:

  Open a new plot window (via `wr_open_device()`) before plotting
  (logical, TRUE or FALSE); a separate window gives a "nice" looking
  plot, but in an R markdown document the plot may then not appear.
  `Default=FALSE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

Returns a list which contains 10 elements element 1: Amplitude matrix
element 2: Power matrix element 3: Phase matrix element 4: AR1_CL matrix
element 5: AR1_Fit matrix element 6: AR1_90_power matrix element 7:
AR1_95_power matrix element 8: AR1_99_power matrix element 9: depth
element 10: y_axis If genplot is `Default=TRUE` then a plot of one of
the elements 1:8 is plotted

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Based on the
[periodogram](https://rdrr.io/pkg/astrochron/man/periodogram.html)
function of the 'astrochron' R package.

## Examples

``` r
# \donttest{
#Conduct a windowed ftt on the magnetic susceptibility record
#of the Sullivan core of Pas et al., (2018).

mag_win_fft <- win_fft(data= mag,
                   padfac = 5,
                   window_size = 12.5,
                   run_multicore = FALSE,
                   genplot = FALSE,
                   x_lab = c("depth (m)"),
                   y_lab = c("frequency cycle/metre"),
                   plot_res = 1,
                   perc_vis = 0.5,
                   freq_max = 5,
                   freq_min = 0.001,
                   palette_name ="rainbow",
                   color_brewer= "grDevices",
                   keep_editable=FALSE,
                   verbose=FALSE,
                   dev_new=FALSE)
# }
```
