# Calculate sum of maximum spectral power for sedimentation rates for a wavelet spectra

The `sum_power_sedrate` function is used calculate the sum of maximum
spectral power for a list of astronomical cycles from a wavelet spectra.
The data is first normalised using the average spectral power curves for
a given percentile based on results of the
[`model_red_noise_wt`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md)
function

## Usage

``` r
sum_power_sedrate(
  red_noise = NULL,
  wavelet = NULL,
  percentile = NULL,
  sedrate_low = NULL,
  sedrate_high = NULL,
  spacing = NULL,
  cycles = c(NULL),
  x_lab = "depth",
  y_lab = "sedrate",
  run_multicore = FALSE,
  n_cores = NULL,
  genplot = FALSE,
  plot_res = 1,
  keep_editable = FALSE,
  palette_name = "viridis",
  color_brewer = "viridis",
  verbose = FALSE,
  colour_brewer = NULL
)
```

## Arguments

- red_noise:

  Red-noise power curves (required): the output of
  [`model_red_noise_wt`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md)
  for the same `wavelet`, a matrix with one row per period of `wavelet`
  and one column per red-noise simulation. `Default=NULL`.

- wavelet:

  Wavelet object (required) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
  of a record with depth in metres. `Default=NULL`.

- percentile:

  Percentile of the red-noise power (numeric, between 0 and 1) by which
  the wavelet power is weighted: at every period the threshold is the
  `percentile` quantile of a normal distribution with the mean and
  standard deviation of the `red_noise` runs. This accounts for the
  distortion of the power spectrum by the analytical technique and by
  random red noise. Required. `Default=NULL`.

- sedrate_low:

  Lowest sedimentation rate tested, in cm/kyr (numeric, \> 0). Required.
  `Default=NULL`.

- sedrate_high:

  Highest sedimentation rate tested, in cm/kyr (numeric, \>
  `sedrate_low`). Required. `Default=NULL`.

- spacing:

  Step between the tested sedimentation rates, in cm/kyr (numeric, \>
  0). Required. `Default=NULL`.

- cycles:

  Periods of the astronomical cycles, in kyr (numeric vector, e.g.
  `c(405, 125, 95, 40, 21)`), whose expected periods in depth
  (sedimentation rate x period) are searched for in the scalogram.
  Required. `Default=c(NULL)`.

- x_lab:

  Label of the depth axis of the plot (character). `Default="depth"`.

- y_lab:

  Label of the sedimentation rate axis of the plot (character).
  `Default="sedrate"`.

- run_multicore:

  Run the calculation for the positions along the record in parallel on
  several cores (logical). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer, or `NULL`).
  `NULL` uses an automatic, safe number: physical cores minus one, at
  most 4, and limited by the free memory; set
  `options(waverider.n_cores = ...)` to change it. A value \> 1 also
  switches on parallel computing. Results do not depend on the number of
  cores. `Default=NULL`.

- genplot:

  Plot the result (selected with `plot_res`) as a depth versus
  sedimentation rate map with its average curve (logical).
  `Default=FALSE`.

- plot_res:

  Which result to plot when `genplot = TRUE` (integer): `1` the
  normalised sum of maximum spectral power, `2` the number of
  astronomical components found. `Default=1`.

- keep_editable:

  Keep the graphical parameters of the plot after the function returns,
  so that extra elements can be added (logical); only used with
  `genplot = TRUE`. `Default=FALSE`.

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

- verbose:

  Show the progress of the calculation (logical, TRUE or FALSE).
  `Default=FALSE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

Returns a list which contains 4 elements element 1: sum of maximum
spectral power element 2: number of cycles used in the sum of maximum
spectral power element 3: y-axis values of the matrices which is
sedimentation rate element 4: x-axis values of the matrices which is
depth

If `Default="TRUE"` a plot is created with 3 subplots. Subplot 1 is plot
in which the the sum of maximum spectral power for a given sedimentation
rate or nr of cycles is plotted for each depth given depth. Subplot 2 is
a plot in which the average sum of maximum spectral power is plotted fro
each sedimentation Subplot 3 is a color scale for subplot 1.

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

Acycle: Time-series analysis software for paleoclimate research and
education, Mingsong Li, Linda Hinnov, Lee Kump, Computers &
Geosciences,Volume 127,2019,Pages 12-22,ISSN 0098-3004,
\<doi:10.1016/j.cageo.2019.02.011\>

Tracking variable sedimentation rates and astronomical forcing in
Phanerozoic paleoclimate proxy series with evolutionary correlation
coefficients and hypothesis testing, Mingsong Li, Lee R. Kump, Linda A.
Hinnov, Michael E. Mann, Earth and Planetary Science Letters,Volume 501,
T2018,Pages 165-179,ISSN 0012-821X,\<doi:10.1016/j.epsl.2018.08.041\>

## Author

Based on the [asm](https://rdrr.io/pkg/astrochron/man/asm.html) and
[eAsm](https://rdrr.io/pkg/astrochron/man/eAsm.html) functions of the
'astrochron' R package and the 'eCOCO' and 'COCO' functions of the
'Acycle' software

## Examples

``` r
# \donttest{
#estimate sedimentation rate for the the magnetic susceptibility record
# of the Sullivan core of Pas et al., (2018).

mag_wt <- analyse_wavelet(data = mag[1:200, ],
dj = 1/10,
lowerPeriod = 0.1,
upperPeriod = 32,
verbose = FALSE,
omega_nr = 10)

#increase n_simulations to better define the red noise spectral power curve
mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
n_simulations = 3,
run_multicore=FALSE,
verbose=FALSE)

sedrates <- sum_power_sedrate(red_noise=mag_wt_red_noise,
wavelet=mag_wt,
percentile=0.75,
sedrate_low = 0.5,
sedrate_high = 4,
spacing = 1,
cycles = c(2376,1600,1180,696,406,110),
x_lab="depth",
y_lab="sedrate",
run_multicore=FALSE,
genplot = FALSE,
plot_res=1,
keep_editable=FALSE,
palette_name = "rainbow",
color_brewer="grDevices",
verbose=FALSE)
# }

```
