# Fit linear models to spectral peaks extracted from the wavelet spectra to astronomical cycles multiplied by sedimentation rate x

The `flmw` function is used calculate the linear correlation for a list
of astronomical cycles transformed using a range of sedimentation rates
and then compared to spectral peaks of a wavelet spectra

## Usage

``` r
flmw(
  wavelet = NULL,
  sedrate_low = NULL,
  sedrate_high = NULL,
  spacing = NULL,
  cycles = c(NULL),
  x_lab = "depth",
  y_lab = "sedrate",
  run_random = FALSE,
  rand_simulations = 1000,
  run_multicore = FALSE,
  n_cores = NULL,
  genplot = FALSE,
  palette_name = "viridis",
  color_brewer = "viridis",
  plot_res = 2,
  keep_editable = FALSE,
  verbose = FALSE,
  method = c("classic", "weighted", "template"),
  tolerance = 0.1,
  whiten = TRUE,
  colour_brewer = NULL
)
```

## Arguments

- wavelet:

  Wavelet object of a depth record (depth in metres), the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (superlet output of class `"analyze.superlet"` is also accepted).
  Required (`Default=NULL`).

- sedrate_low:

  Lowest sedimentation rate to test, in cm/kyr (numeric, \> 0). Required
  (`Default=NULL`).

- sedrate_high:

  Highest sedimentation rate to test, in cm/kyr (numeric, \>
  `sedrate_low`). Required (`Default=NULL`).

- spacing:

  Step between the tested sedimentation rates, in cm/kyr (numeric, \>
  0); the trial rates are
  `seq(sedrate_low, sedrate_high, by = spacing)`. Required
  (`Default=NULL`).

- cycles:

  Periods of the astronomical cycles in kyr (numeric vector, e.g.
  `c(405, 125, 95, 37, 20)`). For each trial sedimentation rate they are
  converted to periods in metres (`cycles * sedrate / 100`) and compared
  with the wavelet spectrum at every depth. Required
  (`Default=c(NULL)`).

- x_lab:

  Label of the x-axis (depth) of the plot (character).
  `Default="depth"`.

- y_lab:

  Label of the y-axis (sedimentation rate) of the plot (character).
  `Default="sedrate"`.

- run_random:

  Only for `method = "classic"`: also run Monte Carlo simulations in
  which the spectral peak periods are randomly perturbed, to give the
  percentile of each fit statistic against this null hypothesis
  (logical, TRUE or FALSE); without it the percentile elements (5-8) of
  the output are `NA`. `Default=FALSE`.

- rand_simulations:

  Number of random simulations used when `run_random = TRUE` (integer,
  \> 0; more simulations give more stable percentiles but take longer).
  `Default=1000`.

- run_multicore:

  Fit the depth levels in parallel on several cores (logical, TRUE or
  FALSE); not used by `method = "template"`. Results do not depend on
  it. `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer or `NULL`).
  `NULL` gives an automatic, safe number (physical cores minus one, at
  most 4, limited by the free memory); change this default with
  `options(waverider.n_cores = ...)`. A value \> 1 also switches on
  parallel computing. Results do not depend on the number of cores.
  `Default=NULL`.

- genplot:

  Plot the result matrix (sedimentation rate against depth) with its
  depth-averaged profile and a colour key (logical, TRUE or FALSE).
  `Default=FALSE`.

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

- plot_res:

  Which result matrix is plotted when `genplot = TRUE` (integer 1-8;
  only used with `method = "classic"`, the other methods always plot the
  score): 1 slope coefficient, 2 r squared, 3 number of matched
  components, 4 difference to the origin, 5 slope coefficient
  percentile, 6 r squared percentile, 7 number of components percentile,
  8 difference to the origin percentile (5-8 need `run_random = TRUE`).
  `Default=2`.

- keep_editable:

  Keep the graphical parameters (`par`) of the plot after the function
  returns, so that features can be added to it (logical, TRUE or FALSE);
  `FALSE` restores the previous settings. `Default=FALSE`.

- verbose:

  Show a progress bar over the depth levels (logical, TRUE or FALSE).
  `Default=FALSE`.

- method:

  Method used to compare the cycles with the spectrum (character), one
  of:

  - `"classic"`: the original peak-based linear-fit method (output
    identical to earlier versions);

  - `"weighted"`: matches the predicted cycles one-to-one to spectral
    peaks, weighting by peak power and using a relative (log period)
    misfit;

  - `"template"`: correlates the spectrum with Gaussian peaks at the
    predicted periods (no peak picking, fastest and usually most
    robust).

  `Default=c("classic", "weighted", "template")`, i.e. `"classic"`.

- tolerance:

  Relative period tolerance used by the `"weighted"` and `"template"`
  methods: the standard deviation of the Gaussian misfit in natural-log
  period (numeric, \> 0; 0.1 is about 10 percent). `Default=0.1`.

- whiten:

  For the `"weighted"` and `"template"` methods only: divide the wavelet
  power of every period by its average over depth before the comparison,
  which removes the red-noise background (logical, TRUE or FALSE).
  `Default=TRUE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

Returns a list. For `method = "classic"` it has the 10 elements below
plus `score_mat` (a heuristic combined score between 0 and 1),
`best_sedrate` (sedimentation rate with the highest score per depth) and
`method`. For the other methods the list has `score_mat`,
`best_sedrate`, `depth`, `y_axis` and `method`. Use
[`flmw_path`](https://stratigraphy.github.io/WaverideR/reference/flmw_path.md)
to turn `score_mat` into a smooth sedimentation rate curve with
uncertainty.

Elements of the classic result: element 1: slope coefficient element 2:
r squared element 3: nr of components element 4: difference to the
origin element 5: slope coefficient percentile element 6: r squared
percentile element 7: nr of components percentile, element 8: difference
to the origin percentile element 9: y-axis values of the matrices which
is sedimentation rate element 10: x-axis values of the matrices which is
depth

## Choosing a method

`"classic"` keeps the earlier results reproducible. In a synthetic test
with strong noise and a sedimentation rate that changes with depth,
`"weighted"` and `"template"` followed by
[`flmw_path`](https://stratigraphy.github.io/WaverideR/reference/flmw_path.md)
recovered the true curve far better than the classic score (see the NEWS
file). `"template"` needs no peak picking and is the fastest.

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

Based on the [eAsm](https://rdrr.io/pkg/astrochron/man/eAsm.html)
function of the 'astrochron' R package and the 'eCOCO' and 'COCO'
function of the 'Acycle' software

## Examples

``` r
# \donttest{
#estimate sedimentation rate for the magnetic susceptibility record
# of the Sullivan core of Pas et al., (2018).

mag_wt <- analyse_wavelet(data = mag[1:200, ],
dj = 1/10,
lowerPeriod = 0.1,
upperPeriod = 32,
verbose = FALSE,
omega_nr = 10)

sedrates <- flmw(wavelet = mag_wt,
    sedrate_low = 0.5,
    sedrate_high = 4,
    spacing = 1,
    cycles = c(2376,1600,1180,696,406,110),
    x_lab = "depth",
    y_lab = "sedrate",
    run_random = FALSE,
    rand_simulations = 50, # increase to get better constrainted resutls
    run_multicore = FALSE,
    genplot = FALSE,
    palette_name = "rainbow",
    color_brewer = "grDevices",
    plot_res = 2,
    keep_editable=FALSE,
    verbose=FALSE)

# template method (no peak picking) followed by a smooth sedimentation rate curve
sed_t <- flmw(wavelet = mag_wt, sedrate_low = 0.5, sedrate_high = 4, spacing = 0.5,
    cycles = c(406, 110), method = "template")
sed_path <- flmw_path(sed_t)
# }
```
