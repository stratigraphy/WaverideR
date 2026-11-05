# Extract a signal using standard deviation

Extract signal from a wavelet spectra in the depth domain using a the
standard deviation of the omega (number of cycles) as boundaries. The
uncertainty is based on the Gabor uncertainty principle applied to the
continuous wavelet transform using a Morlet wavelet. The calculated
uncertainty is the underlying analytical uncertainty which is the result
of applying the Gabor uncertainty principle to the continuous wavelet
transform using a Morlet wavelet.

## Usage

``` r
extract_signal_standard_deviation(
  wavelet = NULL,
  tracked_cycle_curve = NULL,
  multi = 1,
  extract_cycle = NULL,
  tracked_cycle_period = NULL,
  add_mean = TRUE,
  tune = FALSE,
  genplot_uncertainty_wt = FALSE,
  genplot_extracted = FALSE,
  keep_editable = FALSE,
  palette_name = "viridis",
  color_brewer = "viridis",
  colour_brewer = NULL
)
```

## Arguments

- wavelet:

  Wavelet object (list) created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  the signal is reconstructed from its coefficients and its number of
  cycles `omega_nr` sets the Gabor (Morlet) uncertainty. Required
  (`Default=NULL`).

- tracked_cycle_curve:

  Tracked period of the cycle along the record (matrix or data frame:
  depth/time in column 1, period in column 2, in the depth/time unit of
  the wavelet), e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  completed with
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md).
  Rows with `NA` are dropped and the curve is linearly interpolated (and
  extrapolated) to the depth/time axis of the wavelet. Required
  (`Default=NULL`).

- multi:

  Multiple of the standard deviation of the analytical frequency
  uncertainty of the Morlet wavelet that sets the band boundaries,
  `1 / (f0 -/+ multi * sd)` around the frequency `f0` of the (rescaled)
  tracked period (numeric, \> 0); larger values give a wider band.
  `Default=1`.

- extract_cycle:

  Duration of the cycle to extract (numeric, \> 0), in the same unit as
  `tracked_cycle_period`; use the same value to extract the tracked
  cycle itself. The tracked curve is multiplied by
  `extract_cycle / tracked_cycle_period`. Required (`Default=NULL`).

- tracked_cycle_period:

  Duration in time of the tracked cycle, e.g. 405 (kyr) for the long
  eccentricity cycle (numeric, \> 0), in the same unit as
  `extract_cycle`; also used for the tuning when `tune=TRUE`. Required
  (`Default=NULL`).

- add_mean:

  Add the mean of the original data series to the extracted cycle
  (logical, TRUE or FALSE); `FALSE` returns it centred on zero.
  `Default=TRUE`.

- tune:

  Convert the extracted cycle (and the plots) from the depth to the time
  domain with
  [`curve2tune`](https://stratigraphy.github.io/WaverideR/reference/curve2tune.md),
  using `tracked_cycle_curve` and `tracked_cycle_period` (logical, TRUE
  or FALSE). `Default=FALSE`.

- genplot_uncertainty_wt:

  Plot the wavelet scalogram (with
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md))
  with the tracked curve (black) and the upper (red) and lower (blue)
  boundaries set by `multi` times the standard deviation of the
  uncertainty superimposed (logical, TRUE or FALSE). `Default=FALSE`.

- genplot_extracted:

  Plot the data with the extracted cycle (blue) on top of it (logical,
  TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the plot settings active after the call so that extra features
  can be added to the plots (logical, TRUE or FALSE); `FALSE` restores
  the previous graphical parameters. `Default=FALSE`.

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

Signal extracted from the wavelet spectra. Output is a matrix with the
first column being depth/time and the second column is the astronomical
cycle extracted from the proxy record

If `genplot_uncertainty_wt=TRUE` then a wavelet spectra will be plotted
with the uncertainty superimposed on top of it. In the plot the red
curve and blue curves are the upper and lower bounds based on the
`multi` parameter.The black curve is the `Default=tracked_cycle_curve`
curve. If `genplot_extracted=TRUE` plot with the data set and the
extracted cycle on top of it will be plotted.

## References

Angi Roesch and Harald Schmidbauer (2018). WaveletComp: Computational
Wavelet Analysis. R package version 1.1.
<https://CRAN.R-project.org/package=WaveletComp>

Gouhier TC, Grinsted A, Simko V (2021). R package biwavelet: Conduct
Univariate and Bivariate Wavelet Analyses. (Version 0.20.21),
<https://github.com/tgouhier/biwavelet>

Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet
Analysis. Bulletin of the American Meteorological Society 79:61-78.
<https://paos.colorado.edu/research/wavelets/bams_79_01_0061.pdf>

Gabor, Dennis. "Theory of communication. Part 1: The analysis of
information." Journal of the Institution of Electrical Engineers-part
III: radio and communication engineering 93, no. 26 (1946): 429-441.
[doi:10.1049/ji-3-2.1946.0074](https://doi.org/10.1049/ji-3-2.1946.0074)

Russell, Brian, and Jiajun Han. "Jean Morlet and the continuous wavelet
transform. " CREWES Res. Rep 28 (2016): 115.
<https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf>

Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard. "Wave
propagation and sampling theory—Part I: Complex signal and scattering in
multilayered media. " Geophysics 47, no. 2 (1982): 203-221.

J. Morlet, G. Arens, E. Fourgeau, D. Giard; Wave propagation and
sampling theory; Part II, Sampling theory and complex waves. Geophysics
1982 47 (2): 222–236.

## Author

Code based on the reconstruct function of the 'WaveletComp' R package
which is based on the wavelet 'MATLAB' code written by Christopher
Torrence and Gibert P. Compo (1998). The assignment of the standard
deviation of the uncertainty of the wavelet is based on the work of
Gabor (1946) and Russell et al., (2016)

## Examples

``` r
# \donttest{
#Extract the 405 kyr eccentricity cycle from the magnetic susceptibility
#record of the Sullivan core of Pas et al., (2018) and use the Gabor
# uncertainty principle to define the mathematical uncertainty of the
# analysis and use a factor of that standard deviation to define
# boundaries

# perform the CWT
mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)

#Track the 405 kyr eccentricity cycle in a wavelet spectra

#mag_track <- track_period_wavelet(astro_cycle = 405,
#                                   wavelet=mag_wt,
#                                   n.levels = 100,
#                                   periodlab = "Period (metres)",
#                                   x_lab = "depth (metres)",
#                                   palette_name="rainbow",
#                                   color_brewer="grDevices")

#Instead of tracking, the tracked solution data set mag_track_solution is used
mag_track <- mag_track_solution

mag_track_complete <- completed_series(
  wavelet = mag_wt,
  tracked_curve = mag_track,
  period_up = 1.2,
  period_down = 0.8,
  extrapolate = TRUE,
  genplot = FALSE
)

# smooth the tracking of the 405 kyr eccentricity cycle
mag_track_complete <- loess_auto(time_series = mag_track_complete,
genplot = FALSE, print_span = FALSE)

# extract the 405 kyr eccentricity cycle from the wavelet spectrum and use
# the Gabor uncertainty principle to define the mathematical uncertainty of
# the analysis and use a multiple of the derived standard deviation to define boundaries

mag_405_ecc <- extract_signal_standard_deviation(
wavelet = mag_wt,
tracked_cycle_curve = mag_track_complete,
multi = 1,
extract_cycle = 405,
tracked_cycle_period = 405,
add_mean = TRUE,
tune = FALSE,
genplot_uncertainty_wt = FALSE,
genplot_extracted = FALSE,
keep_editable=FALSE,
palette_name="rainbow",
color_brewer="grDevices"
)
# }
```
