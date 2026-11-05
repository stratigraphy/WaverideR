# Calculate the uncertainty associated with the wavelet analysis based on the Gabor uncertainty principle

The `wavelet_uncertainty` function is used to calculate uncertainties
associated with the wavelet analysis based on the Gabor uncertainty
principle applied to the continuous wavelet transform using a Morlet
wavelet. The calculated uncertainty is the underlying analytical
uncertainty which is the result of applying the Gabor uncertainty
principle to the continuous wavelet transform using a Morlet wavelet.

## Usage

``` r
wavelet_uncertainty(
  tracked_cycle = NULL,
  period_of_tracked_cycle = NULL,
  wavelet = NULL,
  multi = 1,
  verbose = FALSE,
  genplot_time = FALSE,
  genplot_uncertainty = FALSE,
  genplot_uncertainty_wt = FALSE,
  keep_editable = FALSE,
  palette_name = "viridis",
  color_brewer = "viridis",
  colour_brewer = NULL
)
```

## Arguments

- tracked_cycle:

  Curve of the tracked cycle (data frame with depth or time in column 1
  and the period of the tracked cycle, in the same unit, in column 2;
  further columns are ignored), e.g. from
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md)
  completed with
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md)
  and smoothed with
  [`loess_auto`](https://stratigraphy.github.io/WaverideR/reference/loess_auto.md).
  Must be given (the `NULL` default fails). `Default=NULL`.

- period_of_tracked_cycle:

  Duration of the tracked cycle in time (numeric, \> 0, in kyr, e.g. 405
  for the long eccentricity cycle), passed to
  [`curve2time`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md)
  to convert the period curves into time. Must be given. `Default=NULL`.

- wavelet:

  Wavelet object created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  (class `analyze.wavelet`); its `omega_nr` sets the width of the Gabor
  uncertainty, and it is plotted with `genplot_uncertainty_wt`. Must be
  given. `Default=NULL`.

- multi:

  Multiple of the standard deviation of the Gabor uncertainty (in
  frequency) used for the upper and lower curves (numeric, \> 0; e.g. 2
  for about 95 percent). `Default=1`.

- verbose:

  Print the total duration, the uncertainty of the duration and the
  relative uncertainty per cycle (logical, TRUE or FALSE).
  `Default=FALSE`.

- genplot_time:

  Plot the depth-time curves (black: tracked period, red and blue:
  period plus and minus `multi` standard deviations of the Gabor
  uncertainty of the Morlet wavelet) (logical, TRUE or FALSE).
  `Default=FALSE`.

- genplot_uncertainty:

  Plot the tracked period against depth with the upper (red) and lower
  (blue) period curves of `multi` standard deviations (logical, TRUE or
  FALSE). `Default=FALSE`.

- genplot_uncertainty_wt:

  Plot the wavelet spectrum
  ([`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md))
  with the tracked period (black), the period plus (red) and minus
  (blue) the uncertainty of `multi` standard deviations, and the area in
  between shaded in grey (logical, TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters after plotting, so that extra features
  can be added to the last plot (logical, TRUE or FALSE); `FALSE`
  restores the previous settings. `Default=FALSE`.

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

Results pertaining to the uncertainty calculated based on the Gabor
uncertainty principle.  
If the genplot_time is TRUE then a depth time plot will be plotted with
3 lines, the mean age,age plus x times the standard deviation and age
minus x times the standard deviation .  
If the genplot_uncertainty is TRUE then a curve will be plotted with the
mean period, the tracked period plus x times the standard deviation and
the tracked period minus x times the standard deviation.  
If the genplot_uncertainty_wt is TRUE a wavelet spectra will be plotted
with the tracked period, the tracked period plus x times the standard
deviation,the tracked period minus x times the standard deviation and
the area in between will be shaded in grey.  

Returns a matrix with 8 columns.  
The first column is called "depth" eg. depth  
The second column is "period" of the originally tracked period.  
The third column is "frequency" of the originally tracked period.  
The fourth column "uncertainty in frequency FWHM" is the uncertainty in
frequency based on the Gabor uncertainty principle defined as (FWHM)
full width at half maximum.  
The fifth column "uncertainty in frequency x_times SD" is the
uncertainty in frequency based on the Gabor uncertainty principle
defined as times x standard deviations.  
The sixth column "time mean" is the mean time based on the tracked
period.  
The seventh column "time plus x_times sd" is the time based on the
tracked period plus x times the standard deviation.  
The eight column "time min x_times sd" is the time based on the tracked
period min x times the standard deviation.  

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

Gabor, Dennis. Theory of communication. Part 1: The analysis of
information. Journal of the Institution of Electrical Engineers-part
III: radio and communication engineering 93, no. 26 (1946): 429-441.
[doi:10.1049/ji-3-2.1946.0074](https://doi.org/10.1049/ji-3-2.1946.0074)

Russell, Brian, and Jiajun Han.Jean Morlet and the continuous wavelet
transform. CREWES Res. Rep 28 (2016): 115.
<https://www.crewes.org/Documents/ResearchReports/2016/CRR201668.pdf>

Morlet, Jean, Georges Arens, Eliane Fourgeau, and Dominique Glard. Wave
propagation and sampling theory—Part I: Complex signal and scattering in
multilayered media. Geophysics 47, no. 2 (1982): 203-221.

J. Morlet, G. Arens, E. Fourgeau, D. Giard; Wave propagation and
sampling theory; Part II, Sampling theory and complex waves. Geophysics
1982 47 (2): 222–236.

## Author

Code based on the "analyze.wavelet" function of the 'WaveletComp' R
package and "wt" function of the 'biwavelet' R package which are based
on the wavelet 'MATLAB' code written by Christopher Torrence and Gibert
P. Compo (1998). The assignment of the standard deviation of the
uncertainty of the wavelet is based on the work of Gabor (1946) and
Russell et al., (2016)

## Examples

``` r
# \donttest{
#calculate the Gabor uncertainty derived mathematical uncertainty of the
#magnetic susceptibility record of the Sullivan core of Pas et al., (2018)

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
#                                palette_name="rainbow",
#                                color_brewer= "grDevices")

#Instead of tracking, the tracked solution data set mag_track_solution is used
mag_track <- mag_track_solution

mag_track_complete <- completed_series(
  wavelet = mag_wt,
  tracked_curve = mag_track,
  period_up = 1.2,
  period_down = 0.8,
  extrapolate = FALSE,
  genplot = FALSE,
  keep_editable=FALSE
)

mag_track_complete <- loess_auto(time_series = mag_track_complete,
 genplot = FALSE, print_span = FALSE,keep_editable=FALSE)

uncertainty <- wavelet_uncertainty(
  tracked_cycle = mag_track_complete,
  period_of_tracked_cycle = 405,
  wavelet = mag_wt,
  multi=1,
  verbose = FALSE,
  genplot_time = FALSE,
  genplot_uncertainty = FALSE,
  genplot_uncertainty_wt = FALSE,
  keep_editable=FALSE,
  palette_name="rainbow",
  color_brewer= "grDevices"
)
# }
```
