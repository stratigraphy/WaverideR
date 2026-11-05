# Re-track cycles using a Monte-Carlo simulation

When analysing multi-proxy records an age-model can be created for each
proxy. These age-models can be in general agreement but might also
indicate conflicting deposition rates. Picking one age-model out of the
all multi-proxy age-models and stating that, that age-model is the best
overlooks the information contained within the other proxies and hence a
degree of error remains the age-model exists.To combine the multiple
age-models all the age models can be averaged out and the uncertainty
can be calculated by means of the standard deviation. The result is an
age-model which takes into account all the age-models from the proxy
records. The averaged out age-model does not take into account any small
user errors during the creation of the individual age-models nor does
the averaging take into account the differences between the age-models
and how the initial age-model of a certain proxy might be off in certain
intervals. the retrack_wt_MC mitigates these problems by re-tracking
periods of astronomical cycles in the wavelet spectra. The re-tracking
is based on the information provided by the age-models constructed from
the different proxy records. First a synthetic tracked curve is created
by adding up fractions (0-1) of the tracked periods of the different
proxy records. This synthetic curve is then used to re-track the
period/spectral peaks of an astronomical cycle in a randomly select
wavelet scalogram. This process is repeated x times. The result x
tracked curves which take into account all the original age-models. From
the re-tracked curves one can calculate the mean period and the standard
deviation. The resulting standard deviation is a good indicator of the
quality of the imprint of of astronomical cycles in the proxy records. A
small standard deviation indicates that given the input of the different
tracked cycles similar periods keep on being tracked indicating the an
astronomical is well recorded in the proxy records and as such the
age-model is very reliable in set interval. A high standard deviation on
the other hand means that the tracking results in vastly different
periods of the tracked astronomical cycle, as such the quality of the
imprint of the astronomical cycle proxy records is poor and hence the
age-model is less-reliable in this interval.

## Usage

``` r
retrack_wt_MC(
  wt_list = NULL,
  data_track = NULL,
  x_axis = NULL,
  smoothing = c("auto"),
  nr_simulations = 50,
  seed_nr = 1337,
  verbose = FALSE,
  genplot = FALSE,
  keep_editable = FALSE,
  create_GIF = FALSE,
  plot_GIF = FALSE,
  width_plt = 600,
  height_plt = 450,
  period_up = 1.5,
  period_down = 0.5,
  plot.COI = TRUE,
  n.levels = 100,
  palette_name = "viridis",
  color_brewer = "viridis",
  periodlab = "Period (metres)",
  x_lab = "depth (metres)",
  add_avg = FALSE,
  time_dir = TRUE,
  file_name = NULL,
  run_multicore = FALSE,
  n_cores = NULL,
  output = 1,
  n_imgs = 50,
  plot_horizontal = TRUE,
  empty_folder = FALSE,
  colour_brewer = NULL
)
```

## Arguments

- wt_list:

  Wavelet scalograms of the different proxy records (required): a `list`
  of objects created with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
  one per proxy, in the same order as the columns of `data_track`. In
  every simulation one of them is picked at random and the cycle is
  re-tracked in it. `Default=NULL`.

- data_track:

  Tracked periods of the cycle in each proxy record (required): numeric
  matrix with one column per proxy (same order as `wt_list`) and one row
  per value of `x_axis`; contains only the period values, not the depth
  axis. Build it with [`cbind`](https://rdrr.io/r/base/cbind.html) from
  tracked curves that were interpolated to the same depths.
  `Default=NULL`.

- x_axis:

  Depth/time values of the rows of `data_track` (numeric vector, evenly
  spaced, same length as `nrow(data_track)`). Required. `Default=NULL`.

- smoothing:

  Smoothing of every re-tracked curve (character vector of length 1 or
  2):

  - `"auto"`: loess smoothing with an automatically selected span
    ([`loess_auto`](https://stratigraphy.github.io/WaverideR/reference/loess_auto.md)),

  - `c("loess", span)`: loess smoothing with the given span (0 to 1;
    `smooth` of
    [`astrochron::noLow`](https://rdrr.io/pkg/astrochron/man/noLow.html)),

  - `c("window", length)`: moving average with the given window length
    in depth/time units (`win` of
    [`astrochron::mwStats`](https://rdrr.io/pkg/astrochron/man/mwStats.html)).

  Any other value applies no smoothing. `Default=c("auto")`.

- nr_simulations:

  Number of Monte Carlo simulations (re-tracked curves) (integer, \> 0).
  `Default=50`.

- seed_nr:

  Seed of the random number generator (integer), so that the random
  choice of scalogram and of the weights of the synthetic curve can be
  reproduced. `Default=1337`.

- verbose:

  Show the progress of the simulations (logical, TRUE or FALSE).
  `Default=FALSE`.

- genplot:

  Plot the mean re-tracked period (green) with the mean minus (red) and
  plus (blue) one standard deviation of the frequency (logical).
  `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters of the plot after the function returns,
  so that extra elements can be added (logical); only used with
  `genplot = TRUE`. `Default=FALSE`.

- create_GIF:

  Save one image per simulation, showing the re-tracked curve on the
  selected scalogram (drawn with
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)),
  in the folder `file_name` and join them into an animated GIF (logical;
  needs the magick package). `Default=FALSE`.

- plot_GIF:

  Display the GIF after it is created (logical); only used with
  `create_GIF = TRUE`. `Default=FALSE`.

- width_plt:

  Width of the GIF images in pixels (integer). `Default=600`.

- height_plt:

  Height of the GIF images in pixels (integer). `Default=450`.

- period_up:

  Factor (numeric, \> 1) by which the synthetic tracked curve is
  multiplied to get the upper boundary of the search window for the
  nearest spectral peak (passed to
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md)).
  If no peak is found within the window the interpolated value is used.
  `Default=1.5`.

- period_down:

  Factor (numeric, between 0 and 1) by which the synthetic tracked curve
  is multiplied to get the lower boundary of the search window for the
  nearest spectral peak (passed to
  [`completed_series`](https://stratigraphy.github.io/WaverideR/reference/completed_series.md)).
  `Default=0.5`.

- plot.COI:

  Draw the cone of influence on the scalograms of the GIF images
  (logical); only used with `create_GIF = TRUE`. `Default=TRUE`.

- n.levels:

  Number of colour levels of the scalograms of the GIF images (integer,
  \> 1). `Default=100`.

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

- periodlab:

  Label of the period axis of the plots (character).
  `Default="Period (metres)"`.

- x_lab:

  Label of the depth/time axis of the plots (character).
  `Default="depth (metres)"`.

- add_avg:

  Draw the average wavelet power next to the scalograms of the GIF
  images (logical; passed to
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md));
  only used with `create_GIF = TRUE`. `Default=FALSE`.

- time_dir:

  Direction of the depth/time axis in the GIF images (logical): `TRUE`
  if time increases with increasing depth/time values (e.g. a core that
  gets older downwards), `FALSE` if it decreases (e.g. a section log
  where 0 m is the base); only used with `create_GIF = TRUE`.
  `Default=TRUE`.

- file_name:

  Name (character) of the folder created in the working directory for
  the GIF images, also used as the base name of the images
  (`file_name_<simulation number>.jpeg`) and of the GIF
  (`file_name.gif`). Required when `create_GIF = TRUE`. `Default=NULL`.

- run_multicore:

  Run the Monte Carlo simulations in parallel on several cores
  (logical). `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore = TRUE` (integer, or `NULL`).
  `NULL` uses an automatic, safe number: physical cores minus one, at
  most 4, and limited by the free memory; set
  `options(waverider.n_cores = ...)` to change it. A value \> 1 also
  switches on parallel computing. Results do not depend on the number of
  cores. `Default=NULL`.

- output:

  What is returned (integer, 1 to 7):

  - 1: list of (1) a matrix with the x-axis, the mean tracked frequency
    and its standard deviation (columns `depth`, `mean_freq`, `sd`), (2)
    a matrix with all the re-tracked periods (one column per simulation)
    and (3) the GIF (`NULL` without `create_GIF = TRUE`),

  - 2: list of the mean/sd matrix and the matrix of re-tracked periods,

  - 3: list of the mean/sd matrix and the GIF,

  - 4: list of the matrix of re-tracked periods and the GIF,

  - 5: the mean/sd matrix only,

  - 6: the matrix of re-tracked periods only,

  - 7: the GIF only.

  `Default=1`.

- n_imgs:

  Number of images (simulations) put in the GIF (integer, \> 0; at most
  `nr_simulations`). Many images take long to process and give a large
  GIF. `Default=50`.

- plot_horizontal:

  Orientation of the scalograms of the GIF images (logical): `TRUE` puts
  depth/time on the x-axis, `FALSE` on the y-axis. `Default=TRUE`.

- empty_folder:

  Delete all files in the folder `file_name` before the new images are
  written (logical); only used with `create_GIF = TRUE`.
  `Default=FALSE`.

- colour_brewer:

  British spelling of `color_brewer`, with the same options; if given
  (not `NULL`) it overrides `color_brewer`. `Default=NULL`.

## Value

The output depends on the output setting If genplot = TRUE a plot will
be generated in which the mean period and standard deviation is plotted
if plot_GIF = TRUE a GIF with n number of n_imgs will be plotted in
which the retraced curve is plotted in a wavelet scalogram If output =
1, output is a list which contain 3 objects. object 1 is a matrix with
the x-axis and the mean tracked frequency and standard deviation. object
2 is a matrix with all the tracked periods. Object 3 is a GIF in which
all the tracked periods are plotted. If output = 2, output is a list
which contain 2 objects. object 1 is a matrix with the x-axis and the
mean tracked frequency and standard deviation. object 2 is a matrix with
all the tracked periods. If output = 3, output is a list which contain 2
objects. object 1 is a matrix with the x-axis and the mean tracked
frequency and standard deviation. Object 2 is a GIF in which all the
tracked periods are plotted. If output = 4, output is a list which
contain 3 objects. Object 1 is a matrix with all the tracked periods.
Object 2 is a GIF in which all the tracked periods are plotted. If
output = 4 output is a list which contain 3 objects. Object 1 is a
matrix with all the tracked periods. Object 2 is a GIF in which all the
tracked periods are plotted. If output = 5 a matrix with the x-axis and
the mean tracked period and standard deviation is returned. If output =
6, a matrix with all the tracked periods is returned. If output = 7, a
GIF in which all the tracked periods are plotted is returned

## Examples

``` r
# \donttest{
# Re-track the 110kyr eccentricity cycle in the wavelet scalogram
# from the XRF record of the Bisciaro data set of Arts (2014)

Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
Bisciaro_al <- astrochron::sortNave(Bisciaro_al,verbose=FALSE,genplot=FALSE)
#>  * Duplicates found
#>  * Duplicates values will be averaged.
#>  * Number of samples post-averaging= 287 
Bisciaro_al <- astrochron::linterp(Bisciaro_al, dt = 0.01,verbose=FALSE,genplot=FALSE)
Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2, ]

Bisciaro_al_wt <-
 analyse_wavelet(
   data = Bisciaro_al,
   dj = 1 /50 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

# Bisciaro_al_wt_track <-
#   track_period_wavelet(
#     astro_cycle = 110,
#     wavelet = Bisciaro_al_wt,
#     n.levels = 100,
#     periodlab = "Period (metres)",
#     x_lab = "depth (metres)"
#   )
#
# Bisciaro_al_wt_track <- completed_series(
#   wavelet = Bisciaro_al_wt,
#   tracked_curve = Bisciaro_al_wt_track,
#   period_up = 1.2,
#   period_down = 0.8,
#   extrapolate = TRUE,
#   genplot = FALSE,
#   keep_editable = FALSE
# )
#
# Bisciaro_al_wt_track <-
#   loess_auto(
#     time_series = Bisciaro_al_wt_track,
#     genplot = FALSE,
#     print_span = FALSE,
#     keep_editable = FALSE
#   )

Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
Bisciaro_ca <- astrochron::sortNave(Bisciaro_ca,verbose=FALSE,genplot=FALSE)
#>  * Duplicates found
#>  * Duplicates values will be averaged.
#>  * Number of samples post-averaging= 287 
Bisciaro_ca <- astrochron::linterp(Bisciaro_ca, dt = 0.01,verbose=FALSE,genplot=FALSE)
Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2, ]

Bisciaro_ca_wt <-
 analyse_wavelet(
   data = Bisciaro_ca,
   dj = 1 /50 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

# Bisciaro_ca_wt_track <-
#   track_period_wavelet(
#     astro_cycle = 110,
#     wavelet = Bisciaro_ca_wt,
#     n.levels = 100,
#     periodlab = "Period (metres)",
#     x_lab = "depth (metres)"
#   )
#
# Bisciaro_ca_wt_track <- completed_series(
#   wavelet = Bisciaro_ca_wt,
#   tracked_curve = Bisciaro_ca_wt_track,
#   period_up = 1.2,
#   period_down = 0.8,
#   extrapolate = TRUE,
#   genplot = FALSE,
#   keep_editable = FALSE
# )
#
# Bisciaro_ca_wt_track <-
#   loess_auto(
#     time_series = Bisciaro_ca_wt_track,
#     genplot = FALSE,
#     print_span = FALSE,
#     keep_editable = FALSE)

Bisciaro_sial <- Bisciaro_XRF[,c(1,64)]
Bisciaro_sial <- astrochron::sortNave(Bisciaro_sial,verbose=FALSE,genplot=FALSE)
#>  * Duplicates found
#>  * Duplicates values will be averaged.
#>  * Number of samples post-averaging= 287 
Bisciaro_sial <- astrochron::linterp(Bisciaro_sial, dt = 0.01,verbose=FALSE,genplot=FALSE)
Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2, ]

Bisciaro_sial_wt <-
 analyse_wavelet(
   data = Bisciaro_sial,
   dj = 1 /50 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

# Bisciaro_sial_wt_track <-
#   track_period_wavelet(
#     astro_cycle = 110,
#     wavelet = Bisciaro_sial_wt,
#     n.levels = 100,
#     periodlab = "Period (metres)",
#     x_lab = "depth (metres)"
#   )
#
#
# Bisciaro_sial_wt_track <- completed_series(
#   wavelet = Bisciaro_sial_wt,
#   tracked_curve = Bisciaro_sial_wt_track,
#   period_up = 1.2,
#   period_down = 0.8,
#   extrapolate = TRUE,
#   genplot = FALSE,
#   keep_editable = FALSE
# )
#
# Bisciaro_sial_wt_track <-
#   loess_auto(
#     time_series = Bisciaro_sial_wt_track,
#     genplot = FALSE,
#     print_span = FALSE,
#     keep_editable = FALSE
#   )


Bisciaro_Mn <- Bisciaro_XRF[,c(1,46)]
Bisciaro_Mn <- astrochron::sortNave(Bisciaro_Mn,verbose=FALSE,genplot=FALSE)
#>  * Duplicates found
#>  * Duplicates values will be averaged.
#>  * Number of samples post-averaging= 287 
Bisciaro_Mn <- astrochron::linterp(Bisciaro_Mn, dt = 0.01,verbose=FALSE,genplot=FALSE)
Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2, ]

Bisciaro_Mn_wt <-
 analyse_wavelet(
   data = Bisciaro_Mn,
   dj = 1 /50 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

# Bisciaro_Mn_wt_track <-
#   track_period_wavelet(
#     astro_cycle = 110,
#     wavelet = Bisciaro_Mn_wt,
#     n.levels = 100,
#     periodlab = "Period (metres)",
#     x_lab = "depth (metres)"
#   )
#
#
# Bisciaro_Mn_wt_track <- completed_series(
#   wavelet = Bisciaro_Mn_wt,
#   tracked_curve = Bisciaro_Mn_wt_track,
#   period_up = 1.2,
#   period_down = 0.8,
#   extrapolate = TRUE,
#   genplot = FALSE,
#   keep_editable = FALSE
# )
# Bisciaro_Mn_wt_track <-
#   loess_auto(
#     time_series = Bisciaro_Mn_wt_track,
#     genplot = FALSE,
#     print_span = FALSE,
#     keep_editable = FALSE
#   )

Bisciaro_Mg <- Bisciaro_XRF[,c(1,71)]
Bisciaro_Mg <- astrochron::sortNave(Bisciaro_Mg,verbose=FALSE,genplot=FALSE)
#>  * Duplicates found
#>  * Duplicates values will be averaged.
#>  * Number of samples post-averaging= 287 
Bisciaro_Mg <- astrochron::linterp(Bisciaro_Mg, dt = 0.01,verbose=FALSE,genplot=FALSE)
Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2, ]

Bisciaro_Mg_wt <-
 analyse_wavelet(
   data = Bisciaro_Mg,
   dj = 1 /50 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

# Bisciaro_Mg_wt_track <-
#   track_period_wavelet(
#     astro_cycle = 110,
#     wavelet = Bisciaro_Mg_wt,
#     n.levels = 100,
#     periodlab = "Period (metres)",
#     x_lab = "depth (metres)"
#   )
#
#
# Bisciaro_Mg_wt_track <- completed_series(
#   wavelet = Bisciaro_Mg_wt,
#   tracked_curve = Bisciaro_Mg_wt_track,
#   period_up = 1.2,
#   period_down = 0.8,
#   extrapolate = TRUE,
#   genplot = FALSE,
#   keep_editable = FALSE
# )
#
# Bisciaro_Mg_wt_track <-
#   loess_auto(
#     time_series = Bisciaro_Mg_wt_track,
#     genplot = FALSE,
#     print_span = FALSE,
#     keep_editable = FALSE)




wt_list_bisc <- list(Bisciaro_al_wt,
               Bisciaro_ca_wt,
               Bisciaro_sial_wt,
               Bisciaro_Mn_wt,
               Bisciaro_Mg_wt)

#Instead of tracking, the tracked solution data sets Bisciaro_al_wt_track,
#Bisciaro_ca_wt_track, Bisciaro_sial_wt_track, Bisciaro_Mn_wt_track,
# Bisciaro_Mn_wt_track and Bisciaro_Mg_wt_track are used

data_track_bisc <- cbind(Bisciaro_al_wt_track[,2],
                     Bisciaro_ca_wt_track[,2],
                     Bisciaro_sial_wt_track[,2],
                     Bisciaro_Mn_wt_track[,2],
                     Bisciaro_Mg_wt_track[,2])

x_axis_bisc <- Bisciaro_al_wt_track[,1]


bisc_retrack <- retrack_wt_MC(wt_list = wt_list_bisc,
             data_track = data_track_bisc,
             x_axis = x_axis_bisc,
             nr_simulations = 20,
             seed_nr = 1337,
             verbose = FALSE,
             genplot = FALSE,
             keep_editable = FALSE,
             create_GIF = FALSE,
             plot_GIF = FALSE,
             width_plt =  600,
             height_plt = 450,
            period_up  =  1.5,
             period_down = 0.5,
             plot.COI = TRUE,
             n.levels = 100,
             palette_name = "rainbow",
             color_brewer = "grDevices",
             periodlab = "Period (metres)",
             x_lab = "depth (metres)",
             add_avg = FALSE,
             time_dir = TRUE,
             file_name = NULL,
             run_multicore = FALSE,
             output = 1,
             n_imgs = 50,
             plot_horizontal = TRUE,
             empty_folder = FALSE)


# }
```
