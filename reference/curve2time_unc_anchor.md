# Convert the re-tracked curve results to a depth time space with uncertainty

Converts the re-tracked curve results from
[`retrack_wt_MC`](https://stratigraphy.github.io/WaverideR/reference/retrack_wt_MC.md)
function to a depth time space using an anchor date while also taking
into account the uncertainty of the tracked astronomical cycle

## Usage

``` r
curve2time_unc_anchor(
  age_constraint = NULL,
  tracked_cycle_curve = NULL,
  tracked_cycle_period = NULL,
  tracked_cycle_period_unc = NULL,
  tracked_cycle_period_unc_dist = "n",
  n_simulations = 20,
  gap_constraints = NULL,
  proxy_data = NULL,
  cycles_check = NULL,
  uncer_cycles_check = NULL,
  max_runs = 1000,
  run_multicore = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  genplot = FALSE,
  keep_nr = 2,
  keep_all_time_curves = FALSE,
  dj = 1/200,
  lowerPeriod = 1,
  upperPeriod = 4600,
  omega_nr = 6,
  seed_nr = 1337,
  dir = TRUE
)
```

## Arguments

- age_constraint:

  Age constraints (anchor points, e.g. radiometric dates) for the
  modelling run: a data frame with 7 columns and one row per constraint.
  Column 1: ID name (character); column 2: age (numeric, usually in kyr,
  in the unit of the output ages); column 3: age uncertainty in the same
  unit (the standard deviation for a normal distribution, the full width
  of the range for a uniform one); column 4: distribution of the age,
  `"n"` (normal) or `"u"` (uniform); column 5: depth of the constraint
  (in the depth unit of `tracked_cycle_curve`); column 6: depth
  uncertainty (full width of the range for a uniform distribution);
  column 7: distribution of the depth, `"u"` (uniform), `"n"` (normal)
  or `"t"` (trapezoidal, with column 5 then holding the string
  `"min mode1 mode2 max"`). Required (`Default=NULL` gives an error).

- tracked_cycle_curve:

  Re-tracked cycle with its uncertainty (matrix or data frame, three
  columns): depth in column 1 (metres), the mean tracked frequency in
  cycles per metre in column 2 and its standard deviation (1 sigma) in
  column 3, e.g. the output of
  [`retrack_wt_MC`](https://stratigraphy.github.io/WaverideR/reference/retrack_wt_MC.md).
  Rows with `NA` are dropped and the curve is interpolated to its median
  depth spacing. Required (`Default=NULL` gives an error).

- tracked_cycle_period:

  Duration of the tracked cycle in kyr (single numeric value, \> 0, e.g.
  405 for the long eccentricity cycle), used to turn the tracked period
  into a sedimentation rate with
  [`curve2time`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md).
  Required (`Default=NULL` gives an error).

- tracked_cycle_period_unc:

  Uncertainty of `tracked_cycle_period` in kyr (single numeric value,
  \>= 0): the standard deviation when
  `tracked_cycle_period_unc_dist="n"`, the half-width of the range when
  `"u"`. Required (`Default=NULL` gives an error).

- tracked_cycle_period_unc_dist:

  Distribution from which the duration of the tracked cycle is drawn in
  every simulation (character): `"n"` normal or `"u"` uniform.
  `Default="n"`.

- n_simulations:

  Number of anchored age models to simulate (integer, \> 0).
  `Default=20`.

- gap_constraints:

  Hiatuses (gaps) to insert in the age models: a data frame with 4
  columns and one row per gap. Column 1: duration of the gap (kyr);
  column 2: its uncertainty (kyr; the standard deviation for a normal
  distribution, the half-width of the range for a uniform one; negative
  draws are set to 0); column 3: depth of the gap (depth unit of
  `tracked_cycle_curve`), below which all ages are shifted by the drawn
  duration; column 4: distribution, `"n"` (normal) or `"u"` (uniform).
  `NULL` inserts no gaps. `Default=NULL`.

- proxy_data:

  Proxy records used to check that the astronomical cycles are preserved
  after tuning (list of matrices or data frames, each with depth in
  column 1, in the unit of `tracked_cycle_curve`, and the proxy in
  column 2); in every run one record is picked at random, tuned with the
  simulated age model and analysed with
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md).
  Only used when `run_multicore=TRUE`. `Default=NULL`.

- cycles_check:

  Periods in kyr of the astronomical cycles (numeric vector) that must
  appear as peaks in the average wavelet spectrum of the tuned
  `proxy_data` for a simulation to be accepted. Only used when
  `run_multicore=TRUE`. `Default=NULL`.

- uncer_cycles_check:

  Tolerance of each period in `cycles_check`, in kyr (numeric vector of
  the same length): a peak counts if it lies within
  `cycles_check +/- uncer_cycles_check`. Only used when
  `run_multicore=TRUE`. `Default=NULL`.

- max_runs:

  Number of rejected simulations (in which the anchored astronomical
  ages do not match the age constraints) after which the age constraint
  that disagrees most often is dropped (integer, \> 0). `Default=1000`.

- run_multicore:

  Run the simulations in parallel on several cores (logical, TRUE or
  FALSE). Note that only the parallel version also checks the tuned
  proxy data for the cycles in `cycles_check`. `Default=FALSE`.

- n_cores:

  Number of cores used when `run_multicore=TRUE` (integer); `NULL` uses
  an automatic, safe number (physical cores minus one, at most 4, and
  limited by the free memory), which can be changed with
  `options(waverider.n_cores = ...)`. Ignored when
  `run_multicore=FALSE`. Results do not depend on the number of cores.
  `Default=NULL`.

- verbose:

  Show a progress bar (logical, TRUE or FALSE). `Default=FALSE`.

- genplot:

  Plot the age constraints (with their 2-sigma range) and the mean
  simulated age model with its 2-sigma envelope against depth (logical,
  TRUE or FALSE). `Default=FALSE`.

- keep_nr:

  Minimum number of age constraints that are kept when constraints are
  dropped after `max_runs` rejections (integer, \>= 1). `Default=2`.

- keep_all_time_curves:

  Also return all generated age curves, including the ones rejected
  during the modelling run, as the last element of the output list
  (logical, TRUE or FALSE). `Default=FALSE`.

- dj:

  Spacing of the periods of the wavelet transform of the tuned proxy
  data, in octaves (numeric, \> 0; e.g. `1/200` gives 200 periods per
  doubling); smaller values give a finer period axis but take longer.
  Passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  only used when `run_multicore=TRUE`. `Default=1/200`.

- lowerPeriod:

  Shortest period of the wavelet transform of the tuned proxy data, in
  kyr (numeric, \> 0; at least twice the sample spacing in time). Passed
  to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  only used when `run_multicore=TRUE`. `Default=1`.

- upperPeriod:

  Longest period of the wavelet transform of the tuned proxy data, in
  kyr (numeric, larger than `lowerPeriod`). Passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  only used when `run_multicore=TRUE`. `Default=4600`.

- omega_nr:

  Central angular frequency (omega_0) of the Morlet wavelet used on the
  tuned proxy data (numeric, \> 0; larger values give a better period
  resolution and a poorer depth resolution). Passed to
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  only used when `run_multicore=TRUE`. `Default=6`.

- seed_nr:

  Seed of the random number generator (integer), so that the Monte Carlo
  simulations can be reproduced. `Default=1337`.

- dir:

  Direction of time with depth (logical): `TRUE` time increases with
  depth (younger at the top), `FALSE` time decreases with depth (the age
  model is reversed). `Default=TRUE`.

## Value

The output is a list of 3 or 4 elements if keep_all_time_curves is set
to TRUE then the list consist of the x-axis, all the fitted curves in a
matrix format, the astrochronologically fitted age of the anchor, all
the generated depth time curves if keep_all_time_curves is set to TRUE
then the list consists of the x-axis, all the fitted curves in a matrix
format and the astrochronologically fitted age of the anchor If
`genplot=TRUE` then 3 plots stacked on top of each other will be
plotted. Plot 1: the original data set. Plot 2: the depth time plot.
Plot 3: the data set in the time domain. \#'

## References

Routines for astrochronologic testing, astronomical time scale
construction, and time series analysis
\<doi:10.1016/j.earscirev.2018.11.015\>

## Author

Part of the code is based on the
[sedrate2time](https://rdrr.io/pkg/astrochron/man/sedrate2time.html)
function of the 'astrochron' R package

## Examples

``` r
if (FALSE) { # \dontrun{


Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
Bisciaro_al <-
 astrochron::sortNave(Bisciaro_al, verbose = FALSE, genplot = FALSE)
Bisciaro_al <-
 astrochron::linterp(Bisciaro_al,
                     dt = 0.01,
                     verbose = FALSE,
                     genplot = FALSE)
Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2, ]

Bisciaro_al_wt <-
 analyse_wavelet(
   data = Bisciaro_al,
   dj = 1 / 200 ,
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
Bisciaro_ca <-
 astrochron::sortNave(Bisciaro_ca, verbose = FALSE, genplot = FALSE)
Bisciaro_ca <-
 astrochron::linterp(Bisciaro_ca,
                     dt = 0.01,
                     verbose = FALSE,
                     genplot = FALSE)
Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2, ]

Bisciaro_ca_wt <-
 analyse_wavelet(
   data = Bisciaro_ca,
   dj = 1 / 200 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )


#
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
#     keep_editable = FALSE
#   )


Bisciaro_sial <- Bisciaro_XRF[, c(1, 64)]
Bisciaro_sial <-
 astrochron::sortNave(Bisciaro_sial, verbose = FALSE, genplot = FALSE)
Bisciaro_sial <-
 astrochron::linterp(Bisciaro_sial,
                     dt = 0.01,
                     verbose = FALSE,
                     genplot = FALSE)
Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2, ]

Bisciaro_sial_wt <-
 analyse_wavelet(
   data = Bisciaro_sial,
   dj = 1 / 200 ,
   lowerPeriod = 0.01,
   upperPeriod = 50,
   verbose = FALSE,
   omega_nr = 8
 )

#Bisciaro_sial_wt_track <-
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


Bisciaro_Mn <- Bisciaro_XRF[, c(1, 46)]
Bisciaro_Mn <-
 astrochron::sortNave(Bisciaro_Mn, verbose = FALSE, genplot = FALSE)
Bisciaro_Mn <-
 astrochron::linterp(Bisciaro_Mn,
                     dt = 0.01,
                     verbose = FALSE,
                     genplot = FALSE)
Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2, ]

Bisciaro_Mn_wt <-
 analyse_wavelet(
   data = Bisciaro_Mn,
   dj = 1 / 200 ,
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

Bisciaro_Mg <- Bisciaro_XRF[, c(1, 71)]
Bisciaro_Mg <-
 astrochron::sortNave(Bisciaro_Mg, verbose = FALSE, genplot = FALSE)
Bisciaro_Mg <-
 astrochron::linterp(Bisciaro_Mg,
                     dt = 0.01,
                     verbose = FALSE,
                     genplot = FALSE)
Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2, ]

Bisciaro_Mg_wt <-
 analyse_wavelet(
   data = Bisciaro_Mg,
   dj = 1 / 200 ,
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
#     keep_editable = FALSE
#   )




wt_list_bisc <- list(Bisciaro_al_wt,
                    Bisciaro_ca_wt,
                    Bisciaro_sial_wt,
                    Bisciaro_Mn_wt,
                    Bisciaro_Mg_wt)


data_track_bisc <- cbind(
 Bisciaro_al_wt_track[, 2],
 Bisciaro_ca_wt_track[, 2],
Bisciaro_sial_wt_track[, 2],
 Bisciaro_Mn_wt_track[, 2],
 Bisciaro_Mg_wt_track[, 2]
)

x_axis_bisc <- Bisciaro_al_wt_track[, 1]

bisc_retrack <- retrack_wt_MC(
wt_list = wt_list_bisc,
 data_track = data_track_bisc,
 x_axis = x_axis_bisc,
 nr_simulations = 500,
 seed_nr = 1337,
 verbose = TRUE,
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
file_name = "TEST",
run_multicore = TRUE,
output = 5,
 n_imgs = 50,
 plot_horizontal = TRUE,
 empty_folder = FALSE
)

proxy_list_bisc <- list(Bisciaro_al,
                    Bisciaro_ca,
                    Bisciaro_sial,
                    Bisciaro_Mn,
                    Bisciaro_Mg)




id <- c("CCT18_322", "CCT18_315", "CCT18_311")
ages <- c(20158, 20575, 20857)
ageSds <- c(28, 40, 34)
ages_unc_dist <- c("n", "n", "n")
position <- c(13.3, 7.25, 3.2)
anchor_thick <- c(0.2, 0.1, 0.1)
anchor_thick_unc_dist <- c("u", "u", "u")

ash_Bisc <-
 as.data.frame(
   cbind(
     id,
     ages,
     ageSds,
     ages_unc_dist,
     position,
     anchor_thick,
     anchor_thick_unc_dist
   )
 )

gap_dur = c(10, 20)
gap_unc = c(3, 10)
gap_depth = c(2.5, 9)
gap_unc_dist = c("n", "n")


gap_constraints_Bisc <-
 as.data.frame(cbind(gap_dur, gap_unc, gap_depth, gap_unc_dist))

cycles_checks <- c(110,40,22)
uncer_cycles_checks <- c(20,5,7)

curve2time_unc_anchor_res <-
 curve2time_unc_anchor(
 age_constraint =  ash_Bisc,
  tracked_cycle_curve =  bisc_retrack,
  tracked_cycle_period = 110,
   tracked_cycle_period_unc = ((135 - 110) + (110 - 95)) / 2,
  tracked_cycle_period_unc_dist = "n",
   n_simulations = 20,
   gap_constraints = gap_constraints_Bisc,
   proxy_data = proxy_list_bisc,
   cycles_check = NULL,
   uncer_cycles_check = NULL,
   cycles_check = cycles_checks,
   uncer_cycles_check = uncer_cycles_checks,
  max_runs = 1000,
   run_multicore = FALSE,
   verbose = FALSE,
   genplot = FALSE,
   keep_nr = 2,
   keep_all_time_curves = FALSE,
   dj = 1/200,
   lowerPeriod =1,
   upperPeriod =2500,
   omega_nr = 6,
   seed_nr=1337,
   dir=TRUE
 )

} # }
```
