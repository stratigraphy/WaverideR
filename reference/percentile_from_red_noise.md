# Calculate average spectral power from red noise curves for a given percentile

The `percentile_from_red_noise` function is used to generate and average
spectral power curve based on a set percentile based. To generate the
percentile curve the results of the
[`model_red_noise_wt`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md)
function are used.

## Usage

``` r
percentile_from_red_noise(red_noise = NULL, wavelet = NULL, percentile = NULL)
```

## Arguments

- red_noise:

  Average spectral power curves of red-noise simulations, the output of
  [`model_red_noise_wt`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md)
  (numeric matrix, one row per period of `wavelet` and one column per
  simulation). Required (`Default=NULL`).

- wavelet:

  Wavelet object, the output of
  [`analyse_wavelet`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
  the same object that was given to `model_red_noise_wt`, its periods
  form the first output column. Required (`Default=NULL`).

- percentile:

  Percentile of the red-noise power that is returned for every period
  (numeric, between 0 and 1, e.g. 0.95); computed from a normal
  distribution with the mean and standard deviation of the simulations.
  Required (`Default=NULL`).

## Value

Returns a matrix with 2 columns.  
The first column is the period (m).  
The second column is the spectral power at percentile x based on  
the red noise modelling runs.  

## Examples

``` r
# \donttest{
#'#generate average spectral power curves based on red noise curves which are
# based on the magnetic susceptibility record of the Sullivan core of Pas et al., (2018)

mag_wt <- analyse_wavelet(data = mag,
dj = 1/100,
lowerPeriod = 0.1,
upperPeriod = 254,
verbose = FALSE,
omega_nr = 10)

#increase n_simulations to better define the red noise spectral power curve
mag_wt_red_noise <- model_red_noise_wt(wavelet=mag_wt,
n_simulations=10, # Increase number for a better constrained result
run_multicore=FALSE,
verbose=FALSE)


prob_curve <- percentile_from_red_noise(
red_noise = mag_wt_red_noise,
wavelet = mag_wt,
percentile = 0.9)

# }
```
