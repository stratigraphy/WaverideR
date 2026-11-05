# Detect and filter out all minima in a signal

The `min_detect` function is used to detect and filter out local minima
in a sinusoidal signal

## Usage

``` r
min_detect(data = NULL, pts = 3)
```

## Arguments

- data:

  Signal in which the local extremes are detected: a two-column matrix
  or data frame with depth or time in column 1 (in increasing order) and
  the proxy value in column 2. Required (`Default=NULL`).

- pts:

  Number of points to each side that a point is compared with (integer,
  \>= 2): a point is a local minimum when it is above (maximum) or below
  (minimum) the value `pts` points further on and the value `pts - 1`
  points back. Of consecutive extremes of the same type only the most
  extreme one is kept. Larger values give more robust extremes but may
  miss minor ones. `Default=3`.

## Value

\#Returns a matrix with 2 columns first column is depth/time the second
column are local minima values

## Examples

``` r
#Example in which the ~210yr de Vries cycle is extracted from the Total Solar
#Irradiance data set of Steinhilber et al., (2012)
#after which all minima are extracted

TSI_wt <-
analyse_wavelet(
data = TSI,
dj = 1/200,
lowerPeriod = 16,
upperPeriod = 8192,
   verbose = FALSE,
   omega_nr = 6
 )

de_Vries_cycle <- extract_signal_stable(wavelet=TSI_wt,
cycle=210,
period_up =1.25,
period_down = 0.75,
add_mean=TRUE,
plot_residual=FALSE)


min_de_Vries_cycle <- min_detect(de_Vries_cycle)
```
