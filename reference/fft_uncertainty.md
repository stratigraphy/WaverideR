# Calculate the uncertainty associated with a tracked cycle based on the frequency uncertainty of the FFT

Fourier counterpart of
[`wavelet_uncertainty`](https://stratigraphy.github.io/WaverideR/reference/wavelet_uncertainty.md).
A tracked cycle (depth, period) gets an uncertainty in frequency that
equals the bin spacing of an FFT of the analysed record (\\1/T\\, the
practical resolution) or the Heisenberg-Gabor limit (\\1/(4 \pi T)\\),
see
[`fft_frequency_uncertainty`](https://stratigraphy.github.io/WaverideR/reference/fft_frequency_uncertainty.md).
In contrast to the wavelet the uncertainty does not depend on the
frequency, only on the length of the record (or of the window of an
evolutive spectrum). The uncertain frequencies are converted to curves
of period (sedimentation rate) and with
[`curve2time`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md)
to uncertain time models.

## Usage

``` r
fft_uncertainty(
  tracked_cycle = NULL,
  period_of_tracked_cycle = NULL,
  signal_length = NULL,
  type = c("practical", "gabor"),
  multi = 1,
  verbose = FALSE,
  genplot_time = FALSE,
  genplot_uncertainty = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- tracked_cycle:

  Curve of the tracked cycle, e.g. the output of
  [`track_period_wavelet`](https://stratigraphy.github.io/WaverideR/reference/track_period_wavelet.md):
  a matrix or data frame with depth in column 1 and the period of the
  tracked cycle (in the same depth unit, e.g. metres) in column 2;
  further columns are ignored. Required (`Default=NULL`).

- period_of_tracked_cycle:

  Duration of one cycle of the tracked cycle in time (numeric, \> 0,
  e.g. 405 for the 405-kyr eccentricity cycle); passed to
  [`curve2time`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md)
  and sets the time unit of the output. Required (`Default=NULL`).

- signal_length:

  Length of the record (or of the window of an evolutive spectrum) that
  sets the frequency uncertainty, in the depth unit of `tracked_cycle`
  (numeric, \> 0). `NULL` uses the depth range of `tracked_cycle`.
  `Default=NULL`.

- type:

  Definition of the frequency uncertainty (character): `"practical"`
  uses the FFT bin spacing \\1/T\\; `"gabor"` uses the Heisenberg-Gabor
  limit \\1/(4 \pi T)\\, which is about 12.6 times smaller.
  `Default=c("practical", "gabor")`, i.e. `"practical"`.

- multi:

  Multiple of the frequency uncertainty that is added to and subtracted
  from the frequency of the tracked cycle (numeric, \> 0; e.g. 2 for
  twice the uncertainty). `Default=1`.

- verbose:

  Print the total duration, its uncertainty and the frequency
  uncertainty (logical, TRUE or FALSE). `Default=FALSE`.

- genplot_time:

  Plot the time-depth curves of the mean, lower and upper uncertainty
  (logical, TRUE or FALSE). `Default=FALSE`.

- genplot_uncertainty:

  Plot the period-depth curves of the tracked cycle with its lower and
  upper uncertainty (logical, TRUE or FALSE). `Default=FALSE`.

- keep_editable:

  Keep the graphical parameters (`par`) of the plot after the function
  returns, so that features can be added to it (logical, TRUE or FALSE);
  `FALSE` restores the previous settings. `Default=FALSE`.

## Value

A data frame with the columns "depth", "period", "frequency",
"uncertainty in frequency", "uncertainty in frequency x_times", "time
mean", "time plus x_times" and "time min x_times" (the same layout as
[`wavelet_uncertainty`](https://stratigraphy.github.io/WaverideR/reference/wavelet_uncertainty.md),
with the Fourier uncertainty), and the attribute `"type"`.

## Examples

``` r
tr <- data.frame(depth = seq(0, 100, by = 1), period = 10 + 0 * seq(0, 100, by = 1))
unc <- fft_uncertainty(tr, period_of_tracked_cycle = 100, type = "practical")
tail(unc, 2)
#>     depth period frequency uncertainty in frequency
#> 100    99     10       0.1                     0.01
#> 101   100     10       0.1                     0.01
#>     uncertainty in frequency x_times time mean time min x_times
#> 100                             0.01       990              891
#> 101                             0.01      1000              900
#>     time plus x_times
#> 100              1089
#> 101              1100
```
