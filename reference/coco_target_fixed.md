# Fixed-weight astronomical target spectrum (Fixed-target COCO)

Builds the fixed-relative-weight astronomical target used by Acycle's
Fixed-target COCO: nine (or other) orbital sine terms with the
historical amplitudes 1.0 (eccentricity), 0.8 (obliquity) and 0.6
(precession) (cocoFixedTargetWeights.m), sampled on the data grid mapped
to time at sedimentation rate `sr`, linearly detrended, and turned into
a rectangular-window periodogram of length max(pad, n) that is linearly
interpolated at the requested temporal frequencies (0 outside the target
band). Unlike the adaptive target, the amplitudes never depend on the
data spectrum. Translated from Acycle v3.1 freq2targetFixed.m and
cocoFixedTargetWeights.m (Mingsong Li).

## Usage

``` r
coco_target_fixed(
  data = NULL,
  pad = 5000,
  orbit9 = c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.682, 22.3758,
    18.9519),
  freq = NULL,
  sr = NULL
)
```

## Arguments

- data:

  Record whose sampling grid the target is synthesised on (matrix or
  data frame; only the depth column 1, in metres, is used): the target
  has as many samples as `data`, at the median depth spacing converted
  to time with `sr`. `NULL` (or `sr=NULL`) uses `5 * max(orbit9)`
  samples spaced 1 kyr apart. `Default=NULL`.

- pad:

  Length (number of points) of the zero-padded periodogram of the target
  (integer, \> 0; at least the number of samples is used); a larger
  value gives a finer frequency grid. `NULL` or invalid values are
  replaced by 5000. `Default=5000`.

- orbit9:

  Astronomical periods of the target sinusoids, in kyr (numeric vector,
  \> 0; usually nine values); their amplitudes are fixed by
  [`coco_fixed_target_weights`](https://stratigraphy.github.io/WaverideR/reference/coco_fixed_target_weights.md)
  (1.0 eccentricity, 0.8 obliquity, 0.6 precession). The default is the
  present-day set.
  `Default=c(405.6912, 130.6979, 123.8532, 98.8517, 94.8856, 40.9897, 23.6820, 22.3758, 18.9519)`.

- freq:

  Frequencies at which the target power is returned, in cycles/kyr
  (numeric vector); values outside the periodogram grid get 0. `NULL`
  returns an empty vector. `Default=NULL`.

- sr:

  Sedimentation rate used to convert the depth grid of `data` to time,
  in cm/kyr (single numeric value, \> 0); `NULL` uses the 1-kyr grid
  described under `data`. `Default=NULL`.

## Value

Numeric vector of target power at `freq`.

## References

Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking variable
sedimentation rates and astronomical forcing in Phanerozoic paleoclimate
proxy series with evolutionary correlation coefficients and hypothesis
testing. Earth and Planetary Science Letters 501, 165-179.

## Author

Michiel Arts (R translation); original MATLAB code: Mingsong Li

## Examples

``` r
# \donttest{
depth <- seq(0, 20, by = 0.05)
f <- seq(0, 0.06, by = 0.0005)
la <- coco_target_fixed(cbind(depth, 0), pad = 5000, freq = f, sr = 5)
plot(f, la, type = "l", xlab = "Frequency (cycles/kyr)", ylab = "Power")

# }
```
