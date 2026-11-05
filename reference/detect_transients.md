# Transient (impulsive) features of a record by transient extraction

Finds the positions in a record where energy over many periods arrives
at the same instant: the signature of an impulsive feature such as a
hiatus (a phase jump of all cycles), an erosion surface, a turbidite or
a spike. The wavelet coefficients are reassigned along the record to
their group delay (the time-reassigned synchrosqueezing /
transient-extracting transform, He et al. 2019; Yu et al. 2023), only
the coefficients that lie on their own group delay are kept, and their
energy is summed over the periods. The sum is compared with a running
median and peaks above `threshold` times the median are reported, which
can be passed as `jump_depths` to
[`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
or
[`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md).

## Usage

``` r
detect_transients(
  data,
  lowerPeriod,
  upperPeriod,
  dj = 1/20,
  omega_nr = 6,
  threshold = 4,
  min_distance = NULL,
  window = NULL
)
```

## Arguments

- data:

  Record to scan (matrix or data frame): depth in column 1 and the proxy
  in column 2; incomplete rows are dropped, the rows sorted, the record
  linearly interpolated to its median spacing and standardised.
  Required, no default.

- lowerPeriod:

  Shortest period of the Morlet wavelet transform, in depth units
  (numeric, \> 0; at least twice the sample spacing). Required, no
  default.

- upperPeriod:

  Longest period of the transform, in depth units (numeric, larger than
  `lowerPeriod`); it also sets the defaults of `min_distance` and
  `window`. Required, no default.

- dj:

  Spacing of the periods of the transform, in octaves (numeric, \> 0;
  `1/20` gives 20 periods per doubling). `Default=1/20`.

- omega_nr:

  Central angular frequency (omega_0) of the Morlet wavelet (numeric, \>
  0); small values give short wavelets, which localise transients best.
  `Default=6`.

- threshold:

  Detection threshold (numeric, \> 1): positions where the reassigned
  (transient) energy exceeds this multiple of the running median of the
  energy are reported as transients. `Default=4`.

- min_distance:

  Minimum distance between two reported transients, in depth units
  (numeric, \>= 0); of candidates closer together only the strongest is
  kept. `NULL` uses `upperPeriod / 4`. `Default=NULL`.

- window:

  Length of the running median that forms the baseline of the energy, in
  depth units (numeric, \> 0). `NULL` uses `2 * upperPeriod`.
  `Default=NULL`.

## Value

A list with `profile` (data frame: depth, energy, ratio to the running
median) and `transients` (data frame: depth, ratio of the detected
peaks, strongest first).

## Details

The transform cannot see a jump of the sedimentation rate itself (the
cycles stay continuous); it sees the phase jumps that a gap in time
causes.

## References

He, D., Cao, H., Wang, S. and Chen, X. (2019). Time-reassigned
synchrosqueezing transform: the algorithm and its applications in
mechanical signal processing. Mech. Syst. Signal Process. 117, 255-279.

Yu, G. et al. (2023). Synchro-transient-extracting transform for the
analysis of signals with both harmonic and impulsive components.
arXiv:2306.01320.

## Examples

``` r
set.seed(2)
t <- seq(0, 2000, by = 0.5)
x_t <- sin(2 * pi * t / 23.7) + 0.6 * sin(2 * pi * t / 41) + 0.4 * cos(2 * pi * t / 405)
keep <- !(t > 800 & t < 900) & !(t > 1500 & t < 1550)        # two gaps in time (hiatuses)
z <- seq(0, by = 2, length.out = sum(keep))
x <- x_t[keep] + rnorm(length(z), sd = 0.5)
tr <- detect_transients(data.frame(z, x), lowerPeriod = 20, upperPeriod = 400)
tr$transients
#> [1] depth ratio
#> <0 rows> (or 0-length row.names)
```
