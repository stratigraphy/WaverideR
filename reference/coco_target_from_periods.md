# Target periodogram from a list of orbital periods (period2spectrum.m)

Sum of sinusoids at the orbital periods, sampled from t1 to t2 (kyr),
and its zero-padded periodogram between f1 and f2 (cycles/kyr).

## Usage

``` r
coco_target_from_periods(
  orbit = NULL,
  age_ma = NULL,
  t1 = 3e+05,
  t2 = 308000,
  step = 1,
  f1 = 0,
  f2 = NULL,
  pad = 10000,
  weights = NULL
)
```

## Arguments

- orbit:

  Astronomical periods of the target sinusoids, in kyr (numeric vector,
  \> 0; sorted longest first internally), e.g.
  `c(405, 125, 95, 41, 23, 22, 19)`. `NULL` uses the nine periods at
  `age_ma`. `Default=NULL`.

- age_ma:

  Age in Ma (single numeric value, 0 to 4429) at which the nine periods
  are taken from
  [`orbit9_at_age`](https://stratigraphy.github.io/WaverideR/reference/orbit9_at_age.md)
  (La2004 g and s frequencies with the AstroGeo22 precession constant);
  only used when `orbit` is `NULL`, and one of the two must be given.
  `Default=NULL`.

- t1, t2:

  Start and end of the time axis on which the synthetic sum of sinusoids
  is sampled, in kyr (numeric, `t2 > t1`); the length `t2 - t1` sets the
  width of the target peaks. `Default=300000` and `Default=308000`.

- step:

  Time step of the synthetic series, in kyr (numeric, \> 0).
  `Default=1`.

- f1:

  Lowest frequency of the target that is returned, in cycles/kyr
  (numeric, \>= 0). `Default=0`.

- f2:

  Highest frequency of the target that is returned, in cycles/kyr
  (numeric, \> `f1`); `NULL` uses 1.2 times the highest orbital
  frequency, as in the Acycle GUI. `Default=NULL`.

- pad:

  Length (number of points) to which the synthetic series is zero-padded
  for the periodogram (integer; at least the series length is used); a
  larger value gives a finer frequency grid. `Default=10000`.

- weights:

  Amplitudes of the sinusoids, one per period in `orbit` after sorting
  longest first (numeric vector). `NULL` uses 1.0 for eccentricity, 0.8
  for obliquity and 0.6 for precession (for nine periods via
  [`coco_fixed_target_weights`](https://stratigraphy.github.io/WaverideR/reference/coco_fixed_target_weights.md);
  otherwise the three longest get 1.0, the fourth 0.8 and the rest 0.6,
  as in Acycle period2spectrum.m). `Default=NULL`.

## Value

two-column matrix: frequency (cycles/kyr), power

## Details

The returned target carries attribute "orbit" so coco_classic() uses the
same periods to count contributing parameters.
