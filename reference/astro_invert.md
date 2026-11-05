# Joint inversion of several measured terms, with Monte Carlo uncertainty

Joint inversion of several measured terms, with Monte Carlo uncertainty

## Usage

``` r
astro_invert(
  windows,
  nsim = 5000,
  obl_deg = AC$obl0_deg,
  k_range = c(10, 34),
  seed = 1,
  gs = LA2004
)
```

## Arguments

- windows:

  Measured period ranges of the k-bearing (obliquity/precession) terms,
  in kyr (named list; each element a numeric vector `c(low, high)`).
  Allowed names: `"k+s3"` (obliquity), `"k+g5"`, `"k+g2"`, `"k+g4"` and
  `"k+g3"` (precession), e.g. `list("k+s3" = c(29.0, 30.4))`. In each
  draw one period is sampled uniformly in every window, converted into a
  precession frequency, and the frequencies of all windows are averaged.
  Required, no default.

- nsim:

  Number of Monte Carlo draws (integer, \> 0); draws whose k falls
  outside `k_range` are discarded afterwards. `Default=5000`.

- obl_deg:

  Obliquity in degrees used in the precession relation when k is
  converted into the Earth-Moon state (numeric); see
  [`astro_obliquity`](https://stratigraphy.github.io/WaverideR/reference/astro_obliquity.md)
  to obtain it for a given age. `Default=AC$obl0_deg`, i.e. today's
  obliquity, 23.2545 degrees.

- k_range:

  Accepted range of the precession constant k in kyr (numeric vector
  `c(min, max)`); draws outside it are rejected. `Default=c(10, 34)`.

- seed:

  Seed of the random number generator (integer), passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) so that the Monte
  Carlo results can be reproduced. `Default=1`.

- gs:

  Fundamental secular frequencies in arcsec/yr (list with numeric
  vectors `g` and `s` of length 9, as in Laskar et al. 2004); g2-g5 and
  s3 are perturbed with the Laskar et al. (2011) uncertainties in each
  draw. `Default=LA2004`, the La2004 values stored inside the function.

## Value

list(k = summary of k, state = summary of the Earth-Moon state, draws =
the raw draws)
