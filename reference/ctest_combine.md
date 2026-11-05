# Combine several CTest results by multiplying their densities

Reproduces the combination step of mmc8
(reconstruct_parameters_Devonian.R), using a base-R weighted quantile
instead of Hmisc::wtd.quantile.

## Usage

``` r
ctest_combine(results, what = "k", n = 1000)
```

## Arguments

- results:

  CTest results to combine (list, ideally named): each element is a
  `ctest` object from
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)
  or a numeric vector of samples (e.g. precession constants in years).
  The kernel densities of all elements are multiplied. Required, no
  default.

- what:

  Column of the `samples` data frame of the `ctest` objects that is
  combined (character), e.g. `"k"` (precession constant, yr), `"SR"`
  (sedimentation rate, m/kyr), `"obl_ks3"` or `"g2"`; ignored for
  numeric vectors. `Default="k"`.

- n:

  Number of points of the grid on which the densities are evaluated and
  multiplied (integer, \> 1). `Default=1000`.

## Value

list: grid \`x\`, combined density, mean (median of the combined
density) and sd (half the inter-quartile span, as in the paper)
