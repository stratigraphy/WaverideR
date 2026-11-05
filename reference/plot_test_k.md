# Comparison figure of the precession-constant test

Plots the result of
[`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md):
the MTM spectrum with the AstroGeo22 periods, search windows and the
peaks used (A), k, Earth-Moon distance and length of day against
AstroGeo22 and the Devonian composite (B-D), AstroGeo22 through time
with the estimate and optional published constraints (E-G), and the
pre-scan of the record.

## Usage

``` r
plot_test_k(res, constraints = NULL, age_ma = res$spec$age_ma)
```

## Arguments

- res:

  Result to plot (required, no default): the list returned by
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md).

- constraints:

  Path to a CSV file of published constraints drawn in panels E-G
  (character, or `NULL` for none); see `constraints` in
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md)
  for the columns. `Default=NULL`.

- age_ma:

  Age of the record in Ma (numeric), where the estimate is drawn in
  panels E-G and at which AstroGeo22 is evaluated.
  `Default=res$spec$age_ma`, the age used in
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md).

## Value

Called for its plot.
