# Plot an astro_sedrate_test result

Main panel: the significance spectra (-log10 of the same-rate p) of
TimeOpt (blue), the tuned-band peaks (green) and the joint score
(black), each with its family-wise 5 % level as a dashed horizontal
line; below it the raw TimeOpt r2_opt and the noise-standardised
z-scores. Untestable rates are shaded.

## Usage

``` r
# S3 method for class 'astro_sedrate_test'
plot(x, ...)
```

## Arguments

- x:

  Result to plot: an object of class `astro_sedrate_test`, as returned
  by
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md).
  Required, no default.

- ...:

  Not used; present for compatibility with the generic
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

A plot.
