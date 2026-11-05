# Plot DYNOT results

The `plot_dynot` function plots the input data and the DYNOT noise model
(median and confidence envelopes, reversed y-axis) as in the Acycle
DYNOT window. Translated from the plotting part of Acycle v3.1 DYNOS.m
(Mingsong Li).

## Usage

``` r
plot_dynot(dynot_res = NULL, plot_data = TRUE, x_lab = "time (kyr)")
```

## Arguments

- dynot_res:

  DYNOT result to plot: the list returned by
  [`dynot`](https://stratigraphy.github.io/WaverideR/reference/dynot.md)
  (its `data`, `percentile`, `interval_count` and `median_index`
  elements are used). Required (`Default=NULL` gives an error).

- plot_data:

  Also plot the input data in an upper panel above the noise model
  (logical, TRUE or FALSE); with `FALSE` only the noise panel (median
  and up to five confidence envelopes, reversed y-axis) is drawn.
  `Default=TRUE`.

- x_lab:

  Label of the time axis of both panels (character).
  `Default="time (kyr)"`.

## Value

No return value; a plot is drawn.

## References

Li, M., Hinnov, L.A., Huang, C., Ogg, J.G., 2018. Sedimentary noise and
sea levels linked to land-ocean water exchange and obliquity forcing.
Nature Communications 9, 1004.
[doi:10.1038/s41467-018-03454-y](https://doi.org/10.1038/s41467-018-03454-y)

## Author

Michiel Arts (R translation); original MATLAB code: Mingsong Li
(Acycle).

## Examples

``` r
# \donttest{
set.seed(2)
t <- cumsum(runif(900, 1, 3))
x <- sin(2 * pi * t / 405) + sin(2 * pi * t / 41) + rnorm(900)
res <- dynot(data = cbind(t, x), n_sim = 20, n_out = 200, seed = 1)
plot_dynot(res)

# }
```
