# Convert a proxy record to the time domain using anchor points

Convert a proxy record to the time domain using anchor points made using
the
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)
function.

## Usage

``` r
anchor2time(
  anchor_points = NULL,
  data = NULL,
  genplot = FALSE,
  keep_editable = FALSE
)
```

## Arguments

- anchor_points:

  Anchor points: the output of
  [`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md),
  or a matrix or data frame with depth in column 1 and time in column 2
  (at least two rows, sorted by increasing depth). Time is linearly
  interpolated between the anchor points; above the first and below the
  last anchor point it is extrapolated with the sedimentation rate of
  the first/last anchored interval. Required (`Default=NULL` gives an
  error).

- data:

  Proxy record to convert from the depth to the time domain: a
  two-column matrix or data frame with depth in column 1 (same unit as
  column 1 of `anchor_points`, sorted by increasing depth) and the proxy
  value in column 2. Required (`Default=NULL` gives an error).

- genplot:

  Plot the result (logical, TRUE or FALSE). If `TRUE`, 3 plots are
  stacked on top of each other: 1) the original record in the depth
  domain, 2) the depth-time plot of the (extended) anchor points, 3) the
  record in the time domain. `Default=FALSE`.

- keep_editable:

  Keep the plot layout and graphical parameters after plotting so that
  extra features can be added to the plot (logical, TRUE or FALSE); if
  `FALSE` the previous [`par()`](https://rdrr.io/r/graphics/par.html)
  settings are restored. Only used when `genplot=TRUE`. `Default=FALSE`.

## Value

The output is a matrix with 2 columns. The first column is time. The
second column sedimentation proxy value.

If `genplot=TRUE` then 3 plots stacked on top of each other will be
plotted. Plot 1: the original data set. Plot 2: the depth time plot.
Plot 3: the data set in the time domain.

## Examples

``` r
# Use the age_model_zeeden example anchor points of Zeeden et al., (2013)
#to anchor the grey data set of Zeeden et al., (2013) in the time domain.

grey_time <- anchor2time(anchor_points=age_model_zeeden,
data=grey,
genplot=FALSE,
keep_editable=FALSE)
```
