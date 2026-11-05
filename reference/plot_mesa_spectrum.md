# Plot a maximum entropy (MESA) spectrum

Plots the power spectrum of
[`mesa_spectrum`](https://stratigraphy.github.io/WaverideR/reference/mesa_spectrum.md)
on a logarithmic power axis, against period or frequency, with the
strongest peaks labelled.

## Usage

``` r
plot_mesa_spectrum(
  x,
  xaxis = c("period", "frequency"),
  peaks = TRUE,
  n_label = 5,
  xlim = NULL,
  log_x = FALSE,
  col = "black",
  main = NULL,
  ...
)
```

## Arguments

- x:

  Maximum entropy spectrum to plot: the list returned by
  [`mesa_spectrum`](https://stratigraphy.github.io/WaverideR/reference/mesa_spectrum.md).
  Required, no default.

- xaxis:

  Quantity on the x axis (character): `"period"` (in the depth/time unit
  of the input record) or `"frequency"` (cycles per unit).
  `Default="period"` (the first element of `c("period", "frequency")`).

- peaks:

  Mark the spectral peaks of `x$peaks` that fall inside `xlim` with red
  triangles and label the strongest ones with their period (logical,
  TRUE or FALSE). `Default=TRUE`.

- n_label:

  Number of strongest peaks (inside `xlim`) that are labelled with their
  period (integer, \>= 0). `Default=5`.

- xlim:

  Range of the x axis, in the unit of `xaxis` (numeric vector of length
  2); `NULL` uses the full range of positive frequencies (or the
  corresponding periods). `Default=NULL`.

- log_x:

  Draw the x axis on a logarithmic scale (logical, TRUE or FALSE); the
  power axis is always logarithmic. `Default=FALSE`.

- col:

  Colour of the spectrum line (a colour name, hexadecimal string or
  number). `Default="black"`.

- main:

  Title of the plot (character); `NULL` gives a title with the AR order
  and the order-selection criterion. `Default=NULL`.

- ...:

  Further graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) (e.g. `lwd`,
  `ylim`).

## Value

Called for its plot; invisibly returns `x`.

## Examples

``` r
set.seed(1)
t <- seq(0, 600, by = 2)
y <- sin(2 * pi * t / 100) + rnorm(length(t), sd = 0.5)
plot_mesa_spectrum(mesa_spectrum(data.frame(t, y)))
```
