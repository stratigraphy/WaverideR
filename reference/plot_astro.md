# Six-panel summary of an inversion

Six-panel summary of an inversion

## Usage

``` r
plot_astro(
  fit,
  peaks = NULL,
  age_ma = NULL,
  file = NULL,
  width = 11,
  height = 7,
  res = 300
)
```

## Arguments

- fit:

  Inversion result to plot: the list returned by
  [`astro_invert`](https://stratigraphy.github.io/WaverideR/reference/astro_invert.md)
  (with elements `k`, `state`, `draws` and `windows`). Required, no
  default.

- peaks:

  Measured spectral peak periods in kyr (numeric vector). If given,
  panel (b) shows the
  [`k_scan`](https://stratigraphy.github.io/WaverideR/reference/k_scan.md)
  of these peaks (number of terms matched by each precession constant
  k); `NULL` leaves panel (b) empty. `Default=NULL`.

- age_ma:

  Age of the record in Ma (numeric, single value). If given, the
  AstroGeo22 (Farhat et al. 2022) value of k at that age is added to
  panel (c) and the fitted Earth-Moon distance and length of day are
  plotted at that age in panels (d) and (e); `NULL` skips these
  comparisons. `Default=NULL`.

- file:

  Path of an output file (character). A name ending in `".png"` writes a
  PNG file, any other name a PDF file; `NULL` draws on the current
  device. `Default=NULL`.

- width:

  Width of the output file in inches (numeric, \> 0); only used when
  `file` is given. `Default=11`.

- height:

  Height of the output file in inches (numeric, \> 0); only used when
  `file` is given. `Default=7`.

- res:

  Resolution of the PNG output in dpi (numeric, \> 0); only used when
  `file` ends in `".png"`. `Default=300`.
