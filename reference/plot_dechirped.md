# Plot a de-chirped scalogram

Plot a de-chirped scalogram

## Usage

``` r
plot_dechirped(
  x,
  style = c("outline", "significance", "plain", "tuned"),
  show = NULL,
  level = 0.95,
  palette = "viridis",
  main = NULL,
  xlab = "Depth",
  ylab = "Period",
  scale = "auto",
  outline_lwd = 0.6
)
```

## Arguments

- x:

  Result of
  [`analyse_dechirped`](https://stratigraphy.github.io/WaverideR/reference/analyse_dechirped.md)
  (list of class `wr_dechirped`). Required, no default.

- style:

  What to draw (character), one of:

  - `"outline"`: the de-chirped depth map with the cells that are
    significant against the warped surrogates outlined in white (needs
    `n_sim > 0` in
    [`analyse_dechirped`](https://stratigraphy.github.io/WaverideR/reference/analyse_dechirped.md));

  - `"significance"`: the ratio to the 95th percentile of the warped
    surrogates (log10, diverging colour scale; needs `n_sim > 0`);

  - `"plain"`: the de-chirped depth map;

  - `"tuned"`: the transform in the time of the track (periods in kyr).

  Without surrogates `"outline"` and `"significance"` fall back to
  `"plain"` with a message. `Default="outline"` (the first element of
  `c("outline", "significance", "plain", "tuned")`).

- show:

  Astronomical periods in kyr to draw as dotted lines (numeric vector,
  e.g. `c(405, 100, 41)`): in the depth styles at the period they have
  at the tracked rate (period x local rate), in the `"tuned"` style as
  horizontal lines. `NULL` draws none. `Default=NULL`.

- level:

  Significance level of the outline (numeric, between 0 and 1): cells
  with a p-value of at most `1 - level` are outlined. Only used when
  `style="outline"`. `Default=0.95`.

- palette:

  Name of a [`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html)
  palette for the colour scale, a character string, for example
  `"viridis"`, `"Plasma"`, `"Inferno"`, `"YlOrRd"`, `"Spectral"`,
  `"Blue-Red"`, `"Zissou 1"`; run
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  to list all 115 names. Low values get the first colour of the palette.
  Names of the viridis family (`"viridis"`, `"magma"`, `"inferno"`,
  `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`) give the
  colours of the viridis package, as in the wavelet plots.
  `Default="viridis"`.

- main:

  Title of the plot (character string); `NULL` gives an automatic title
  that names the transform and the style. `Default=NULL`.

- xlab, ylab:

  Labels of the depth and period axes (character strings) of the depth
  styles; the `"tuned"` style always uses "Time of the track (kyr)" and
  "Period (kyr)". Defaults: `xlab`: `Default="Depth"`; `ylab`:
  `Default="Period"`.

- scale:

  Colour scale as in all image plots of the package (character):
  `"auto"`, `"log"`, `"quantile"`, `"sqrt"` or `"linear"` (see
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)).
  The significance style always uses a fixed diverging scale centred on
  the 95th percentile of the surrogates. `Default="auto"`.

- outline_lwd:

  Line width of the significance outline (numeric, \> 0). `Default=0.6`.

## Value

Invisibly, the matrix shown.
