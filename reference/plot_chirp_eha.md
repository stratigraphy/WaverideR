# Plot a chirp-matched EHA map, with or without the surrogate significance

Two styles of "enhanced spectrum" for a
[`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md)
result: `"significance"` shows the line-to-continuum ratio divided by
the 95th percentile of the warped surrogates (`Power_sig`, from
[`chirp_eha_significance`](https://stratigraphy.github.io/WaverideR/reference/chirp_eha_significance.md);
above 1, i.e. above 0 on the log scale, is significant), `"outline"`
shows the de-chirped map itself (`Power_norm`) with the significant
cells outlined, and `"plain"` the de-chirped map alone.

## Usage

``` r
plot_chirp_eha(
  ce,
  style = c("outline", "significance", "plain"),
  track = NULL,
  show = NULL,
  level = 0.95,
  palette = "viridis",
  main = NULL,
  xlab = "Depth",
  ylab = "Period",
  scale = "auto"
)
```

## Arguments

- ce:

  Chirp-matched EHA result to plot: the list returned by
  [`analyse_chirp_eha`](https://stratigraphy.github.io/WaverideR/reference/analyse_chirp_eha.md);
  for `style = "outline"` or `"significance"` it must first be passed
  through
  [`chirp_eha_significance`](https://stratigraphy.github.io/WaverideR/reference/chirp_eha_significance.md)
  (which adds `Power_sig` and `Power_pval`). Required, no default.

- style:

  Plot style (character), one of:

  - `"outline"`: the de-chirped map (`Power_norm`) with the cells
    significant at `level` outlined in white;

  - `"significance"`: log10 of the ratio of the observed map to the 95th
    percentile of the warped surrogates (`Power_sig`), with a diverging
    blue-red scale centred on 0 (above 0 is significant);

  - `"plain"`: the de-chirped map alone (no significance needed).

  `Default="outline"` (the first element of
  `c("outline", "significance", "plain")`).

- track:

  Sedimentation-rate track used to draw the expected periods of the
  astronomical cycles in `show`: a `wr_track` object or a
  [`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md)
  result (its `track` element is used), i.e. a list or data frame with
  `depth` and `sedrate` (depth units per kyr). Nothing is drawn when
  `track` or `show` is `NULL`. `Default=NULL`.

- show:

  Astronomical periods in kyr (numeric vector, e.g.
  `c(405, 100, 40, 20)`) drawn as dotted lines at period = `show` times
  `track$sedrate` along the record; only used with `track`.
  `Default=NULL`.

- level:

  Confidence level of the significance outline of `style = "outline"`
  (numeric, between 0 and 1): cells with surrogate p-value
  `Power_pval <= 1 - level` are outlined. `Default=0.95`.

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

  Title of the plot (character); `NULL` gives a title that describes the
  chosen `style`. `Default=NULL`.

- xlab:

  Label of the depth axis (character). `Default="Depth"`.

- ylab:

  Label of the period axis (character). `Default="Period"`.

- scale:

  Spacing of the colour breaks of the colour scale, as in all image
  plots of the package (character): `"quantile"` (histogram-equalised),
  `"log"`, `"sqrt"`, `"linear"` or `"auto"` (logarithmic for mostly
  empty maps or values spanning more than three orders of magnitude,
  quantiles otherwise); see
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md).
  Not used for `style = "significance"`, which has a fixed diverging
  scale. `Default="auto"`.

## Value

Invisibly, the matrix shown.
