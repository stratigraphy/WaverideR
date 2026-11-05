# Plot a continuous track on its scalogram

Plots the scalogram with the track(s) of
[`track_ridge_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_ridge_hmm.md)
or
[`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
and their credible bands, the confidence of the track and, when
available, the sedimentation rate and the age model with its
uncertainty.

## Usage

``` r
plot_track(
  track,
  wavelet,
  show = NULL,
  flags = NULL,
  palette = "viridis",
  keep_editable = FALSE,
  scale = "auto"
)
```

## Arguments

- track:

  The track to draw (data frame of class `wr_track`): the output of
  [`track_ridge_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_ridge_hmm.md),
  [`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
  or
  [`track_chirp_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_chirp_hmm.md)
  (for
  [`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md)
  pass its element `track`). Required, no default.

- wavelet:

  The scalogram the track was made from, or another one on the same
  depth axis to display (list of class `analyze.wavelet`, e.g. the
  element `scalogram` of
  [`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md)).
  Required, no default.

- show:

  Periods (in time units, numeric vector, e.g. `c(405, 41, 23.7)`) of
  the astronomical cycles whose depth periods are drawn for a
  sedimentation-rate track; each is matched within 10 percent to the
  periods stored in the track (attribute `astro_periods`). For a track
  without stored periods (e.g. from
  [`track_chirp_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_chirp_hmm.md))
  the given periods are drawn as they are. `NULL` draws all stored
  periods (or the tracked period of a
  [`track_ridge_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_ridge_hmm.md)
  result). `Default=NULL`.

- flags:

  Optional intervals to shade in all panels: a data frame with the
  columns `from` and `to` (depths), e.g. the element `flags` of
  [`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md),
  or `TRUE` to use the attribute `flags` of `track`; anything else draws
  no shading. `Default=NULL`.

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

- keep_editable:

  Keep the graphical parameters (`par`) and the panel layout after
  plotting, so that more elements can be added to the last panel
  (logical, TRUE or FALSE); `FALSE` restores the previous settings.
  `Default=FALSE`.

- scale:

  Colour scale of the scalogram, as in all image plots of the package
  (see
  [`plot_wavelet`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)),
  one of `"auto"`, `"log"`, `"quantile"`, `"sqrt"`, `"linear"`. `"auto"`
  uses a log scale for synchrosqueezed (mostly empty) maps and quantiles
  otherwise; the power of a non-synchrosqueezed scalogram is first
  divided by its mean per period. `Default="auto"`.

## Value

Invisibly the track.
