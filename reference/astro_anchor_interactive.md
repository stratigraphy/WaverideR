# Interactively tie a proxy record to an astronomical solution (random order, tie and untie, inspect)

A more flexible alternative to
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md).
The proxy record (top) and the astronomical solution (bottom) are shown
in one window with their peaks (or troughs). Tie points can be set,
removed and changed in any order: click a peak of the proxy record and
then a peak of the astronomical solution to tie them (or first the
astronomical peak and then the proxy peak); click a tied peak again to
untie it. The windows can be zoomed and moved with the buttons above the
plot, so that long records can be tied piece by piece. Press *Tune and
inspect* at any moment: a second window shows the record tuned to the
current tie points together with a cycle extracted from it, the
astronomical solution and the sedimentation rate between the tie points,
so that the result can be judged while tie points are still being added
or removed.

## Usage

``` r
astro_anchor_interactive(
  astro_solution = NULL,
  proxy_signal = NULL,
  proxy_min_or_max = "max",
  astrosolution_min_or_max = "max",
  pts = 3,
  ties = NULL,
  time_dir = TRUE,
  peaks_in_view = 12,
  inspect_period = NULL,
  inspect_band = 0.25,
  sedrate_factor = 1,
  inspect_wavelet = FALSE,
  locator_fun = NULL
)
```

## Arguments

- astro_solution:

  Astronomical solution to tie to (matrix or data frame with two
  columns: time/age (e.g. kyr) in column 1, value in column 2; rows with
  `NA` are removed and rows are sorted on column 1). Required, no
  default (`Default=NULL` gives an error).

- proxy_signal:

  Proxy record to be tied (matrix or data frame with two columns:
  depth/time in column 1, proxy value in column 2; rows with `NA` are
  removed). Use
  [`extract_signal`](https://stratigraphy.github.io/WaverideR/reference/extract_signal.md)
  first if a single cycle of the record should be tied. Required, no
  default (`Default=NULL` gives an error).

- proxy_min_or_max:

  Tie the maxima or the minima of the proxy record (character): `"max"`
  or `"min"`. `Default="max"`.

- astrosolution_min_or_max:

  Tie the maxima or the minima of the astronomical solution (character):
  `"max"` or `"min"`. `Default="max"`.

- pts:

  Number of points to the left and right of a value that the peak
  detection
  ([`max_detect`](https://stratigraphy.github.io/WaverideR/reference/max_detect.md),
  [`min_detect`](https://stratigraphy.github.io/WaverideR/reference/min_detect.md))
  compares it with (integer, \>= 1); larger values give fewer, more
  certain peaks. `Default=3`.

- ties:

  Existing tie points to start from (data frame with the columns
  `data_x` (proxy depth/time) and `insolation_x` (astronomical age),
  otherwise the first two columns are used; e.g. the output of
  [`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)).
  They are snapped to the nearest peaks. `NULL` starts without tie
  points. `Default=NULL`.

- time_dir:

  Direction of time in the proxy record (logical): `TRUE` if time
  increases with increasing depth/time values, `FALSE` if it decreases
  (the proxy axis is then reversed). `Default=TRUE`.

- peaks_in_view:

  Approximate number of peaks shown in each window at the start
  (numeric, \> 0). `Default=12`.

- inspect_period:

  Period of the cycle extracted from the tuned record in the inspection
  window, in the time unit of the astronomical solution (numeric, \> 0);
  `NULL` uses the median spacing of the peaks of the astronomical
  solution. `Default=NULL`.

- inspect_band:

  Relative half-width of the extracted band (numeric, between 0 and 1):
  `0.25` extracts 0.75 to 1.25 times `inspect_period`. `Default=0.25`.

- sedrate_factor:

  Factor applied to the sedimentation rate (depth units per time unit)
  shown in the inspection window (numeric, \> 0), for example 100 to go
  from m/kyr to cm/kyr. `Default=1`.

- inspect_wavelet:

  Also show the wavelet scalogram of the tuned record in the inspection
  window at the start (logical, TRUE or FALSE); it can be switched on
  and off with the *Wavelet* button. `Default=FALSE`.

- locator_fun:

  Function without arguments that returns the next mouse click in
  normalised device coordinates as a list with `x` and `y`, like
  [`locator`](https://rdrr.io/r/graphics/locator.html)`(1)`; returning
  `NULL` finishes. `NULL` uses `graphics::locator(1)`. Another function
  (for example one that replays a list of clicks) allows scripted use
  and testing. When no interactive device is available and `ties` is
  given, the snapped tie points are returned without opening a window.
  `Default=NULL`.

## Value

A data frame with the tie points sorted by depth/time: `data_x`,
`data_y` (proxy), `insolation_x`, `insolation_y` (astronomical
solution).

## Details

Buttons (click them with the left mouse button): *Tune and inspect*
(open or update the inspection window), *Undo* (remove the last change),
*Clear* (remove all tie points), *Zoom +*/*Zoom -* (both windows), *P
\<*, *P \>* (move the proxy window), *A \<*, *A \>* (move the
astronomical window), *Fit* (show the whole records), *Sync* (move the
astronomical window to the partners of the tie points in the proxy
window, or the other way round if there are none), *Wavelet* (switch the
wavelet scalogram of the tuned record in the inspection window on or
off) and *Done*. A right click also finishes. Tie points that are in the
wrong order (the tuned time would run backwards) are drawn in orange,
and the tuning refuses to run until they are fixed.

The result has the same four columns as the output of
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)
and can be used directly by
[`plot_astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md)
and
[`anchor2time`](https://stratigraphy.github.io/WaverideR/reference/anchor2time.md).

## See also

[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md),
[`plot_astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md),
[`anchor2time`](https://stratigraphy.github.io/WaverideR/reference/anchor2time.md)

## Examples

``` r
# \donttest{
# Tie points can be given as a starting point and are snapped to the nearest peaks;
# without an interactive graphics device (as here) they are returned as they are.
astrosignal <- na.omit(astrosignal_example)
tie_start <- data.frame(data_x = c(155, 160, 165), insolation_x = c(5040, 5400, 5560))
proxy <- grey[grey[, 1] < 175, ]
ties <- astro_anchor_interactive(astro_solution = astrosignal,
                                 proxy_signal = proxy,
                                 ties = tie_start)
ties
#>   data_x   data_y insolation_x insolation_y
#> 1 155.04 362.4783         5036    1.0902778
#> 2 160.01 380.1990         5393    0.6553615
#> 3 164.99 426.3102         5550    2.1563294
# in an interactive session, leave out 'ties' to tie the peaks with the mouse
# }
```
