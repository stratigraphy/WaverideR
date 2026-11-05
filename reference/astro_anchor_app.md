# Interactive Shiny app to tie a proxy record to an astronomical solution (tune live while tying)

A browser-based version of
[`astro_anchor_interactive`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor_interactive.md).
One window shows three linked rows: the proxy record, the astronomical
solution and, underneath, the record tuned to the current tie points
together with a cycle extracted from it (and optionally the astronomical
solution and the sedimentation rate). The tuned row is updated while tie
points are made, moved or removed, so the effect of every change on the
fit can be judged on the go. Tie points are numbered and colour coded in
all rows and joined by lines between the two records.

## Usage

``` r
astro_anchor_app(
  astro_solution = NULL,
  proxy_signal = NULL,
  proxy_min_or_max = "max",
  astrosolution_min_or_max = "max",
  pts = 3,
  ties = NULL,
  time_dir = TRUE,
  peaks_in_view = 12,
  tie_on = c("raw", "cycle"),
  proxy_flow = NULL,
  proxy_fhigh = NULL,
  proxy_roll = 3,
  inspect_flow = NULL,
  inspect_fhigh = NULL,
  inspect_roll = 3,
  sedrate_factor = 1,
  launch_browser = TRUE
)
```

## Arguments

- astro_solution:

  Astronomical solution to tie to (matrix or data frame with two
  columns: time/age (e.g. kyr) in column 1, value in column 2), or a
  named list of such solutions to choose from in the app (the first one
  is shown at the start). In the app the target can also be switched to
  a metronome (a cosine with a period, phase and time interval you set),
  or to a Laskar solution (eccentricity, obliquity, climatic precession
  or ETP of La2004, La2010a-d, La2011 eccentricity, or 65 N insolation)
  that is downloaded on request with
  [`getLaskar`](https://rdrr.io/pkg/astrochron/man/getLaskar.html) when
  you press the download button. Required, no default (`Default=NULL`
  gives an error).

- proxy_signal:

  Proxy record to be tied (matrix or data frame with two columns:
  depth/time in column 1, proxy value in column 2; rows with `NA` are
  removed). Required, no default (`Default=NULL` gives an error).

- proxy_min_or_max:

  Peaks offered when snapping on the proxy record (character): `"max"`
  (maxima) or `"min"` (minima). `Default="max"`.

- astrosolution_min_or_max:

  Peaks on which the tie points of the astronomical solution are placed
  (character): `"max"` (maxima) or `"min"` (minima). `Default="max"`.

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
  The astronomical positions are moved to the nearest maximum/minimum,
  the proxy positions are used as given. `NULL` starts without tie
  points. `Default=NULL`.

- time_dir:

  Direction of time in the proxy record (logical): `TRUE` if time
  increases with increasing depth/time values, `FALSE` if it decreases.
  `Default=TRUE`.

- peaks_in_view:

  Approximate number of peaks shown in each window at the start
  (numeric, \> 0). `Default=12`.

- tie_on:

  Proxy record on which the tie points are placed at the start
  (character): `"raw"` (the raw record) or `"cycle"` (a cycle extracted
  from it with the Taner filter set by `proxy_flow`, `proxy_fhigh` and
  `proxy_roll`). Can be changed in the app. `Default=c("raw", "cycle")`,
  i.e. `"raw"`.

- proxy_flow, proxy_fhigh:

  Low and high cut-off frequencies of the Taner filter that extracts the
  proxy cycle, in cycles per depth/time unit of the proxy record
  (numeric, 0 \< `proxy_flow` \< `proxy_fhigh`; give both). `NULL` uses
  the strongest period of the wavelet spectrum of the proxy -/+ 25
  percent. They can be changed in the app. `Default=NULL`.

- proxy_roll:

  Roll-off of the proxy Taner filter as the exponent `x` of
  `roll = 10^x` (numeric; see
  [`taner`](https://rdrr.io/pkg/astrochron/man/taner.html)).
  `Default=3`.

- inspect_flow, inspect_fhigh:

  Low and high cut-off frequencies of the Taner filter that extracts the
  cycle from the tuned record, in cycles per time unit of the
  astronomical solution (numeric, 0 \< `inspect_flow` \<
  `inspect_fhigh`; give both). `NULL` uses the median spacing of the
  peaks of the astronomical solution -/+ 25 percent. `Default=NULL`.

- inspect_roll:

  Roll-off of the Taner filter of the tuned record as the exponent `x`
  of `roll = 10^x` (numeric). `Default=3`.

- sedrate_factor:

  Factor applied to the sedimentation rate (depth units per time unit)
  that is shown (numeric, \> 0), for example 100 to go from m/kyr to
  cm/kyr. `Default=1`.

- launch_browser:

  Open the app in the default web browser (logical, TRUE or FALSE);
  `FALSE` prints the local address (or uses the RStudio viewer when R is
  run there). Passed to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html) as
  `launch.browser`. `Default=TRUE`.

## Value

A data frame with the tie points sorted by depth/time: `data_x`,
`data_y` (proxy), `insolation_x`, `insolation_y` (astronomical
solution). Without an interactive session and with `ties` the starting
tie points are returned without starting the app (those on the
astronomical solution moved to the nearest peak).

## Details

**What is tied.** Tie points on the astronomical solution always sit on
a maximum or minimum (see `astrosolution_min_or_max`). On the proxy side
they can be placed on the *raw* record or on a cycle that is extracted
from it with a Taner filter
([`taner`](https://rdrr.io/pkg/astrochron/man/taner.html); choose *Cycle
extracted from the proxy* and give `flow`, `fhigh` and the exponent of
the roll-off, `roll = 10^x`). The tuning is always applied to the raw
record. Proxy tie points do not have to lie on a maximum or minimum:
with *Snap to peaks* switched off they are placed exactly where you
click, and every proxy tie point can be moved left and right with the
buttons (fine and coarse steps), by typing its position, by typing its
offset from the nearest maximum/minimum, or by dragging. The cycle shown
in the tuned record is extracted with the same kind of Taner filter.

**Tie** a point: click one record (red) and then the other (in any
order). **Select** a tie point by clicking it. **Move** it by clicking a
free position of one of the records afterwards, by dragging over the
record from the tie point to the new position, by the arrow buttons, or
by editing the numbers. **Untie** it by double clicking it or with the
*Untie selected* button. Tie points in the wrong order (the tuned time
would run backwards) are marked with a black cross; the tuned row waits
until they are fixed.

**Target and windows.** The card *Astronomical solution* selects the
target (a provided solution, a metronome, or a Laskar solution), the
time interval it spans and whether maxima or minima are tied. Two
sliders above the plot show the visible window (from, to) of the proxy
record and of the astronomical solution; type new numbers to set them.
With the mouse, the *wheel* zooms the record under the pointer around
the pointer and *Shift + drag* (or dragging with the right or middle
mouse button) moves it from left to right. *Sync* brings the
astronomical window to the partners of the tie points in the proxy
window, *Show both records in full* resets both.

The sidebar also has buttons to undo changes and to finish (*Done*).
Closing the browser tab also finishes and returns the tie points.

The app runs on the local computer only (a Shiny app started with
[`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html)); nothing is sent
to other servers. The package shiny has to be installed.

The result has the same four columns as the output of
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)
and can be used by
[`plot_astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md)
and
[`anchor2time`](https://stratigraphy.github.io/WaverideR/reference/anchor2time.md).
`data_y` and `insolation_y` are the values of the records at the tie
positions.

## See also

[`astro_anchor_interactive`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor_interactive.md),
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md),
[`plot_astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md),
[`anchor2time`](https://stratigraphy.github.io/WaverideR/reference/anchor2time.md)

## Examples

``` r
# \donttest{
# Without an interactive session (as here) the starting tie points are returned as they are.
astrosignal <- na.omit(astrosignal_example)
tie_start <- data.frame(data_x = c(155, 160, 165), insolation_x = c(5040, 5400, 5560))
proxy <- grey[grey[, 1] < 175, ]
astro_anchor_app(astro_solution = astrosignal, proxy_signal = proxy, ties = tie_start)
#>   data_x   data_y insolation_x insolation_y
#> 1    155 373.5227         5036    1.0902778
#> 2    160 373.6639         5393    0.6553615
#> 3    165 392.9241         5550    2.1563294
# in an interactive session, leave out 'ties' to start the app and tie the records with the mouse
# }
```
