# Precession-constant test for a tuned record

Tests a tuned (time-domain) record with the precession constant k.
Significant multi-taper (MTM) spectral peaks of the record, optionally
pooled with other proxies on the same time scale, are assigned to the
obliquity (k+s3, k+s6) and precession (k+g5, k+g2, k+g4) terms. Every
combination of candidate peaks is tried over a wide range of k (no
window is centred on a reference value), and the combinations in which
all terms give one consistent k are kept. Because the terms share one k,
a correct assignment gives the same k from every term. The result is
compared with AstroGeo22 (Farhat et al. 2022) and a Devonian range, and
combined with two further checks of the tuned record: the spectral peaks
of the eccentricity, obliquity and precession terms, and the amplitude
modulation (405 kyr from the 100 kyr envelope, 100 kyr from the
precession envelope) against phase-randomised surrogates. This is the
updated version of the test (version 3), which works better for
Silurian/Devonian data than a search around the AstroGeo22 value.

## Usage

``` r
test_tuned(
  URAN_time,
  others = NULL,
  tbw = 2,
  pad = 5,
  siglevel = 0.9,
  ftest = FALSE,
  tol = 0.04,
  scan = TRUE,
  scan_band = c(15, 23),
  scan_win = 200,
  scan_thr = 0.5,
  scan_min = 100,
  n_surr = 200,
  seed = 42,
  plot = TRUE,
  verbose = FALSE,
  age_ma = 386,
  k_ref = NULL,
  k_ref_range = NULL,
  k_dev = c(20.4, 21.6),
  k_search = c(19.5, 23.5),
  constraints = NULL
)
```

## Arguments

- URAN_time:

  Tuned record to test (required, no default): data frame or matrix with
  time in kyr in column 1 and the proxy value in column 2, e.g. the
  output of a tuning function. It is sorted, linearly interpolated to 1
  kyr, detrended and standardised before the analysis.

- others:

  Other tuned proxies on the same time scale (list of data
  frames/matrices with the same two-column layout as `URAN_time`, or
  `NULL`); they are cut to the time range of `URAN_time` and their peaks
  are pooled with those of `URAN_time` (a term found in several proxies
  adds support). `Default=NULL`.

- tbw:

  Time-bandwidth product of the MTM spectrum (numeric, \> 0; `tbw` of
  [`astrochron::mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)). It
  sets the resolution: the precession terms k+g5 and k+g2 only separate
  when `tbw <= (record length in kyr) / 1600` (about 3 for 5.5 Myr). A
  larger `tbw` merges them into one broad peak and k cannot be measured
  (a warning is printed). `Default=2`.

- pad:

  Zero-padding factor of the MTM spectrum (numeric, \>= 1; `padfac` of
  [`astrochron::mtm`](https://rdrr.io/pkg/astrochron/man/mtm.html)).
  `Default=5`.

- siglevel:

  Confidence level a spectral peak must reach to be used (numeric,
  between 0 and 1): the AR1 confidence level at the peak (and, with
  `ftest = TRUE`, the harmonic F-test level). Only peaks with periods
  between 12 and 40 kyr are considered. `Default=0.9`.

- ftest:

  Also require the harmonic F-test confidence level to reach `siglevel`
  (logical). The default `FALSE` uses the AR1 level only, because the
  F-test tests for a pure sinusoid and rejects broad, strong peaks of a
  tuned record. `Default=FALSE`.

- tol:

  Half-width of the search window around each predicted term period, as
  a fraction of the period (numeric, e.g. `0.04` = +/- 4%).
  `Default=0.04`.

- scan:

  Remove intervals with little precession imprint (set them to zero,
  with a 10-kyr cosine taper at the edges) before the spectrum is built
  (logical); the pre-scan is controlled by `scan_band`, `scan_win`,
  `scan_thr` and `scan_min`. `Default=TRUE`.

- scan_band:

  Precession band of the pre-scan, in kyr (numeric vector
  `c(shortest, longest)` period); the record is band-passed
  ([`astrochron::taner`](https://rdrr.io/pkg/astrochron/man/taner.html))
  in this band and its Hilbert envelope gives the band variance.
  `Default=c(15, 23)`.

- scan_win:

  Length of the running window over which the precession-band variance
  is averaged, in kyr (numeric, \> 0; at most a third of the record is
  used). `Default=200`.

- scan_thr:

  An interval is weak when its averaged precession-band variance is
  below `scan_thr` times its median (numeric, \>= 0). `Default=0.5`.

- scan_min:

  Shortest piece kept or removed by the pre-scan, in kyr (numeric, \>=
  0): weak pieces shorter than this are kept and kept pieces shorter
  than this are removed. `Default=100`.

- n_surr:

  Number of phase-randomised surrogates for the amplitude-modulation
  test (integer, \> 0). `Default=200`.

- seed:

  Seed of the random number generator (integer), so that the surrogates
  can be reproduced. `Default=42`.

- plot:

  Draw the comparison figure with
  [`plot_test_k`](https://stratigraphy.github.io/WaverideR/reference/plot_test_k.md)
  (logical). `Default=TRUE`.

- verbose:

  Also print the peak tables and all peak assignments of the
  precession-constant test (logical, TRUE or FALSE); the compact report
  is always printed. `Default=FALSE`.

- age_ma:

  Age of the record in Ma (numeric); the AstroGeo22 reference k and its
  envelope are taken at this age
  ([`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md))
  unless `k_ref` and `k_ref_range` are given. `Default=386`.

- k_ref:

  Reference precession constant, as a period in kyr (numeric); used for
  the resolution check of `tbw`, the old search windows
  (`k_search = NULL`) and the plot. `NULL` uses AstroGeo22 at `age_ma`.
  `Default=NULL`.

- k_ref_range:

  Range of the reference precession constant in kyr (numeric vector
  `c(lower, upper)`); the test passes when the 1-sigma interval of the
  estimate overlaps it (or the 20.70 kyr composite of Zeeden et al.,
  2023). `NULL` uses the AstroGeo22 envelope at `age_ma`.
  `Default=NULL`.

- k_dev:

  Plausible range of k for the test verdict, in kyr (numeric vector
  `c(lower, upper)`); the default is the Devonian range.
  `Default=c(20.4, 21.6)`.

- k_search:

  Range of k in which peaks are searched, in kyr (numeric vector
  `c(lower, upper)`), independent of AstroGeo22. The search window of
  each term runs from its period at `k_search[1]` to its period at
  `k_search[2]` (plus or minus `tol`). `NULL` gives the old windows
  around the AstroGeo22 period, which can only return k within about 1
  kyr of AstroGeo22 and so cannot test it. `Default=c(19.5, 23.5)`.

- constraints:

  Path to a CSV file of published constraints (character, or `NULL` for
  none), with columns `set`, `observable`, `age_ma`, `a_RE`, `a_sd`,
  `lod_h`, `lod_sd`, `k_arcsec` and `k_sd`, drawn in panels E-G of the
  figure. `Default=NULL`.

## Value

Invisibly a list with the combined k (`k`, `k_obliquity`,
`k_precession`, `k_ranges`), the peak assignment tables (`k_terms`,
`solutions`), the verdicts and score of the three tests, and the spectra
used for plotting.

## References

Zeeden et al. (2023), see
[`zeeden2023_table1`](https://stratigraphy.github.io/WaverideR/reference/zeeden2023_table1.md).
Farhat, M., Auclair-Desrotour, P., Boue, G., Laskar, J. (2022). The
resonant tidal evolution of the Earth-Moon distance. Astronomy &
Astrophysics 665, L1.

## Examples

``` r
# \donttest{
# synthetic tuned record with a known precession constant (k = 21 kyr)
test_tuned_check(k_true = 21, tbw = 2, len = 3000)
#> true k 21.00 -> obliquity 20.85, precession 20.99, obliquity + precession 20.94
# }
```
