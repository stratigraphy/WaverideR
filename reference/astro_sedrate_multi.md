# Combined astronomical test for several proxies and sections

Runs
[`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md)
on every series (proxy) and fit and combines the evidence: within a
section all series share one sedimentation rate, so their same-rate
evidence is combined rate by rate; sections have their own rates and are
combined at the end.

## Usage

``` r
astro_sedrate_multi(
  data = NULL,
  section = NULL,
  sedmin = 1,
  sedmax = 20,
  numsed = 60,
  fit = c(1, 2),
  combine = c("minp", "stouffer"),
  scores = c("TimeOpt", "peaks"),
  weights = NULL,
  age_ma = NULL,
  targetE = NULL,
  targetP = NULL,
  targetO = NULL,
  n_null = 200,
  null_smooth = 2,
  run_multicore = FALSE,
  n_cores = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  The series to test (named list; each element an evenly spaced
  two-column matrix or data frame with depth in metres in column 1 and
  proxy value in column 2). Series of the same section must share the
  depth scale. Unnamed elements are named `"series1"`, `"series2"`, ...
  Required, no default (`Default=NULL` gives an error).

- section:

  Section label of every series (character or numeric vector, same
  length as `data`); series with the same label share one sedimentation
  rate and are combined rate by rate. `NULL` puts all series in one
  section. `Default=NULL`.

- sedmin, sedmax:

  Lowest and highest trial sedimentation rate in cm/kyr (numeric, \> 0):
  one value for all sections, or one value per section (named by section
  label, or in the order of `unique(section)`). `Default=1` and
  `Default=20`.

- numsed:

  Number of trial rates per section (integer, \> 0), log-spaced unless
  `linLog = 0` is passed through `...`. `Default=60`.

- fit:

  TimeOpt fits run for every series (numeric vector with `1` (precession
  modulated by eccentricity), `2` (short eccentricity modulated by the
  405-kyr cycle) or both); every fit of every series is a separate
  component. `Default=c(1, 2)`.

- combine:

  How the components of one section are combined at each rate
  (character): `"minp"` takes the best component, `"stouffer"` sums
  their correlation-corrected standard-normal scores (see Details).
  `Default=c("minp", "stouffer")`, i.e. `"minp"`.

- scores:

  Evidence used per series and fit, entered as separate components
  (character, one or both of `"TimeOpt"` (the TimeOpt significance
  S_opt) and `"peaks"` (the tuned-band peak score S_peaks)).
  `Default=c("TimeOpt", "peaks")`.

- weights:

  Weights of the series in the Stouffer combination (named numeric
  vector, names as in `data`); series not named get weight 1. `NULL`
  gives all series equal weight. `Default=NULL`.

- age_ma:

  Age of the record in Ma (numeric) for the AstroGeo22 target periods;
  passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md).
  `NULL` uses present-day periods. `Default=NULL`.

- targetE, targetP, targetO:

  Your own eccentricity, precession and obliquity target periods in kyr
  (numeric vectors), passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md);
  `NULL` uses the La2004/AstroGeo22 periods. `Default=NULL`.

- n_null:

  Number of noise series simulated per series and fit (integer, \> 0);
  passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md).
  `Default=200`.

- null_smooth:

  Number of neighbouring rates on each side pooled for the same-rate p
  and smoothing of the noise statistics (integer, \>= 0); passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md)
  and used for the combined score. `Default=2`.

- run_multicore:

  Simulate the noise series in parallel on several cores (logical, TRUE
  or FALSE); passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md).
  `Default=FALSE`.

- n_cores:

  Number of cores for the parallel computation (integer); `NULL` chooses
  automatically (see
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md)).
  `Default=NULL`.

- verbose:

  Print a progress message for every series and fit (logical, TRUE or
  FALSE). `Default=FALSE`.

- ...:

  Further arguments passed to
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md),
  e.g. `linLog`, `coco`, `tune` or `null_fit` (not `keep_null`,
  `sensitivity` or `target_sensitivity`, which are set by this
  function).

## Value

A list of class `astro_sedrate_multi` with

- sections:

  per section: data frame of rate, combined score S_multi (-log10
  same-rate p), same-rate p, family-wise p and the number of components
  used

- best:

  per section: best rate, S_multi and family-wise p

- components:

  per series and fit: best rate and family-wise p of its own score

- p_all:

  Fisher combination of the section family-wise p-values

- runs:

  the `astro_sedrate_test` results

## Details

**Which tests enter.** For every series and fit the TimeOpt significance
(S_opt) and the tuned-band peak score (S_peaks) of `astro_sedrate_test`
enter as separate components (not their joint score, which loses power
on short records); a component enters the combination at a rate only
where it is testable there (enough samples per carrier cycle and enough
modulating cycles after the edge trim). Tests that cannot work on a
series (for example the 405-kyr fit on a record that spans fewer than
two 405-kyr cycles at that rate) therefore do not dilute the evidence of
the others. Weights can be given per series (for example to down-weight
a proxy known to be diagenetically overprinted); they must be chosen
before looking at the results.

**Combination within a section.** `combine = "stouffer"` sums the
standard-normal scores of the same-rate p-values, z = sum(w z_i) /
sqrt(w' R w), where R is the correlation between the components'
z-curves over the testable rates (clamped to 0-1; proxies of one section
share noise, and fit 1 and fit 2 of one proxy share data, so the
components are not independent). This is the most powerful combination
when several proxies carry the signal. `combine = "minp"` (default)
takes the best component at each rate (most powerful when only some
proxies or tests carry the signal, the usual case; it is calibrated with
independent noise series per component and is therefore slightly
conservative for correlated proxies). Either statistic is calibrated
with the noise series of the components (same-rate p, rates within
`null_smooth` steps pooled, exponential tail), and the family-wise p
over the rates of the section uses the minimum-p method.

**Combination of sections.** Section-level family-wise p-values are
combined with Fisher's method (sections are independent records): p_all.

Independent age control (biostratigraphic duration bounds,
magnetostratigraphy, radioisotopic ages) is given through `sedmin` and
`sedmax` per section; a narrower, independently justified range reduces
the number of rates tested and so raises the attainable significance
(Kemp, 2016).

## References

Kemp, D.B. (2016). Optimizing significance testing of astronomical
forcing in cyclostratigraphy. Paleoceanography 31, 1516-1531
\<doi:10.1002/2016PA002963\>

Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. and Williams,
R.M. (1949). The American Soldier, Vol. 1: Adjustment during Army Life.
Princeton University Press.

Westfall, P.H. and Young, S.S. (1993). Resampling-Based Multiple
Testing: Examples and Methods for p-Value Adjustment. Wiley, New York,
ISBN 978-0-471-55761-6.

## Examples

``` r
# \donttest{
set.seed(1)
d1 <- mag[1:300, ]
d2 <- cbind(d1[, 1], d1[, 2] + rnorm(300, sd = sd(d1[, 2])))
res <- astro_sedrate_multi(list(a = d1, b = d2), sedmin = 0.4, sedmax = 1.5, numsed = 12,
    fit = 1, n_null = 10, coco = FALSE)
res$best
#>   section best_sedrate  S_multi      p_fw n_series
#> 1       1     1.330168 1.321453 0.4545455        2
# }
```
