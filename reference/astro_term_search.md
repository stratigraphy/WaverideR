# Which term could a measured period be?

Enumerates single terms and low-order combinations and returns those
whose period matches, so an unidentified peak can be checked against the
theory.

## Usage

``` r
astro_term_search(
  period_kyr,
  age_ma = 0,
  tol = 0.03,
  max_order = 2,
  k_yr = NULL,
  data_dir = NULL
)
```

## Arguments

- period_kyr:

  The measured period of the unidentified peak, in kyr (numeric, single
  value \> 0). Required, no default.

- age_ma:

  Age in Ma (numeric, single value, 0 to 4429) used to take the
  precession constant of the k-bearing terms from
  [`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md)
  (ignored when `k_yr` is given). `Default=0`.

- tol:

  Relative tolerance for a match (numeric, \> 0): a term is returned
  when \|period - `period_kyr`\| / `period_kyr` \<= `tol`, e.g. 0.03 =
  3%. `Default=0.03`.

- max_order:

  Highest order of the candidate terms (numeric): `1` tests single terms
  only, `2` (or more) also their harmonics 2(term) and the sums and
  differences of two terms. `Default=2`.

- k_yr:

  Your own precession constant k, as a period in years (numeric, \> 0),
  used instead of the AstroGeo22 value; `NULL` takes it from
  `astrogeo22(age_ma)$p_arcsec`. `Default=NULL`.

- data_dir:

  Folder that contains the data files `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the copy installed with
  the package. `Default=NULL`.

## Value

data.frame of candidate terms sorted by \|offset from the measurement\|
