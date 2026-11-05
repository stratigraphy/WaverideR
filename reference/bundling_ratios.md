# Precession cycles per eccentricity bundle, predicted

Precession cycles per eccentricity bundle, predicted

## Usage

``` r
bundling_ratios(
  age_ma = 0,
  k_yr = NULL,
  ecc_terms = c("g2-g5", "g3-g2", "g4-g2", "g3-g5", "g4-g5"),
  prec_terms = c("k+g5", "k+g2", "k+g3", "k+g4", "k+s3"),
  data_dir = NULL
)
```

## Arguments

- age_ma:

  Age at which the ratios are predicted, in Ma (single numeric value, 0
  to 4429); the precession constant at that age is taken from the
  AstroGeo22 model
  ([`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md),
  linear interpolation). Ignored when `k_yr` is given. `Default=0`.

- k_yr:

  Your own value of the precession constant, given as its period in
  years (numeric, e.g. about 25700 for today), used instead of the
  AstroGeo22 value; `NULL` uses the model value at `age_ma`.
  `Default=NULL`.

- ecc_terms:

  Eccentricity terms that form the rows of the output (character
  vector). Each is a secular frequency combination such as `"g2-g5"`,
  `"g3-g2"`, combined terms such as `"2(g2-g5)"` or `"(g3-g2)-(g2-g5)"`,
  or an alias such as `"e405"`, `"e125"`, `"e95"` (see
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)).
  `Default=c("g2-g5", "g3-g2", "g4-g2", "g3-g5", "g4-g5")`.

- prec_terms:

  Precession and obliquity terms that form the columns of the output
  (character vector), each of the form `"k+gi"` (precession) or `"k+si"`
  (obliquity), or an alias such as `"prec"` or `"obl"`.
  `Default=c("k+g5", "k+g2", "k+g3", "k+g4", "k+s3")`.

- data_dir:

  Folder that contains the astronomical data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the tables installed with
  the package (or the folder set with `options(cocoR.data = ...)`).
  `Default=NULL`.

## Value

matrix of bundling ratios, eccentricity terms in rows, precession (and
obliquity) terms in columns
