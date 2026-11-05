# Sensitivity of every cycle to each fundamental frequency

Sensitivity of every cycle to each fundamental frequency

## Usage

``` r
term_sensitivity(
  age_ma = 0,
  delta = 0.01,
  k_yr = NULL,
  fundamentals = c("g1", "g2", "g3", "g4", "g5", "s1", "s3", "s4", "s6", "k"),
  data_dir = NULL
)
```

## Arguments

- age_ma:

  Age in Ma (numeric, single value, 0 to 4429) at which the term periods
  are computed; the precession constant is taken from
  [`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md)
  at this age unless `k_yr` is given. `Default=0`.

- delta:

  Fractional change applied to each fundamental frequency in turn
  (numeric, e.g. `0.01` = +1%; may be negative). `Default=0.01`.

- k_yr:

  Your own precession constant k, as a period in years (numeric, \> 0),
  used instead of the AstroGeo22 value at `age_ma`; `NULL` uses the
  AstroGeo22 value (`astrogeo22(age_ma)$p_arcsec`). `Default=NULL`.

- fundamentals:

  Fundamental frequencies to vary, one column of the result each
  (character vector): any of `"g1"` to `"g9"`, `"s1"` to `"s9"` (La2004
  secular frequencies) and `"k"` (precession constant). Names of a named
  vector are used as column names.
  `Default=c("g1", "g2", "g3", "g4", "g5", "s1", "s3", "s4", "s6", "k")`.

- data_dir:

  Folder that contains the data files `AstroGeo22.txt` and
  `La04gisi.csv` (character path); `NULL` uses the copy installed with
  the package. `Default=NULL`.

## Value

matrix of percent change in each term's period (rows) per \`delta\`
change in each fundamental (columns), including the precession constant
k
