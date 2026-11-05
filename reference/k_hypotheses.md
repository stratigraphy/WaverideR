# Which precession term is that peak? One k per hypothesis.

With only one or two precession cycles the assignment matters more than
the measurement. This returns the precession constant and Earth-Moon
state implied by reading a single cycle as each of the candidate terms.

## Usage

``` r
k_hypotheses(
  period_kyr,
  age_ma,
  terms = c("k+g5", "k+g2", "k+g4", "k+g3", "k+s3"),
  method = c("angular_momentum", "astrogeo22"),
  data_dir = NULL
)
```

## Arguments

- period_kyr:

  Observed period of the precession (or obliquity) cycle in kyr
  (numeric, one value \> 0). Required (no default).

- age_ma:

  Approximate age of the record in Ma (numeric, one value); used for the
  obliquity of the angular-momentum inversion and for the AstroGeo22
  model value that the column `dk_vs_model_kyr` is compared with.
  Required (no default).

- terms:

  Candidate astronomical terms the cycle may represent (character
  vector; each must contain k), e.g. `"k+g5"`, `"k+g2"`, `"k+g4"`,
  `"k+g3"`, `"k+s3"`; one output row per term (names of a named vector
  become row names).
  `Default=c("k+g5", "k+g2", "k+g4", "k+g3", "k+s3")`.

- method:

  How the precession constant is converted to Earth-Moon distance and
  length of day (character, passed to
  [`earthmoon_state`](https://stratigraphy.github.io/WaverideR/reference/earthmoon_state.md)):
  `"angular_momentum"` (conservation of Earth-Moon angular momentum) or
  `"astrogeo22"` (inversion of the AstroGeo22 curve of Farhat et al.,
  2022). `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- data_dir:

  Folder that contains the data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path). `NULL` uses `getOption("cocoR.data")`
  or the copies installed with WaverideR. `Default=NULL`.

## Value

data.frame, one row per candidate term
