# Earth-Moon distance, length of day and obliquity for CTest samples

Each sampled precession constant is converted with
precession_frequency() (Farhat et al. 2022, eq. A.7), either by
inverting the AstroGeo22 curve or by assuming Earth-Moon
angular-momentum conservation.

## Usage

``` r
ctest_earthmoon(
  x,
  age_ma = NULL,
  method = c("angular_momentum", "astrogeo22"),
  obliquity_deg = NULL,
  probs = c(0.025, 0.5, 0.975),
  data_dir = NULL
)
```

## Arguments

- x:

  Precession constants to convert: a `ctest` object from
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)
  (its samples of `k` are used) or a numeric vector of precession
  constants given as periods in years. Required, no default.

- age_ma:

  Age of the record in Ma (numeric), used for the obliquity of the
  angular-momentum inversion (when `obliquity_deg` is `NULL`) and for
  the comparison with AstroGeo22 at that age (`model_at_age`). `NULL`
  uses the mean of the age range of a `ctest` object, or no age for a
  numeric vector (then the present-day obliquity is used and no
  comparison is made). `Default=NULL`.

- method:

  How the precession constant is converted to Earth-Moon distance and
  length of day (character): `"angular_momentum"` assumes conservation
  of the Earth-Moon angular momentum at the given obliquity;
  `"astrogeo22"` reads distance, length of day, obliquity and the
  equivalent model age off the AstroGeo22 curve (last 1000 Myr).
  `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- obliquity_deg:

  Obliquity in degrees (numeric) used in the angular-momentum inversion;
  `NULL` takes the AstroGeo22 value at `age_ma` (or the present-day
  23.25 degrees when no age is known). Ignored by `method="astrogeo22"`.
  `Default=NULL`.

- probs:

  Probabilities of the quantiles reported in the summary next to the
  mean (numeric vector, values between 0 and 1).
  `Default=c(0.025, 0.5, 0.975)`.

- data_dir:

  Folder that contains the astronomical data table `AstroGeo22.txt`
  (character path); `NULL` uses the table installed with the package (or
  the folder set with `options(cocoR.data = ...)`). `Default=NULL`.

## Value

list with \`summary\`, \`samples\`, \`method\`, \`age_ma\` and, when an
age is known, \`model_at_age\`
