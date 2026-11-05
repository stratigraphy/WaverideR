# Full Earth-Moon state from any one observable

Give exactly one of: \`p_arcsec\`/\`p_kyr\`/\`k_yr\` (precession),
\`a_RE\`/\`a_km\` (Earth-Moon distance), \`lod_h\` (solar day),
\`months_per_year\`, or an observed \`cycle_kyr\` together with its
\`term\` (e.g. the obliquity cycle).

## Usage

``` r
earthmoon_state(
  p_arcsec = NULL,
  p_kyr = NULL,
  k_yr = NULL,
  a_RE = NULL,
  a_km = NULL,
  lod_h = NULL,
  months_per_year = NULL,
  months_kind = "sidereal",
  cycle_kyr = NULL,
  term = "k+s3",
  method = c("angular_momentum", "astrogeo22"),
  obliquity_deg = NULL,
  age_ma = NULL,
  data_dir = NULL
)
```

## Arguments

- p_arcsec:

  Precession constant (rate of axial precession), in arcsec/yr (numeric,
  \> 0; about 50.47 today). If several precession inputs are given,
  `p_kyr` wins over `k_yr`/`cycle_kyr`, which win over `p_arcsec`.
  `Default=NULL`.

- p_kyr:

  Period of the axial precession, in kyr (numeric, \> 0; about 25.7
  today), converted to `p_arcsec = 360 * 3600 / (p_kyr * 1000)`.
  `Default=NULL`.

- k_yr:

  Period of the axial precession, in years (numeric, \> 0; about 25700
  today), converted to `p_arcsec = 360 * 3600 / k_yr`. `Default=NULL`.

- a_RE:

  Earth-Moon distance, in Earth radii (numeric, \> 0; about 60.14
  today). `Default=NULL`.

- a_km:

  Earth-Moon distance, in km (numeric, \> 0; about 384400 today);
  converted to `a_RE` and takes precedence over it. `Default=NULL`.

- lod_h:

  Length of the solar day, in hours (numeric, \> 0); if no distance is
  given, the Earth-Moon distance is derived from it with
  [`moon_distance_from_lod`](https://stratigraphy.github.io/WaverideR/reference/moon_distance_from_lod.md).
  `Default=NULL`.

- months_per_year:

  Number of lunar months per year (numeric, \> 0), e.g. from fossil
  growth bands; converted to an Earth-Moon distance with
  [`moon_distance_from_months`](https://stratigraphy.github.io/WaverideR/reference/moon_distance_from_months.md)
  (takes precedence over `a_RE` and `a_km`). `Default=NULL`.

- months_kind:

  Type of month counted in `months_per_year` (character): `"sidereal"`
  or `"synodic"`. `Default="sidereal"`.

- cycle_kyr:

  Period of an observed astronomical cycle in the record, in kyr
  (numeric, \> 0), converted to the precession constant with
  [`cycles_to_k`](https://stratigraphy.github.io/WaverideR/reference/cycles_to_k.md)
  using `term`; overrides `k_yr` and `p_arcsec` (but not `p_kyr`).
  `Default=NULL`.

- term:

  Which astronomical term `cycle_kyr` is (character), passed to
  [`cycles_to_k`](https://stratigraphy.github.io/WaverideR/reference/cycles_to_k.md):
  `"k+s3"` (obliquity) or a precession term `"k+g5"`, `"k+g2"`,
  `"k+g4"`, `"k+g3"` (or another term containing k listed by
  [`astro_terms`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md)).
  Only used with `cycle_kyr`. `Default="k+s3"`.

- method:

  How the Earth's spin and the lunar distance are linked when a
  precession observable has to be inverted, passed to
  [`earthmoon_from_precession`](https://stratigraphy.github.io/WaverideR/reference/earthmoon_from_precession.md)
  (character): `"angular_momentum"` (Earth-Moon angular momentum
  conserved) or `"astrogeo22"` (read off the Farhat et al. 2022
  AstroGeo22 curve; then the obliquity and model age of the first
  matching age are returned). Not used when a distance, day length or
  number of months is given.
  `Default=c("angular_momentum", "astrogeo22")`, i.e.
  `"angular_momentum"`.

- obliquity_deg:

  Obliquity of the Earth, in degrees (numeric). `NULL` uses the
  AstroGeo22 value at `age_ma` if `age_ma` is given, otherwise the
  present-day value 23.2545. `Default=NULL`.

- age_ma:

  Age in Ma (numeric), only used to take the AstroGeo22 obliquity when
  `obliquity_deg=NULL`. `Default=NULL`.

- data_dir:

  Directory that contains the AstroGeo22 data files (character, or
  `NULL`); `NULL` uses the data installed with the package.
  `Default=NULL`.

## Value

one-row-per-input data.frame with the precession constant, the
Earth-Moon distance, the solar and sidereal day, days per year, sidereal
and synodic months per year and the obliquity
