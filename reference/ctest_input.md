# Build the CTest input table from interpreted cycle wavelengths

Build the CTest input table from interpreted cycle wavelengths

## Usage

``` r
ctest_input(
  age,
  e405,
  e125,
  e95,
  obl = c(NA, NA),
  p1 = c(NA, NA),
  p2 = c(NA, NA),
  p3 = c(NA, NA)
)
```

## Arguments

- age:

  Age range of the record in Ma (numeric vector `c(min, max)`); the
  secular frequencies are sampled within it in
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md).
  Required, no default.

- e405, e125, e95:

  Wavelength ranges of the 405-kyr, ~125-kyr and ~95-kyr eccentricity
  cycles in the depth-domain spectrum, in metres (numeric vectors
  `c(min, max)`; the order of the two values does not matter). Required,
  no default.

- obl:

  Wavelength range of the obliquity cycle (k+s3), in metres (numeric
  vector `c(min, max)`); `c(NA, NA)` when no obliquity is interpreted
  (then use `mode="P"` in
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)).
  `Default=c(NA, NA)`.

- p1, p2, p3:

  Wavelength ranges of the precession cycles, in metres (numeric vectors
  `c(min, max)`), in the row order of the published spreadsheet, i.e. by
  decreasing period: `p1` the longest (k+g5, ~23.7 kyr today), `p2` k+g2
  (~22.4 kyr), `p3` the shortest (k+g4 and k+g3 together, ~19 kyr);
  `c(NA, NA)` when not interpreted (then use `mode="T"` in
  [`ctest`](https://stratigraphy.github.io/WaverideR/reference/ctest.md)).
  Note: the comments inside the published functions label these rows the
  other way round, but the code tests them in the order used here (row 6
  against k+g5, row 7 against k+g2, row 8 against k+g4 and k+g3).
  `Default=c(NA, NA)` for each.

## Value

the 8x2 matrix used by ctest(), in the published layout (column 1 =
lower frequency = longer wavelength)
