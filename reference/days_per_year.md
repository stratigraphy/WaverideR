# Days in a year for a given solar day length

Days in a year for a given solar day length

## Usage

``` r
days_per_year(lod_h)
```

## Arguments

- lod_h:

  Length of the solar day in hours (numeric vector, \> 0), e.g. from
  [`astrogeo22`](https://stratigraphy.github.io/WaverideR/reference/astrogeo22.md);
  the number of days per (Julian) year is `365.25 * 24 / lod_h`.
  Required, no default.
