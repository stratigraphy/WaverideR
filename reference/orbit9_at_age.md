# Nine orbital periods (kyr) at age_ma

The nine periods used by Acycle's COCO/eCOCO GUI (see
calculate_orbit9()).

## Usage

``` r
orbit9_at_age(age_ma, data_dir = NULL)
```

## Arguments

- age_ma:

  Age at which the nine periods are calculated, in Ma (single finite
  numeric value, 0 to 4429); it is rounded to the nearest whole Ma (see
  [`calculate_orbit9`](https://stratigraphy.github.io/WaverideR/reference/calculate_orbit9.md)).
  Required (no default).

- data_dir:

  Folder that contains the data tables `AstroGeo22.txt` and
  `La04gisi.csv` (character path). `NULL` uses `getOption("cocoR.data")`
  or the copies installed with WaverideR. `Default=NULL`.

## Value

numeric vector of nine periods (kyr).
