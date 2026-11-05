# Prepare a depth series for COCO (cocoPrepareRegularData.m)

Removes non-finite rows, sorts by depth, averages duplicate depths and
linearly interpolates to the median step when the spacing is uneven.

## Usage

``` r
coco_prepare_regular_data(raw, verbose = TRUE, min_points = 4)
```

## Arguments

- raw:

  Record to prepare (matrix or data frame; only the first two columns
  are used): depth in metres in column 1 and the proxy in column 2.
  Required, no default.

- verbose:

  Print a message when the data are interpolated to an even spacing
  (logical, TRUE or FALSE). `Default=TRUE`.

- min_points:

  Minimum number of unique, finite depth points required (integer);
  fewer gives an error. `Default=4`.
