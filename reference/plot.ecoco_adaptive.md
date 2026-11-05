# Plot an adaptive eCOCO result

Maps of correlation, global and local p, contributing periods and
normalised ridge score with the tracked sedimentation-rate ridge.

## Usage

``` r
# S3 method for class 'ecoco_adaptive'
plot(x, ridge = TRUE, reverse_depth = FALSE, ...)
```

## Arguments

- x:

  Adaptive eCOCO result to plot: the output of
  [`ecoco_adaptive`](https://stratigraphy.github.io/WaverideR/reference/ecoco_adaptive.md)
  (an object of class `ecoco_adaptive`). Required, no default.

- ridge:

  Draw the ridge-tracked sedimentation-rate path (red line) on the last
  map, the normalised ridge score (logical, TRUE or FALSE).
  `Default=TRUE`.

- reverse_depth:

  Draw depth increasing downwards, i.e. reverse the depth (y) axis of
  the maps (logical, TRUE or FALSE). `Default=FALSE`.

- ...:

  Not used; present for compatibility with the generic `plot` function.

## Value

invisible NULL
