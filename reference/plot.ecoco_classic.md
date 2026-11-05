# Plot an eCOCO result (Li et al. 2018 method)

Maps of correlation, H0 significance, contributing parameters and
window-normalised score with the tracked sedimentation-rate ridge.

## Usage

``` r
# S3 method for class 'ecoco_classic'
plot(x, ridge = TRUE, reverse_depth = FALSE, ...)
```

## Arguments

- x:

  eCOCO result to plot: the output of
  [`ecoco_classic`](https://stratigraphy.github.io/WaverideR/reference/ecoco_classic.md)
  (an object of class `ecoco_classic`). Required, no default.

- ridge:

  Draw the tracked sedimentation-rate ridge (red line) on the last map,
  the window-normalised score (logical, TRUE or FALSE). `Default=TRUE`.

- reverse_depth:

  Draw depth increasing downwards, i.e. reverse the depth (y) axis of
  the maps (logical, TRUE or FALSE). `Default=FALSE`.

- ...:

  Not used; present for compatibility with the generic `plot` function.

## Value

invisible NULL
