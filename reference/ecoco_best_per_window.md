# Best sedimentation rate per eCOCO window

Row of maximum finite score in every window (bestPerWindowSummary in
Acycle ecoco.m).

## Usage

``` r
ecoco_best_per_window(sr, depth, rho, p, norbit, score)
```

## Arguments

- sr:

  Sedimentation-rate grid of the score map, in cm/kyr (numeric vector,
  one value per row of the matrices). Required, no default.

- depth:

  Depths of the window centres, in metres (numeric vector, one value per
  column of the matrices). Required, no default.

- rho:

  Correlation coefficients (numeric matrix, rates in rows, windows in
  columns), e.g. element `rho` of the output of
  [`ecoco_adaptive`](https://stratigraphy.github.io/WaverideR/reference/ecoco_adaptive.md);
  the value at the best rate is returned. Required, no default.

- p:

  Global Monte Carlo p-values (numeric matrix of the same size as
  `rho`), e.g. element `p_global`; the value at the best rate is
  returned. Required, no default.

- norbit:

  Number of contributing orbital periods (numeric matrix of the same
  size as `rho`). Required, no default.

- score:

  Score used to choose the best rate in each window, the rate with the
  highest finite score (numeric matrix of the same size as `rho`), e.g.
  element `score` (pCOCO times the fraction of contributing periods).
  Windows without a finite score give `NA`. Required, no default.

## Value

data.frame with depth, sr, rho, p_global, norbit, score, sr_low, sr_high
