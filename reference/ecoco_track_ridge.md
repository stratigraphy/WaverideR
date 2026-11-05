# Track a smooth sedimentation-rate ridge through an eCOCO score map

Within-window min-max normalised score, quadratic jump penalty (penalty
0.35, jump scale max(2, 4\*sr step)), Viterbi back-tracking
(trackEcocoRidge in Acycle ecoco.m).

## Usage

``` r
ecoco_track_ridge(
  score,
  sr,
  depth,
  rho,
  p,
  norbit,
  fallback = NULL,
  jump_penalty = 0.35
)
```

## Arguments

- score:

  Score map through which the ridge is tracked (numeric matrix, rates in
  rows, windows in columns), e.g. element `score` of the output of
  [`ecoco_adaptive`](https://stratigraphy.github.io/WaverideR/reference/ecoco_adaptive.md);
  it is min-max normalised within each window before tracking. Required,
  no default.

- sr:

  Sedimentation-rate grid of the score map, in cm/kyr (numeric vector,
  one value per row of `score`). Required, no default.

- depth:

  Depths of the window centres, in metres (numeric vector, one value per
  column of `score`). Required, no default.

- rho:

  Correlation coefficients (numeric matrix of the same size as `score`);
  the values along the ridge are returned. Required, no default.

- p:

  Global Monte Carlo p-values (numeric matrix of the same size as
  `score`); the values along the ridge are returned. Required, no
  default.

- norbit:

  Number of contributing orbital periods (numeric matrix of the same
  size as `score`); the values along the ridge are returned. Required,
  no default.

- fallback:

  Result returned for windows the ridge does not reach, the output of
  [`ecoco_best_per_window`](https://stratigraphy.github.io/WaverideR/reference/ecoco_best_per_window.md)
  (data frame); `NULL` computes it from the other arguments.
  `Default=NULL`.

- jump_penalty:

  Weight of the quadratic penalty on rate jumps between neighbouring
  windows (numeric, \>= 0): each step costs
  `jump_penalty * (|sr jump| / max(2, 4 * sr step))^2`. Larger values
  give a smoother ridge, `0` picks the best rate in every window
  independently. `Default=0.35`.

## Value

data.frame with depth, sr, rho, p_global, norbit, score, sr_low, sr_high
