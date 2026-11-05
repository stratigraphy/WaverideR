# Smooth sedimentation rate curve with uncertainty from a flmw score matrix

Turns the sedimentation rate x depth score matrix of
[`flmw`](https://stratigraphy.github.io/WaverideR/reference/flmw.md)
into a continuous sedimentation rate curve. The sedimentation rate is
not allowed to jump around from one depth to the next: the function
finds the most probable path through the score matrix (Viterbi
algorithm) and also the probability of every sedimentation rate at every
depth (forward-backward algorithm of a hidden Markov model). The
probabilities give credible intervals, and prior bounds or known
hiatuses can be supplied.

## Usage

``` r
flmw_path(
  x,
  jump_sd = 0.2,
  sharpness = 3,
  normalize = TRUE,
  sedrate_range = NULL,
  hiatus = NULL,
  probs = c(0.025, 0.16, 0.84, 0.975),
  genplot = FALSE,
  x_lab = "depth",
  y_lab = "sedrate",
  normalise = NULL
)
```

## Arguments

- x:

  Output of
  [`flmw`](https://stratigraphy.github.io/WaverideR/reference/flmw.md)
  (any `method`), or any list with the elements `score_mat` (matrix,
  trial sedimentation rates in rows, depths in columns), `depth` and
  `y_axis` (trial sedimentation rates). Required (no default).

- jump_sd:

  Smoothness of the curve: standard deviation of the relative
  (natural-log) change in sedimentation rate over 1 metre of depth,
  scaled with the square root of the depth step (random walk) (numeric,
  \> 0). Small values give a stiff curve, large values let the rate
  change quickly. `Default=0.2`.

- sharpness:

  How strongly the score is trusted (numeric, \> 0): the likelihood of a
  sedimentation rate is `exp(sharpness * score)` with scores between 0
  and 1. Lower values give wider credible intervals. `Default=3`.

- normalize:

  Divide the scores by their 99th percentile (and cap them at 1) so that
  `sharpness` means the same for every `flmw` method (logical, TRUE or
  FALSE). `Default=TRUE`.

- sedrate_range:

  Prior bounds `c(min, max)` of the sedimentation rate in cm/kyr
  (numeric vector of length 2): trial rates outside this range are
  excluded. `NULL` uses all trial rates. `Default=NULL`.

- hiatus:

  Depths of hiatuses or other known breaks, in the depth unit of `x`
  (numeric vector). At these depths the sedimentation rate may jump
  freely. `NULL` means no breaks. `Default=NULL`.

- probs:

  Cumulative probabilities of the credible interval limits that are
  returned for every depth (numeric vector, values between 0 and 1).
  `Default=c(0.025, 0.16, 0.84, 0.975)` (95 and 68 percent intervals).

- genplot:

  Plot the posterior probability of the sedimentation rate against depth
  with the most probable path (logical, TRUE or FALSE). `Default=FALSE`.

- x_lab:

  Label of the x-axis (depth) of the plot (character).
  `Default="depth"`.

- y_lab:

  Label of the y-axis (sedimentation rate) of the plot (character).
  `Default="sedrate"`.

- normalise:

  British spelling of `normalize` (logical, TRUE or FALSE); if given
  (not `NULL`) it overrides `normalize`. `Default=NULL`.

## Value

A list with `depth`, `y_axis` (trial sedimentation rates), `path` (most
probable sedimentation rate curve), `mean` (posterior mean), `interval`
(matrix with the credible interval limits for every depth), `posterior`
(matrix of probabilities, rates x depths) and `score_mat`.

## Examples

``` r
# synthetic score matrix with a gradual change in sedimentation rate
depth <- seq(0, 20, by = 0.1)
rates <- seq(1, 6, by = 0.25)
true_rate <- 2 + 3 * depth / 20
score <- outer(rates, true_rate, function(r, t) exp(-0.5 * ((r - t) / 0.4)^2))
set.seed(1)
score <- pmin(pmax(score + matrix(rnorm(length(score), 0, 0.2), nrow(score)), 0), 1)
fake <- list(score_mat = score, depth = depth, y_axis = rates)
res <- flmw_path(fake)
head(res$path)
#> [1] 2 2 2 2 2 2
```
