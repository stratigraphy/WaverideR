# Predict ages at new positions from a modBchron (or modifiedBChron) model

This function is used to predict the age of new stratigraphic points
using the output from an ageModel run.

## Usage

``` r
modBchron_predict(
  model,
  newPositions,
  ids = 1:length(newPositions),
  newPositionThicknesses = rep(0, length(newPositions)),
  probability = 0.95
)
```

## Arguments

- model:

  Age-depth model, the output of
  [`modBchron`](https://stratigraphy.github.io/WaverideR/reference/modBchron.md)
  (a list with at least `model`, `predictPositions` and `burn`).
  Required (no default).

- newPositions:

  Stratigraphic positions whose age is predicted, in the unit of the
  `positions` of the model (numeric vector); they must lie within the
  `predictPositions` range of the original model run, outside it the
  prediction is `NA`. Required (no default).

- ids:

  Names of the new positions, used as row names of the output (character
  or numeric vector, same length as `newPositions`).
  `Default=1:length(newPositions)`.

- newPositionThicknesses:

  Stratigraphic uncertainty of each new position, given as half
  thickness in the unit of `newPositions` (numeric vector, \>= 0, same
  length and order as `newPositions`); in every MCMC iteration the
  position is drawn uniformly within plus or minus this value.
  `Default=rep(0, length(newPositions))` (no positional uncertainty).

- probability:

  Probability of the returned highest density interval (numeric, between
  0 and 1). `Default=0.95`.

## Value

HDI = 95 percent Highest Density Interval for each `newPosition`

raw = Age predictions for each `newPosition` for each MCMC run.
