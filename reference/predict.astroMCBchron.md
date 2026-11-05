# Predict ages at new positions from an astroMCBchron model

Predict ages at new positions from an astroMCBchron model

## Usage

``` r
# S3 method for class 'astroMCBchron'
predict(object, newPositions, side = c("upper", "lower"), ...)
```

## Arguments

- object:

  an `astroMCBchron` object.

- newPositions:

  positions (same type and unit as the model).

- side:

  for positions exactly at a hiatus: "upper" (younger side, default) or
  "lower".

- ...:

  unused.

## Value

list with `raw` (draws x positions, in age units) and `HDI` (position,
2.5, 50 and 97.5 percentiles).
