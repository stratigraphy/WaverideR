# Average spectrum of an evolutive MESA

Averages the spectra of all window positions of an evolutive MESA
([`mesa_evolutive`](https://stratigraphy.github.io/WaverideR/reference/mesa_evolutive.md)):
the average amplitude spectrum, shown beside the spectral image in many
cyclostratigraphic workflows. Amplitude is the square root of the power.

## Usage

``` r
mesa_average(x, window = 1, normalise = TRUE)
```

## Arguments

- x:

  Evolutive MESA, the output of
  [`mesa_evolutive`](https://stratigraphy.github.io/WaverideR/reference/mesa_evolutive.md).
  Required (no default).

- window:

  Which window length of `x` is averaged (numeric): either the window
  length itself (as given to `mesa_evolutive`) or its number, i.e. its
  position in `x$windows`. `Default=1`.

- normalise:

  Divide the spectrum of every window position by its maximum before
  averaging, so that all positions weigh equally (logical, TRUE or
  FALSE). `Default=TRUE`.

## Value

data frame with `frequency`, `period` and `amplitude`.
