# Check the precession-constant test on a synthetic record

Builds synthetic records (five orbital terms, a 405 kyr cycle and AR(1)
noise) with a known precession constant and prints the k recovered by
[`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md).
A check that the chosen `tbw` resolves the terms.

## Usage

``` r
test_tuned_check(
  k_true = c(20.7, 21, 21.2, 21.4),
  tbw = 2,
  noise = 1.5,
  seed = 1,
  len = 5500,
  ...
)
```

## Arguments

- k_true:

  Known precession constants of the synthetic records, in kyr (numeric
  vector); one record is built and tested for each value.
  `Default=c(20.7, 21.0, 21.2, 21.4)`.

- tbw:

  MTM time-bandwidth product passed to
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md)
  (numeric, \> 0); see `tbw` there. `Default=2`.

- noise:

  Standard deviation multiplier of the AR(1) noise (coefficient 0.7)
  added to the orbital signal (numeric, \>= 0); the orbital terms have
  amplitudes up to 1. `Default=1.5`.

- seed:

  Seed of the random number generator (integer), reset before each
  synthetic record so that the phases and noise can be reproduced.
  `Default=1`.

- len:

  Length of the synthetic records in kyr, sampled every 1 kyr (numeric,
  \> 0). `Default=5500`.

- ...:

  Further arguments passed on to
  [`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md)
  (e.g. `age_ma`, `k_search`); `tbw`, `n_surr = 2` and `plot = FALSE`
  are already set.

## Value

Prints one line per `k_true`.
