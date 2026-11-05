# Evaluate Adaptive COCO for a batch of spectra (cocoAdaptiveEvaluate.m)

Evaluate Adaptive COCO for a batch of spectra (cocoAdaptiveEvaluate.m)

## Usage

``` r
coco_adaptive_evaluate(
  data_power,
  data_frequency,
  n_target,
  spacing,
  pad,
  orbit9,
  sr_grid,
  method = "Pearson",
  max_frequency = NULL,
  target_model = c("coherent-nine", "phase-averaged"),
  amplitude_mode = c("four-group-area", "per-orbit", "fixed"),
  observed_p = FALSE
)
```

## Arguments

- data_power:

  One-sided periodograms of the data to evaluate (numeric matrix, one
  spectrum per column, with `floor(max(pad, n_target)/2) + 1` rows on
  the frequency grid `data_frequency`); negative values are set to zero.
  Required, no default.

- data_frequency:

  Frequencies of the rows of `data_power`, in cycles per metre (numeric
  vector, evenly spaced from 0). Required, no default.

- n_target:

  Number of samples of the analysed record or window (integer); the
  target sinusoids are synthesised on this many points. Required, no
  default.

- spacing:

  Depth sampling interval of the analysed record, in metres (numeric, \>
  0). Required, no default.

- pad:

  Length (number of points) to which the record was zero-padded for the
  periodogram (integer); the target spectra are padded to
  `max(pad, n_target)`. Required, no default.

- orbit9:

  Astronomical target periods in kyr (numeric vector; nine distinct
  values are needed for `amplitude_mode="four-group-area"` and for the
  grouped weights of `"fixed"`). Required, no default.

- sr_grid:

  Sedimentation rates to evaluate, in cm/kyr (numeric vector, \> 0).
  Required, no default.

- method:

  Correlation coefficient between the data and target spectra
  (character): `"Pearson"` or `"Spearman"`. `Default="Pearson"`.

- max_frequency:

  Highest frequency of the target that enters the correlation, in
  cycles/kyr (numeric, must be at least `max(1/orbit9)`); `NULL` uses
  `1.2 * max(1/orbit9)`. `Default=NULL`.

- target_model:

  How the target spectrum is synthesised (character): `"coherent-nine"`
  (spectrum of the sum of the sinusoids, with their interference) or
  `"phase-averaged"` (sum of the individual spectra, averaged over sine
  and cosine phase). `Default=c("coherent-nine", "phase-averaged")`,
  i.e. `"coherent-nine"`.

- amplitude_mode:

  How the target amplitudes are set (character): `"four-group-area"`
  (one amplitude per group of long eccentricity, short eccentricity,
  obliquity and precession, fitted to the band areas; Acycle default for
  Adaptive COCO/eCOCO; needs `target_model="coherent-nine"`),
  `"per-orbit"` (each amplitude fitted to the spectral peak at its own
  frequency, Acycle adaptive9a) or `"fixed"` (weights 1.0/0.8/0.6 from
  [`coco_fixed_target_weights`](https://stratigraphy.github.io/WaverideR/reference/coco_fixed_target_weights.md)).
  `Default=c("four-group-area", "per-orbit", "fixed")`, i.e.
  `"four-group-area"`.

- observed_p:

  Also compute the parametric p-value of the correlation for the first
  column of `data_power` at every rate (logical, TRUE or FALSE); used
  for the observed spectrum, not for the Monte Carlo spectra. `FALSE`
  returns `NA` p-values. `Default=FALSE`.
