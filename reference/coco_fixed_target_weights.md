# Fixed sinusoid-amplitude weights for a COCO target (cocoFixedTargetWeights.m)

Fixed sinusoid-amplitude weights for a COCO target
(cocoFixedTargetWeights.m)

## Usage

``` r
coco_fixed_target_weights(orbit)
```

## Arguments

- orbit:

  Target periods in kyr (numeric vector). With nine distinct periods the
  weights follow their rank (longest five, i.e. eccentricity: 1.0;
  sixth, obliquity: 0.8; shortest three, precession: 0.6); otherwise
  periods of 80 kyr and longer get 1.0, 30 to 80 kyr 0.8 and shorter
  ones 0.6. Required, no default.
