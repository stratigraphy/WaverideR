## Release summary

This is version 0.6.0, a larger update of WaverideR 0.5.1. Main changes: the functions of acycle2R
and astroMCBchron are integrated, parallel computing uses base R `parallel` (sequential by default,
cores and memory are capped), plotting code was refactored, and new significance functions and
a vignette were added. See NEWS.md.

## Dependencies

Imports were reduced to astrochron and base packages (foreach, doSNOW, DescTools, Hmisc, Matrix,
matrixStats, reshape2, truncnorm, colorednoise, trapezoid, fANCOVA, DecomposeR, rlist, RColorBrewer,
colorRamps and viridis were replaced by internal base-R code). scico and magick are in Suggests and used conditionally.

## Test environments

* local R CMD check --as-cran (Linux): see results below
* win-builder (R-devel): submitted; results are sent by e-mail to the maintainer

## R CMD check results

0 errors | 0 warnings | 0 notes expected on CRAN machines. Locally on a machine without the
required locale there is one locale warning and notes about URL checking and timing;
these are environmental.

Examples that take longer than a few seconds are in \donttest{}. Tests that run longer are skipped on CRAN.
