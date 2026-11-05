## Resubmission

This is a resubmission. In response to the review:

* The Description no longer starts with "Tools for"; it starts with what the package does.
* References in the Description are now written as authors (year) <doi:...> without spaces, for the methods
  (Torrence and Compo 1998; Moca et al. 2021; Daubechies et al. 2011) and the paper that introduced the package
  (Arts et al. 2024); the package itself is linked as <https://CRAN.R-project.org/package=WaverideR>. Every DOI
  was checked against Crossref and doi.org.
* Graphical parameters and options: every function that changes par() or options() now saves the old settings
  and calls on.exit() on the very next line, before the change (plot.coco_adaptive(), plot.coco_classic(),
  plot.ecoco_adaptive(), plot.ecoco_classic(), plot.ctest(), plot_astro(), plot_win_timeOpt(), extract_signal(),
  astro_anchor(), delpts_tracked_period_wt(), dynamic_extraction(), ID_points(), tuning_pts(), track_chirp_hmm()
  were changed). No function calls setwd(). The examples that change par() reset it at the end. The package
  checks this with a script that runs every example and compares par() and options() before and after each call.

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
