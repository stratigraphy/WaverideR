# Changelog

## Changes in version 0.6.0

#### Continuous tracking of cycles with uncertainty

- New
  [`track_ridge_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_ridge_hmm.md):
  continuous tracking of one cycle through any scalogram (wavelet,
  synchrosqueezed, superlet-SST, ConceFT) with a hidden Markov model
  (forward-backward and Viterbi): bridges gaps of sharp scalograms,
  gives a credible band, a confidence curve and, with
  `period_of_tracked_cycle`, the sedimentation rate and an age model
  with uncertainty from posterior paths. The output works with
  [`curve2time()`](https://stratigraphy.github.io/WaverideR/reference/curve2time.md).
- New
  [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md):
  astronomically constrained tracking of the sedimentation rate using
  all cycles together (eccentricity, obliquity, precession; periods and
  weights adjustable for deep time), optionally inside a rough guide. In
  synthetic tests (2000 kyr, changing sedimentation rate, AR1 noise) the
  age model from
  [`analyse_superlet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet_sst.md) +
  [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
  was within 3-20 kyr at the end of the record (maximum error 8-25 kyr),
  against 55-120 kyr for picking the maximum in a band.
- [`track_ridge_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_ridge_hmm.md)
  and
  [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
  accept `jump_depths` (known hiatuses or abrupt facies changes, where
  abrupt changes of the sedimentation rate are allowed) and
  [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
  `background = "ar1"`. Stress tests with three jumps of the
  sedimentation rate and AR1 noise: largest age error about 2 % of the
  duration (noise half the signal), 3-4 % (noise equal to the signal),
  6-9 % (noise twice the signal), about half of the error of picking the
  strongest power in a band; with noise three times the signal only with
  a rough guide (about 7 %). At high noise the credible bands are too
  narrow when the track follows the wrong cycle.
- New
  [`track_astro_auto()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md):
  the recommended workflow in one call (scalogram that keeps weak
  cycles, astronomically constrained tracking with deep-time periods,
  optional guide and jump depths) plus an independent long-window
  (eCOCO-like) estimate of the sedimentation rate; intervals where the
  two disagree or the confidence is low are flagged. In stress tests
  (3000 kyr, three jumps, AR1 noise) the largest age error was 32-44 kyr
  (noise half the signal), 67-78 kyr (equal), about 190 kyr (twice, MSST
  scalogram); with noise three times the signal the track failed, but
  the flags covered about 90 % of the wrong intervals.
- New
  [`astro_periods_tracking()`](https://stratigraphy.github.io/WaverideR/reference/astro_periods_tracking.md):
  eccentricity, obliquity and precession periods with weights for any
  age (AstroGeo22 precession constant via
  [`astro_terms()`](https://stratigraphy.github.io/WaverideR/reference/astro_terms.md))
  or a given k, e.g. for Silurian/Devonian records;
  [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
  gained `age_ma`, `k_yr`, `long_weight`, `local_weight` and
  `guide_penalty`.
- New
  [`plot_track()`](https://stratigraphy.github.io/WaverideR/reference/plot_track.md)
  (with `show` and `flags`).
- [`track_astro_hmm()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md):
  `consistent_cycles = TRUE` (default) uses only the cycles that stay
  inside the resolved period range over the whole sedimentation-rate
  range (before, rates that moved more cycles into the range were
  favoured, e.g. the ~110-kyr cycle read as obliquity in the
  eccentricity-only Klonk record);
  [`track_astro_auto()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md)
  no longer drops the 405-kyr cycle at the top of the sedimentation-rate
  range.
- Tests on six real records (ODP 1119; Orcadian Basin 12/13-1; KTG-01;
  Klonk; Bisciaro; Pod Barrandovem) led to: evidence normalised by a
  running median along the period axis (`normalise = "median"`, removes
  the red-noise bias), cycles outside the resolved period range no
  longer favoured, `duration` argument of
  [`track_astro_auto()`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md)
  (approximate duration as a wide prior on the mean sedimentation rate;
  without such a prior the tracks were often a factor ~2 off, with it
  the durations agreed within ~15 % with the published or manual age
  models in 9 of 12 runs), flags from the spread of the posterior
  (`confidence10`) instead of the long-window estimate, and
  `time_median` (posterior median age) in the output.

#### Even more wavelets, higher-order and superlet synchrosqueezing

- 11 further wavelet families (22 in total): exactly zero-mean Morlet
  (`"morlet0"`), real Morlet (`"morl"`), Hermitian hat (`"hhat"`),
  Poisson, gammatone, Hann (raised cosine), and the orthogonal Haar,
  Daubechies (`"db1"`-`"db20"`), Symlets (`"sym2"`-`"sym20"`) and
  Battle-Lemarie spline wavelets, computed in the frequency domain from
  their filters (Daubechies and Symlet filters by spectral
  factorisation). PyWavelets style names (`"db4"`, `"sym8"`, `"cgau3"`,
  `"gaus2"`, `"cmor1.5-1"`, `"shan2-1.5"`, `"fbsp2-1-1.5"`) are
  accepted, and `wavelet` can be a user defined function
  `psi_hat(omega)`.
- [`analyse_wavelet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_sst.md):
  `sst_order = 2` (second-order synchrosqueezing, Oberlin et
  al. 2015, 2017) corrects the instantaneous frequency for the chirp
  rate, so ridges follow fast changes of the sedimentation rate also
  with long wavelets; `type = "set"` (synchroextracting transform, Yu et
  al. 2017); several `omega_nr` values give adaptive synchrosqueezing
  (per position the wavelet with the lowest Renyi entropy, Sheu et
  al. 2017).
- [`analyse_wavelet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_sst.md)
  also gained `n_iter` (multisynchrosqueezing, Yu et al. 2019: the
  frequency reassignment is repeated) and `type = "reassign"`
  (reassigned scalogram, Auger and Flandrin 1995: power moved in period
  and in depth/time).
- New
  [`analyse_superlet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet_sst.md):
  synchrosqueezed transforms of wavelets with an increasing number of
  cycles (default omega_0 = 8, 12, 16; second order, 4 reassignment
  iterations) combined by their geometric mean (superlet principle, Moca
  et al. 2021). In synthetic tests (eccentricity-modulated precession
  and obliquity, changing sedimentation rate, AR1 noise) it put 75-90 %
  of the energy within 3 % of the true periods, against 25-40 % for the
  Morlet scalogram, and separated the 19 and 23.7 kyr precession bands
  best.
- New
  [`analyse_conceft()`](https://stratigraphy.github.io/WaverideR/reference/analyse_conceft.md):
  ConceFT multitaper synchrosqueezing with orthonormal generalised Morse
  wavelets (Daubechies, Wang and Wu 2016), and `method = "multitaper"`.
- New
  [`sst_significance()`](https://stratigraphy.github.io/WaverideR/reference/sst_significance.md):
  significance of synchrosqueezed transforms against AR noise
  (optionally with a changing variance) by simulation with the same
  settings; non-significant coefficients are removed.

#### More wavelets, wavelet agnostic functions, SST defaults, FFT uncertainty

- [`analyse_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
  now offers 11 mother wavelets: Morlet, Paul, DOG (Mexican hat = order
  2), Morse, bump, complex Morlet (`"cmor"`), analytic Gaussian
  derivative (`"cgau"`), Shannon, frequency B-spline (`"fbsp"`), Meyer
  and log-Gabor; with aliases (`"gabor"`, `"mexh"`, `"ricker"`, `"gmw"`,
  `"cauchy"`, …).
- [`analyse_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xwavelet.md),
  [`analyse_wavelet_coherence()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_coherence.md),
  [`analyse_Awavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Awavelet.md),
  [`analyse_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md)
  and
  [`analyse_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md)
  gained `wavelet` and `wavelet_order`. In the adaptive wavelet and the
  superlets the number of cycles is converted to the Morlet equivalent
  bandwidth of the chosen wavelet; the Morlet results are numerically
  unchanged.
- [`analyse_wavelet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_sst.md)
  works with every analytic wavelet (all except DOG). Default
  `omega_nr = 12`, noise threshold `"riskshrink"` (universal threshold
  on the finest scales), amplitude linear colour levels in plots.
- New
  [`fft_frequency_uncertainty()`](https://stratigraphy.github.io/WaverideR/reference/fft_frequency_uncertainty.md)
  (practical 1/T and Heisenberg-Gabor 1/(4 pi T) frequency uncertainty,
  relative uncertainty) and
  [`fft_uncertainty()`](https://stratigraphy.github.io/WaverideR/reference/fft_uncertainty.md)
  (uncertainty of a tracked cycle and its time model, Fourier
  counterpart of
  [`wavelet_uncertainty()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_uncertainty.md)).

Wavelets \*
[`analyse_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
has new arguments `wavelet` and `wavelet_order`: besides the Morlet
wavelet (default, unchanged) the Paul, derivative-of-Gaussian (DOG,
e.g. Mexican hat), generalised Morse and bump wavelets can be used. The
period axis is always the Fourier period. The reconstruction constant
and the cone of influence follow the wavelet, so `extract_signal*()`,
[`wavelet_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md)
and the plotting functions work with the new wavelets. The
cross-wavelet, coherence and Awavelet functions are still Morlet only.
\* New
[`analyse_wavelet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_sst.md)
(alias
[`analyze_wavelet_sst()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_sst.md)):
wavelet synchrosqueezing transform (Daubechies et al. 2011; Thakur et
al. 2013). The result has the structure of
[`analyse_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
so
[`plot_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md),
tracking and extraction work on it; ridges of cycles are much thinner
and the transform stays invertible.

New functions \*
[`mesa_spectrum()`](https://stratigraphy.github.io/WaverideR/reference/mesa_spectrum.md)
and
[`plot_mesa_spectrum()`](https://stratigraphy.github.io/WaverideR/reference/plot_mesa_spectrum.md):
maximum entropy spectral analysis (MESA, Burg). Burg recursion in base R
(identical coefficients to
[`stats::ar.burg`](https://rdrr.io/r/stats/ar.html)), AR order chosen
with FPE, AIC, MDL or CAT or set by hand, one-sided power spectral
density with the variance preserved, and peak table with prominence. \*
[`mesa_evolutive()`](https://stratigraphy.github.io/WaverideR/reference/mesa_evolutive.md),
[`track_mesa()`](https://stratigraphy.github.io/WaverideR/reference/track_mesa.md),
[`mesa_average()`](https://stratigraphy.github.io/WaverideR/reference/mesa_average.md)
and
[`plot_mesa_evolutive()`](https://stratigraphy.github.io/WaverideR/reference/plot_mesa_evolutive.md):
evolutive MESA with one or several window lengths (for example 30, 40
and 50 m), padding of the window at the top and bottom (`pad`) so that
the spectrum spans the entire record, peak extraction at every window
position, dynamic-programming tracking of a spectral ridge, average
amplitude spectrum, and spectral-image plot.
[`mesa_spectrum()`](https://stratigraphy.github.io/WaverideR/reference/mesa_spectrum.md)
also accepts a single `window` with padding. \*
[`test_tuned()`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md),
[`plot_test_k()`](https://stratigraphy.github.io/WaverideR/reference/plot_test_k.md)
and
[`test_tuned_check()`](https://stratigraphy.github.io/WaverideR/reference/test_tuned_check.md):
updated precession-constant test for tuned records (version 3). MTM
peaks (optionally of several proxies pooled) are assigned to k+s3, k+s6,
k+g5, k+g2 and k+g4 over a wide k range (default 19.5-23.5 kyr) instead
of windows around AstroGeo22, with a pre-scan that removes intervals
with little precession, ambiguity handling of peak assignments, combined
k with uncertainty, an Earth-Moon distance and length-of-day conversion,
and two further tests (spectral peaks, amplitude modulation against
surrogates). Works better for Silurian/Devonian data.

British English \* WaverideR now uses British spelling as its primary
spelling:
[`analyse_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md),
[`analyse_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md),
[`analyse_Awavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Awavelet.md),
[`analyse_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xwavelet.md),
[`analyse_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md),
[`analyse_wavelet_coherence()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet_coherence.md)
and
[`wr_coloured_noise()`](https://stratigraphy.github.io/WaverideR/reference/wr_coloured_noise.md).
The American names (`analyze_*`,
[`wr_colored_noise()`](https://stratigraphy.github.io/WaverideR/reference/wr_coloured_noise.md))
remain available as identical aliases that share the same help page.
Arguments with American spelling also accept British synonyms
(`colour_brewer`, `colour_grid`, `standardise`, `centre_limits`,
`normalise`); the original argument names keep working.

Removed (overlap with ‘astrochron’) \* `berger_periods()` removed
(overlaps with
[`astrochron::bergerPeriods()`](https://rdrr.io/pkg/astrochron/man/bergerPeriods.html));
[`compare_published()`](https://stratigraphy.github.io/WaverideR/reference/compare_published.md)
no longer contains the Berger et al. row and the data file
`Berger1989Periods.txt` was removed. \* `Hilbert_transform()` is no
longer exported; use
[`astrochron::hilbert()`](https://rdrr.io/pkg/astrochron/man/hilbert.html).
A padded internal helper is still used by
[`extract_amplitude()`](https://stratigraphy.github.io/WaverideR/reference/extract_amplitude.md).

Parallel computing \* foreach/doSNOW replaced by an internal base-R
`parallel` helper. Sequential by default; `run_multicore = TRUE` plus
the new `n_cores` argument switch parallel on. Cores are capped
(physical cores - 1, max 4, tasks, free RAM, CRAN limit) and work is
done in batches with incremental combining so memory stays bounded.
Results are reproducible with
[`set.seed()`](https://rdrr.io/r/base/Random.html) and independent of
the number of cores. \* Forked workers on Linux/macOS, PSOCK workers on
Windows and in RStudio (WaverideR must be installed for PSOCK workers).

Plot windows \* New `wr_open_device()`: no more Windows-only code.
Interactive tracking
([`identify()`](https://rdrr.io/r/graphics/identify.html)) works on
Windows, macOS (XQuartz/Quartz) and Linux (X11) and gives a clear
message when no interactive device is available. \* `dev_new` now
defaults to TRUE only in the RStudio GUI (separate window needed for
raster plots); override with `options(waverider.dev_new = ...)`.

Integrated functions \* All acycle2R functions (astrochronology,
Earth-Moon, COCO/eCOCO, etc.), with data in `inst/extdata`. \*
[`astroMCBchron()`](https://stratigraphy.github.io/WaverideR/reference/astroMCBchron.md)
and
[`modBchron()`](https://stratigraphy.github.io/WaverideR/reference/modBchron.md)
with print/summary/predict/plot methods and a memory guard.

Dependencies \* Removed: foreach, doSNOW, DescTools, Hmisc, Matrix,
matrixStats, reshape2, truncnorm, colorednoise, trapezoid, fANCOVA,
DecomposeR, rlist, RColorBrewer, colorRamps, viridis (base-R
replacements). scico and magick are now optional (Suggests).

Fixes \*
[`track_ecoco()`](https://stratigraphy.github.io/WaverideR/reference/track_ecoco.md):
new interactive tracker for eCOCO results, the counterpart of
[`track_period()`](https://stratigraphy.github.io/WaverideR/reference/track_period.md).
Plot the eCOCO map (`pcoco`, `rho`, …), click ridge points (click again
to deselect) and get the tracked sedimentation rate against depth.
Points can also be given as a table (`points =`) for scripted,
reproducible tracking, and the result can be interpolated to all
windows. \*
[`analyze_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md)
was checked against a direct convolution implementation of Moca et
al. (2021) and corrected: (1) the default settings and
`order_scaling = "none"` multiplied the order by `c1` (and `o = NULL`
used order `c1` instead of 1); (2) `order_scaling = "power"` used the
wrong operator precedence and `"quadratic"` was documented but missing;
(3) records with a non-zero mean or a trend (typical for stratigraphic
proxies) gave power at the longest periods and large edge effects, so a
linear trend (or mean) is now removed and the record standardised
(`detrend`, `standardize`; use `detrend = "none", standardize = FALSE`
for the old behaviour); (4) amplitudes are computed from sums of
logarithms (no underflow for high orders) and phases as circular means;
(5) a phase sign convention that was opposite to
[`analyze_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md);
(6) new output `coi` (cone of influence) and `order`, and a warning when
`lowerPeriod` is below the Nyquist period. (7) the additive cycle rule
(`mult = FALSE`) is now checked against Moca et al. (2021): the i-th
wavelet has `c1 + i - 1` cycles (it was `c1 + i`, so order 1 differed
from the multiplicative case); additive results change accordingly. \*
[`analyze_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_Xsuperlet.md):
`order_scaling` and `order_alpha` were not passed on to
[`analyze_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md);
the default `upperPeriod`/`lowerPeriod` values were swapped; the cross
spectrum was divided by the period although superlet amplitudes are not
scale dependent; the cross phase was taken from a complex geometric
mean, whose phase is not defined (it wraps), and is now the difference
of the circular-mean phases. Records are now trimmed to their common,
complete part. \* New
[`wavelet_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md):
significance levels for wavelet and superlet scalograms.
`method = "analytic"` is the instantaneous Torrence and Compo (1998)
red-noise test (local chi-square levels and the global spectrum with its
degrees of freedom; wavelets only; about 4-5 percent false positives in
simulations with Morlet omega 6, 8 and 10). `method = "surrogate"` runs
noise surrogates through the same transform (wavelets and superlets)
with white, AR(1), signal-corrected AR(1) (candidate periods removed
before the noise is estimated) or power-law backgrounds, gives the
scale-wise 95 percent level outside the cone of influence and the global
spectrum level (after the Weesp analysis), and optionally p-values for
candidate cycles. Surrogates run on several cores with results
independent of the number of cores.
[`plot_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)
and
[`plot_superlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_superlet.md)
have the new argument `add_signif` to outline the significant cells and
draw the global level. \* New
[`wavelet_band_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_band_significance.md):
significance of the scale-averaged wavelet power (Torrence and Compo
1998, section 4d) over a period band such as the obliquity or precession
band, against a white or AR(1) background, instantaneous and only at
points where the whole band lies outside the cone of influence.
`method = "exact"` (default) derives the null distribution from the
eigenvalues of the band-averaged quadratic form, so it is valid for any
`omega_nr`, `dj`, band and background (checked against Monte Carlo and
4.6-5.1 percent false positives for omega 4, 6, 10 and 12).
`method = "torrence"` is the published chi-square form; its constant
(decorrelation length 0.60) only holds for omega 6 and makes the test
too conservative for larger omega (about 2 percent false positives at
omega 10), so it is rescaled to 0.6 \* 6 / omega_nr by default. Includes
`print` and `plot` methods. \* The scalogram plots
([`plot_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md),
[`plot_Awavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Awavelet.md),
[`plot_eha_log2()`](https://stratigraphy.github.io/WaverideR/reference/plot_eha_log2.md),
[`plot_superlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_superlet.md),
[`plot_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xwavelet.md),
[`plot_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xsuperlet.md),
[`plot_wavelet_coherence()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet_coherence.md))
now share one set of panel functions (`R/plot_wavelet_panels.R`) instead
of 8 to 16 hand-copied layout branches each: the seven files shrank from
about 17,700 to about 2,000 lines (plus 400 lines of shared panels), and
the colour palette code that was copied into 17 functions is one helper.
The output was checked pixel by pixel against the previous version for
several hundred option combinations (data, average, MTM panels,
horizontal/vertical, direction, period limits, lines, points, ablines,
p-values, significance contours, arrows); it is identical except for the
bug fixes below. \* Plot bug fixes found while doing this:
`color_brewer = "viridis"` and `"colorRamps"` failed in several plot
functions (a stray space in the generated function name); the
`heat.colors`, `terrain.colors`, `topo.colors` and `cm.colors` palettes
failed (they have no start/end arguments);
`analyze_wavelet(pval = TRUE)` stopped with an error;
`plot_eha_log2(plot_opt = "Power")` used a non-existing field for the
average; the vertical
[`plot_eha_log2()`](https://stratigraphy.github.io/WaverideR/reference/plot_eha_log2.md)
layout without extra panels drew no period axis, frame or user lines;
the vertical
[`plot_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xwavelet.md)
layout with data panel and no average failed; the cross spectrum and
superlet plots with data and average panels drew a stray y axis in the
global spectrum panel; p-value overlays redrew axes.
[`plot_eha_log2()`](https://stratigraphy.github.io/WaverideR/reference/plot_eha_log2.md)
and
[`analyze_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_wavelet.md)
lost their large commented-out old copies. \*
[`analyze_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md):
the cone of influence followed the integer order of the widest wavelet
and therefore had steps; it now follows the fractional order and is
smooth. The wavelets are made per frequency instead of all being kept in
memory, which lowers the memory use for long records, and
`run_multicore`/`n_cores` compute the frequencies on several cores
(identical results). The object also stores `mult`, `detrend` and
`standardize`. \*
[`flmw()`](https://stratigraphy.github.io/WaverideR/reference/flmw.md)
accepts superlet objects. \*
[`plot_superlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_superlet.md),
[`plot_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xsuperlet.md)
and
[`track_period()`](https://stratigraphy.github.io/WaverideR/reference/track_period.md)
now shade the cone of influence of superlet objects (`coi.1`, `coi.2`),
in both orientations, as for wavelets. \*
[`flmw()`](https://stratigraphy.github.io/WaverideR/reference/flmw.md):
classic method rewritten with vectorised code, about 50-60 times faster
with identical output (checked against the previous version, including
`run_random = TRUE`). The old code stopped with “missing value where
TRUE/FALSE needed” for some spectra; this no longer happens. New
arguments `method` (“classic”, “weighted”, “template”), `tolerance` and
`whiten`; results gain `score_mat` and `best_sedrate`. \* New
[`flmw_path()`](https://stratigraphy.github.io/WaverideR/reference/flmw_path.md):
continuity along depth (Viterbi path) and posterior probabilities with
credible intervals, optional sedimentation rate bounds and hiatuses. In
a synthetic test (6 cycles, noise equal to the signal, sedimentation
rate changing 3-6 cm/kyr) the root mean square error fell from about 2.0
cm/kyr (classic, best score per depth) to about 0.15-0.2 cm/kyr
(template or weighted score plus path). \*
[`model_red_noise_wt()`](https://stratigraphy.github.io/WaverideR/reference/model_red_noise_wt.md)
computes the Morlet wavelets once and transforms all scales with one
[`mvfft()`](https://rdrr.io/r/stats/fft.html) call (about 2 to 4 times
faster, identical output).
[`sum_power_sedrate()`](https://stratigraphy.github.io/WaverideR/reference/sum_power_sedrate.md)
is rewritten without per-position loops over data frames: band limits
and local maxima are computed once and window maxima come from a sparse
table, about 20 times faster with identical output (checked on several
settings). \*
[`wavelet_significance()`](https://stratigraphy.github.io/WaverideR/reference/wavelet_significance.md):
the critical level of the global spectrum was `NA` (red dashed line
stopped) at periods without any point outside the cone of influence,
while the global spectrum itself was drawn there. The analytic level now
uses the instantaneous level (2 degrees of freedom) there and the
surrogate level uses all points, so the line is continuous over the
whole period range; `gws` (observed) is still `NA` there. \* Found with
a plot gallery of all examples and layouts:
[`plot_eha_log2()`](https://stratigraphy.github.io/WaverideR/reference/plot_eha_log2.md)
failed when `lowerPeriod`/`upperPeriod` were not given (it looked for a
non-existing `$Period`; they now default to the range of the object).
The phase arrows of
[`plot_Xwavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xwavelet.md),
[`plot_Xsuperlet()`](https://stratigraphy.github.io/WaverideR/reference/plot_Xsuperlet.md)
and
[`plot_wavelet_coherence()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet_coherence.md)
had fixed lengths in depth and period units (`T_unit = 15`,
`P_unit = 0.5`), which is far too large for short records and distorts
the angle; the default is now automatic (`T_unit = NULL`,
`P_unit = NULL`: 0.3 inch whatever the record length, with the screen
angle equal to the phase angle). Colour key labels no longer repeat
(0.001, 0.001, 0.000) for small powers, and the p-value axis of the
global spectrum no longer collides with the colour key. \*
[`analyze_superlet()`](https://stratigraphy.github.io/WaverideR/reference/analyse_superlet.md):
the left and right cone-of-influence polygons overlapped (darker band,
crossing lines) at long periods, where the edge zones are wider than
half the record. The cone of influence is now one polygon, whose sides
meet in the middle of the record like the dome of the wavelet plots
(`coi` itself is unchanged). \*
[`plot_astro_anchor()`](https://stratigraphy.github.io/WaverideR/reference/plot_astro_anchor.md):
the example used `time_dir = FALSE` for anchor points in which time
increases with depth, so all tie lines crossed in one point; it now uses
`time_dir = TRUE` (set `time_dir = FALSE` only for sections where time
decreases with depth). \* New
[`astro_anchor_interactive()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor_interactive.md):
a more flexible version of
[`astro_anchor()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md).
Proxy record and astronomical solution are shown together with their
peaks; tie points are set by clicking a peak in each record (in either
order, in any place of the record), removed by clicking a tied peak
again, and the windows can be zoomed and moved. *Tune and inspect* opens
a second window at any moment with the record tuned to the current tie
points, a cycle extracted from it, the astronomical solution (with the
correlation) and the sedimentation rate between the tie points; tie
points in the wrong order are flagged. Tie lines connect the points in
both windows; partners outside the other window are shown as stubs with
their position, and *Sync* moves one window to the partners of the tie
points in the other. The result has the same columns as
[`astro_anchor()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md).
Existing tie points can be passed in `ties`, and `locator_fun` allows
scripted use. \* New
[`astro_anchor_app()`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor_app.md):
a Shiny app (runs locally) for tying a proxy record to an astronomical
solution. Proxy record, astronomical solution and the record tuned to
the current tie points (with an extracted cycle, the astronomical
solution and the sedimentation rate) are shown in one window, joined by
numbered tie lines, and update while tie points are made, selected,
moved, untied (double click) and undone. Astronomical tie points always
sit on a maximum/minimum; proxy tie points can be placed on the raw
record or on a cycle extracted with a Taner filter (`flow`, `fhigh`,
`roll = 10^x`), need not lie on peaks, and can be moved with fine/coarse
buttons, by typing a position or an offset from the nearest peak, by
clicking or by dragging. The astronomical target can be a provided
solution (or a list of them), a metronome (cosine with chosen period,
phase and interval) or a Laskar solution (La2004
ETP/eccentricity/obliquity/precession, La2010a-d, La2011, insolation;
downloaded on request) over a chosen time interval. Both records have
window sliders, number boxes, mouse-wheel zoom and Shift+drag panning,
and both axes can be flipped and the time direction reversed (reference
frame). Needs the suggested package `shiny`. \* `add_MTM = TRUE` in
[`plot_wavelet()`](https://stratigraphy.github.io/WaverideR/reference/plot_wavelet.md)
(“object ‘MTM_res_1’ not found”) is covered by a regression test. \*
[`retrack_wt_MC()`](https://stratigraphy.github.io/WaverideR/reference/retrack_wt_MC.md):
duplicated sequential branch removed. \* Random draws from truncated
normals use an inverse-CDF method and are not bit-identical to earlier
versions for a given seed.

## Changes in version 0.5.1

CRAN release: 2026-04-06

removed issues resulting in notes Resolved issues Possibly misspelled
words in DESCRIPTION: Found the following (possibly) invalid URLs:
Examples with CPU (user + system) or elapsed time \> 5s

## Changes in version 0.5.0

maintained code added functions

## Changes in version 0.4.0

CRAN release: 2025-04-25

maintained code added functions

## Changes in version 0.3.3

fixed an error in: track_period_wavelet

## Changes in version 0.3.2

CRAN release: 2023-09-06

added functions:

add_wavelet adds a wavelet plot to a larger plot

add_wavelet_avg adds the average wavelet power to a larger plot

curve2time_unc converts a tracked period curve to time taking the
uncertainty of the curve into account

curve2time_unc_anchor converts a tracked period curve to time by taking
the uncertainty of the curve and anchored radiometric date into account

dur_gaps calculates the duration off a hiatus based on the missing part
of an astronomical cycle

geo_col generate standard color codes for the Geological Time Scale

geo_loc generates ages for the boundaries of a geochronological
subdivision

geo_mid generate the mean age of a geological subdivision

lithlog_disc discriticizes lithologs from a block pattern to a continous
signal which allows it to analyzed using spectral analysis techniques

minimal_tuning Creates an age model using the minimal tuning technique.
This means that the distance between 2 peaks of an extracted cycle is
set to duration of the interpreted astronomical cycle

retrack_wt_MC Re-track cycles using a Monte-Carlo simulation Then input
is an already tracked period of an astronomical cycle and the wavelet in
which it was tracked. the already tracked periods are combined to create
a pseudo-curve which is then re-tracked in one of the wavelets (randomly
selected) the result is a bunch of retracked curves

sedrate2tune converts a data set using the sedimentation rate to a tuned
age-model

updated functions:

Hilbert_transform added side padding to avoid edge effects

plot_wavelet added additional color plotting option

plot_win_fft added additional color plotting option

track_period_wavelet added additional color plotting option added the
option to select any point within the plot

since submission: \## CRAN package check Kurt Hornik writes:

Still open …

Please fix before 2023-09-10 to safely retain your package on CRAN.

Dear maintainer,

R-devel will soon find

Found the following (possibly) invalid DOIs: DOI:
1016/j.epsl.2018.02.010 From: man/mag.Rd Status: 404 Message: Not Found
DOI: <doi:10.1190/1.1441328> From: man/analyze_wavelet.Rd
man/plot_wavelet.Rd Status: 404 Message: Not Found DOI:
<doi:10.1190/1.1441329> From: man/analyze_wavelet.Rd Status: 404
Message: Not Found DOI: doi.org/10.1051/0004-6361/201117504 From:
man/astro_anchor.Rd Status: 404 Message: Not Found

for you package: previously DOIs in Rd files were missed in the CRAN
submission checks.

Please fix.

When doing so, please note that for the valid DOIs you do things like

win_fft.Rd:time series analysis \<\>

where you enclose the inside \<…\>: please don’t do this (and drop the
extra doi:). One needs to do e.g.

<doi:10.1016/j.earscirev.2018.11.015>

in the package DESCRIPTION Description, but should do

in Rd.

Please fix before 2023-08-23.

Best -k

The (possibly) invalid links were fixed The way links should be written
in the documentation and documentation is also adjusted to be inline
with CRAN guidelines

## Changes in version 0.3.1

CRAN release: 2023-06-22

Examples contained dontrun{} which were replaced with donttest{} when
examples are executed in \>5 seconds otherwise dontrun{} is removed
completely. Searched all files using CTR:+SHIFT+F to see if dontrun{}
was still present and dontrun{} was not found after replacement

Checked and can now confirm that only 1 core is used by default in
functions which can run using multiple cores by default the parameter to
run multiple cores is set to (run_multicore = FALSE). In examples
run_multicore is set to; (run_multicore = FALSE).

Information messages were still written to the console which is replaced
with the if(verbose)cat(..) option for interactive scripts where
printing text to the console is important text is added to the examples
and arguments of text indicating that functionality is reduced if
verbose is set to FALSE meaning that only basic functionality is
maintained

package names, software names and API are now written in quotes
e.g. ‘—-’.

some plotting functions opened new windows for new plots, these plot
would not show up in Rmarkdown documents. An option was added to turn
off the open plots in new windows.

added support for more color palette options. added support to switch
the x and y-axis when plotting.

## CRAN feedback to version 0.3.0

Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: –\> ‘WaverideR’ Please note that package names are case sensitive.

You still have example wrapped within dontrun{} which I think is not
necessary. Please unwrap the examples if they are executable in \< 5
sec, or replace dontrun{} with donttest{} or explain why dontrun{} is
necessary for these examples.

You still write information messages to the console that cannot be
easily suppressed. It is more R like to generate objects that can be
used to extract the information a user is interested in, and then
print() that object. Instead of print()/cat() rather use
message()/warning() or if(verbose)cat(..) (or maybe stop()) if you
really have to write text to the console. (except for print, summary,
interactive functions)

Please confirm that you do not use more than 2 cores in your examples,
vignettes, etc.

## Changes in version 0.3.0

General changes:

- Added References both in the description file of the package and the
  description texts of individual functions.

- Changed print() statements to if(verbose)cat(..) to allow suppression
  of messages in the console An exception was made for the astro_anchor
  function where cat() is used instead because because “feedback” needs
  to entered into the console

- spell-checked code and removed misspelled words

- Removed dontrun from examples and replaced with donttest when examples
  are executed in \>5 seconds otherwise dontrun is removed and extra
  code is added. Many examples still have donttest added to them as the
  examples of functions use real-world examples of already published
  records. Using real-world example data-sets is preferred over
  synthetic data sets as

- Previous iterations of functions often contained code which changed
  the user’s options. by adding oldpar \<- par(no.readonly = TRUE)
  followed by on.exit(par(oldpar)) to functions this issue is avoided.
  This however makes images/plots uneditable after plotting. To keep
  images/plot editable the parameter “keep_editable” is added to
  functions which generate images/plots, If keep_editable is set to TRUE
  plot/images remain editable but by default keep_editable is set to
  FALSE to avoid changing the user’s options.

- Functions which have the option to run on multiple cores have the
  parameter “run_multicore” added to them. On default run_multicore is
  set to FALSE to avoid running the process on multiple cores
  overloading the user’s computer. In all examples run_multicore is set
  to FALSE meaning that only 1 core is used when a function is executed.

Function specific changes:

- added astrosignal_example data set
- added depth_rank_example data set
- added extract_amplitude function
- added flmw function function
- added lithlog_disc function
- added plot_sed_model function
- added plot_win_fft function
- added win_fft function
- major revision of the max_detect and min_detect functions which
  reduced the possibility of generating erroneous output

## WaverideR 0.2.0

- Added a `NEWS.md` file to track changes to the package.
