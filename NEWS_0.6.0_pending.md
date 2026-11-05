# Changes in version 0.6.0

### astroMCBchron() version 2: outlier mixture, laboratory and depth offsets, tempering
* Date model `date_model = c("mixture", "two_stage", "gaussian")`. The default `"mixture"` gives every date a mixture likelihood, (1 - p) x N(date, sd) + p x outlier component, with the outlier probability p estimated (beta prior `outlier_prior`, default `c(1, 19)`, mean 0.05) and the posterior outlier probability of every date reported (`dates$p_outlier`, `outlier_prob`). The outlier component is set by `outlier_type`: `"scaled"` (Gaussian with `outlier_k` x sd), `"old"` (one-sided: inherited or recycled zircon) or `"young"` (one-sided: residual Pb loss). `"two_stage"` reproduces version 1 (Student-t screening, flagged dates with `outlier_k` x sd); `two_stage = TRUE/FALSE` is still accepted (deprecated). Own date likelihoods (`date_loglik`, e.g. from `upb_grain_loglik()`) are normalised and combined with the same mixture.
* Weighted-mean dates: optional column `dates$mswd`; with `inflate_mswd = TRUE` (default) the sd of dates with MSWD > 1 is multiplied by sqrt(MSWD).
* Systematic uncertainties: optional columns `dates$group` (laboratory or tracer) and `dates$tracer_rel`. Dates of other groups than the anchor get an estimated systematic offset with a prior from the tracer calibration uncertainty of both groups (`tracer_rel`, default 1.5e-4, EARTHTIME ET535/ET2535), so the astrochronology estimates the inter-laboratory offset. The tracer term of the reference group and the 238U decay-constant term (`decay_rel`, default 5.4e-4) are reported as extra envelopes (`summary$sd2_tracer_kyr`, `sd2_total_kyr`, `lower_tracer`, `upper_tracer`, `lower_total`, `upper_total`); they do not affect durations.
* Depth offsets: dates that share `dates$depth_group` are shifted together by one estimated offset (prior sd `depth_offset_sd`), for dates whose depths come from another depth scale than the tracks (outcrop samples correlated into a core, box-label depths).
* Parallel tempering of the astronomical part of the posterior (`temper`, inverse temperatures starting at 1): hot chains can change the number of cycles between dates, so alternative solutions are found and weighted; diagnostics in `tempering`.
* Leave-one-out: the date model is refitted in every run; `loo_ids` restricts the test to selected dates.
* `print()` and `plot()` show the outlier probabilities, group and depth offsets and the tracer and total envelopes. The misfit placement (`misfit`, `misfit_window`) and `band` options of version 1 are kept unchanged. Tests: `test-astroMCBchron-v2.R` (date models, outlier types, offsets) and `test-astroMCBchron-misfit.R`.

### astroMCBchron(): no uncertainty minimum between dates
* New argument `misfit = c("localised", "stretch")` in `astroMCBchron()`. The rigid per-interval stretch (`"stretch"`, earlier behaviour) spreads the posterior departure of an interval's duration over every cycle, so the age between two dates is a fixed interpolation between two uncertain end ages; its variance is convex and was smaller between the dates than at the dates (Tiout: 96-104 kyr between, 112-140 kyr at the dates). With `"localised"` the departure sits at one horizon of unknown position (a miscounted cycle, a condensed level), which makes the variance linear between the dated horizons; a unimodal-majorant guard covers hiatuses inside an interval and Monte Carlo noise in the 95 % band. Same MCMC posterior, mean age model, interval durations and leave-one-out. Stress test (16 random records, 3-7 dates, 0-2 hiatuses): 0 interior minima (stretch: 742), mean 2 sd +3 %. New test `test-astroMCBchron-misfit.R`.
* `misfit = "window"` (new default) with `misfit_window = 100` (kyr of floating time): the interval misfit sits at one random position, as with `"localised"`, but is spread over one contiguous window of about one short-eccentricity cycle instead of one step of the prediction grid; a negative misfit compresses the window (no step loses more than 70 % of its time). The variance is computed exactly for the window (rigid-stretch variance + 2 (E[F] - w) Cov(a_i, delta) + (E[F^2] - w^2) E[delta^2]). Altajme (paper configuration, same MCMC): 2 sd within 3 kyr of `"localised"` (mean between the outer dates 77.0 against 77.9 kyr), the guard raises the 2 sd by at most 0.1 kyr, and the stored curves keep time everywhere (curves with an empty step: 0.1 % against 96 %; largest one-step departure from the mean 16 against 55 kyr, median per curve), so they can be used for rate, phase and spectral analyses. Durations of short events lose the heavy tail of `"localised"` (e.g. Ireviken Datum 1 to Datum 4: 51 +/- 16 against 51 +/- 25 kyr). `"localised"` and `"stretch"` remain available.
* New argument `band = c("analytic", "curves")`: the 95 % band and the median in `summary` are by default the analytic mean +/- 1.96 sd (no Monte Carlo noise in the band); the percentiles of the curves are returned as `median_mc`, `lower_mc`, `upper_mc`, and `guard_raise_kyr` reports how much the guard raised the 2 sd. `test-astroMCBchron-misfit.R` now also tests `"window"` (no interior minimum, same mean, no empty steps, within 15 % of `"localised"`, analytic band) and uses no apply-family calls. Ported to the pending v2 version as well (also the tracer and decay-constant envelopes).

### Continuous tracking of cycles with uncertainty
* New `track_ridge_hmm()`: continuous tracking of one cycle through any scalogram (wavelet, synchrosqueezed, superlet-SST, ConceFT) with a hidden Markov model (forward-backward and Viterbi): bridges gaps of sharp scalograms, gives a credible band, a confidence curve and, with `period_of_tracked_cycle`, the sedimentation rate and an age model with uncertainty from posterior paths. The output works with `curve2time()`.
* New `track_astro_hmm()`: astronomically constrained tracking of the sedimentation rate using all cycles together (eccentricity, obliquity, precession; periods and weights adjustable for deep time), optionally inside a rough guide. In synthetic tests (2000 kyr, changing sedimentation rate, AR1 noise) the age model from `analyse_superlet_sst()` + `track_astro_hmm()` was within 3-20 kyr at the end of the record (maximum error 8-25 kyr), against 55-120 kyr for picking the maximum in a band.
* `track_ridge_hmm()` and `track_astro_hmm()` accept `jump_depths` (known hiatuses or abrupt facies changes, where abrupt changes of the sedimentation rate are allowed) and `track_astro_hmm()` `background = "ar1"`. Stress tests with three jumps of the sedimentation rate and AR1 noise: largest age error about 2 % of the duration (noise half the signal), 3-4 % (noise equal to the signal), 6-9 % (noise twice the signal), about half of the error of picking the strongest power in a band; with noise three times the signal only with a rough guide (about 7 %). At high noise the credible bands are too narrow when the track follows the wrong cycle.
* New `track_astro_auto()`: the recommended workflow in one call (scalogram that keeps weak cycles, astronomically constrained tracking with deep-time periods, optional guide and jump depths) plus an independent long-window (eCOCO-like) estimate of the sedimentation rate; intervals where the two disagree or the confidence is low are flagged. In stress tests (3000 kyr, three jumps, AR1 noise) the largest age error was 32-44 kyr (noise half the signal), 67-78 kyr (equal), about 190 kyr (twice, MSST scalogram); with noise three times the signal the track failed, but the flags covered about 90 % of the wrong intervals.
* New `astro_periods_tracking()`: eccentricity, obliquity and precession periods with weights for any age (AstroGeo22 precession constant via `astro_terms()`) or a given k, e.g. for Silurian/Devonian records; `track_astro_hmm()` gained `age_ma`, `k_yr`, `long_weight`, `local_weight` and `guide_penalty`.
* Duration constraints: `track_astro_hmm()`, `track_ridge_hmm()` and `track_astro_auto()` accept `constraints` (data frame `from`, `to`, `duration` in kyr, optional `tolerance`, default 2.5 %), e.g. from radiometric ages or well-dated stage boundaries. The posterior is conditioned on the constraints (paths drawn from a tilted posterior, only those within the tolerance kept) and the best path comes from an exact constrained Viterbi (elapsed time as an extra state dimension). On the Klonk record (Přídolí, 24.8-49.5 m = 4434 kyr) the published astrochronology was reproduced within about 0.3 Myr everywhere in the Přídolí; without the constraint no setting did. Outside the constrained interval (Lochkovian) the record does not decide between the published and an about three times longer duration: the credible band contains both and the confidence is low.
* `track_astro_auto()`: `scalogram = "cwt"` (plain Morlet scalogram with `omega_nr`, e.g. 5-6), which followed the ridges more continuously than the synchrosqueezed maps where only the eccentricity is resolved. `astro_periods_tracking(include = "eccentricity")` gives the 405-kyr cycle weight 1 (instead of 0.5).
* New `analyse_chirplet_sst()` (alias `analyze_chirplet_sst()`): synchrosqueezed wavelet-chirplet transform (Chen et al. 2023; Li et al. 2026). At every scale and position the chirped Morlet wavelet whose linear frequency sweep matches the record is selected (matched chirp), the instantaneous period is estimated with it (unbiased for a locally linear chirp, i.e. a smoothly changing sedimentation rate) and the coefficient is squeezed there; several `omega_nr` are combined as superlets. On a record whose sedimentation rate doubles it put 99 % of the energy within 3 % of the true periods (superlet-SST: 95 %), and in the stress tests (three jumps of the sedimentation rate, AR1 noise equal to or twice the signal) the age models tracked from it had about 35 % lower largest errors than from `analyse_superlet_sst()`, at four times the computing time. The output also holds `Chirp`, the local rate of change of the period along the record. `track_astro_auto()` accepts `scalogram = "chirplet_sst"`.
* New `track_astro_iterative()`: demodulated (warped) tracking. The record is converted to time with a first track, transformed again in that nearly stationary time domain with sharper wavelets, the residual stretch of the time axis is tracked and the age models are composed (the demodulation principle of Meignen et al. 2017 and Jiang et al. 2022, where the age model is the demodulation). Evaluated honestly: one warped iteration lowered the largest age error by about 10-15 % in the stress tests and improved the match with the published Klonk age model inside the Přídolí from about 0.29 to 0.19 Myr; a second iteration gave nothing more, and the iteration cannot repair a first track that follows the wrong cycles (the warped record shows the same ambiguity), so prior information (`duration`, `constraints`, `jump_depths`) remains more important. The last warped scalogram (periods in kyr: the tuned scalogram) is returned.
* `track_astro_hmm()`: optional L1 smoothing of the tracked path (`smooth_path = "tv"` for a piecewise constant, `"trend"` for a piecewise linear log sedimentation rate; majorisation-minimisation, Wang et al. 2023). Evaluated: total variation made the tracks worse in the stress tests (the fraction within 5 % of the true rate fell from 75 to 61-72 %), the trend filter was neutral; both are off by default.
* New `detect_transients()`: time-reassigned (transient-extracting) transform that reports positions where the energy of many periods is reassigned to one position (spikes, event beds, sharp steps). Evaluated: it finds a spike and a step in synthetic records, but it does not see hiatuses (a gap in time is a phase jump of continuous cycles, not an impulse; the ratio at the hiatuses was indistinguishable from the median), so it cannot supply `jump_depths` automatically.
* Benchmark on the Cyclostratigraphy Intercomparison Project cases (Sinnesael et al. 2019). Case 2 (Pleistocene Imbrie-model record, sedimentation rate ramping from 1 to 1.5 cm/kyr, hiatus): every scalogram tracked the ramp within 3-5 % (chirplet-SST and MSST best, 100 % of the record within 10 %); the 200-kyr hiatus is, as expected, invisible to any tracker and is the whole remaining age error. Case 3 (Late Devonian HadSM3 snow cover, 8.7 Myr, very red noise, sedimentation rate 4.5 + 1.5 sin(2 pi t / 1.2 Myr) cm/kyr): the oscillation is recovered with the plain Morlet scalogram with omega 10-12 (`scalogram = "cwt", omega_nr = 10`; correlation with the true rate 0.89, amplitude 82 % of the true one, largest age error 72 kyr over 8.7 Myr, duration 8677 vs 8702 kyr), while the synchrosqueezed maps were worse in this very noisy record (chirplet-SST 0.65, MSST 0.51, superlet-SST 0.39: squeezing a noisy precession band gives speckle, not ridges). Two settings mattered more than the scalogram: the sedimentation-rate range must not extend below the rate at which the precession cycles fall under the sampling resolution (`consistent_cycles = TRUE` then dropped them silently and the oscillation was lost completely), and the astronomical periods must be the ones in the record (the HadSM3 forcing has precession at 18.0 and 16.5 kyr and obliquity at 31.7 kyr; with the AstroGeo22 periods for 370 Ma, 7 % longer, the duration was 8 % too long).
* New `plot_astro_evidence()`: the map that `track_astro_hmm()` tracks, i.e. the summed evidence of all astronomical cycles against sedimentation rate (over a wide, unknown range) and depth, and with `type = "posterior"` the posterior of the rate of a track (new attribute `posterior` of `track_astro_hmm()` results). This is the picture to look at before choosing `sedrate`: on CIP Case 3 with the rate unknown over 1-20 cm/kyr, the posterior of a Morlet omega-16 scalogram shows a ridge that varies between 3 and 6 cm/kyr by eye (nothing about the form of the variation is assumed or needed) (the plain track over that range: correlation 0.51 with the true rate, the track straying in two intervals where the posterior is visibly spread); tracking inside the range read off the map gives 0.84-0.89. The synchrosqueezed maps are not usable for this at this noise level: squeezing the red-noise-dominated long periods produces long horizontal streaks (the slow wander of the noise reassigned to a few periods, not cycles) that look like ridges, give a bias towards the edges of the rate range, and the short-period ridges are speckle; AR1 prewhitening removes the streaks but leaves nothing above the noise threshold. `plot_wavelet()` gained `whiten = TRUE` (power divided by the mean of its period) to show weak short-period ridges in red-noise records.
* Chirp-matched analysis, the method that recovers the oscillating sedimentation rate of CIP Case 3 with nothing assumed: new `analyse_chirp_eha()` (windowed multitaper spectra, Thomson 1982 / Meyers 2012, in which every window is first de-chirped by a linear change of the rate, u = ln(1 + g (z - z0)) / g, so that cycles whose period changes inside the window become lines again; with `criterion = "template"` the de-chirped spectrum of every window is matched to the astronomical template over a grid of rates and gradients, as in eCOCO but with a gradient) and new `track_chirp_hmm()` (hidden Markov tracking with the pair rate + gradient as the hidden state: the rate moves by the gradient, the gradient drifts slowly, so the track is steered along the slope the de-chirped evidence supports). On CIP Case 3 with the rate unknown over 1-20 cm/kyr and no prior of any kind: correlation with the true rate 0.97, 85 % of the record within 10 %, 88 % of the amplitude of the 1.2-Myr variation recovered, largest age error 131 kyr over 8.7 Myr; the tracked gradient alone correlates 0.94 with the true one. On CIP Case 2 (ramp, 0.3-5 cm/kyr unknown): 0.91, 100 % within 10 %. Two things were needed: the continuum against which a line is measured must be a smooth fit of the whole window spectrum (`continuum = "lowess"`; a running median was pulled up inside the precession triplet, so that an alias at 0.45 times the rate, where the obliquity template falls on the real precession, scored better than the truth), and the window must hold 15-25 precession cycles. The de-chirped map (`Power_norm`) with the gradient of a track (`steer`) shows the cycles as lines that follow the rate, but noise features inside a warped window are warped as well, so wavy structure in such a map is not by itself evidence. The per-window chirp without the template (`criterion = "ftest"`, `"peak"`, `"entropy"`) carried almost no information at the noise level of Case 3 (the true gradient ranked 9th of 25 candidates); only the template pooled over all cycles and the continuity of the track did. Three limits found on the way, and what was done about them: (1) abrupt jumps of the rate: a jump is not a gradient and a window straddles it, so on the synthetic records with three jumps the chirp-matched track alone was far worse than `track_astro_hmm()`; `track_chirp_hmm()` now takes `jump_prob` (the rate may jump, the gradient restarts) and `scalogram` (the local evidence of `track_astro_hmm()` from a scalogram is added to the window evidence: the sharp lines of the one, the locality of the other), which brought these records to 81-86 % within 5 % and age errors of 46-290 kyr at noise equal to the signal (the plain `track_astro_hmm()`: 75 %, 142 kyr), while at twice the noise the plain tracker is still better; (2) the warped noise of the de-chirped maps: new `chirp_eha_significance()` runs the same windows, gradients or steering on AR1 surrogates of the record and reports the line-to-continuum ratio against the 95th percentile of the warped surrogates (`Power_sig`) and a p-value map; on Case 3, 66-68 % of the positions along the tracked precession lines are significant against 2-24 % along lines shifted by 0.3 octave and 13 % of all cells; (3) the choice of windows and gradient grid: `window = NULL` takes two windows of about 12 and 25 cycles of the shortest cycle at the geometric mean of the rate range and averages their template scores, `chirp = "auto"` a gradient grid up to a factor 3 change within the longest window.
* `track_astro_auto()` now runs both trackers by default (`method = "both"`): the chirp-matched track is returned as `track`, the HMM track as `track_hmm`, with both rates and their disagreement in `check` and the fraction of agreement within 10 % in `agreement`. Tested on all records: the chirp-matched track with the scalogram evidence added was better or equal everywhere (CIP Case 3 blind: 0.97 against 0.08; Case 2: tie; synthetic jumps: 0.86-0.94 against 0.00-0.92). `jump_prob` tested as the default: 0 kept (1e-3 halved the age errors of two synthetic records and doubled those of two others, and broke Case 3). New `plot_chirp_eha()` with the two styles of enhanced spectrum that were compared (the ratio to the warped surrogates, or the de-chirped map with the significant cells outlined) and new `rate_gram()`: the template score against sedimentation rate per window, a purely spectral display in which the variation of the rate is visible without any tracking (the cyclostratigraphic counterpart of the harmonic product spectrum of pitch detection); on Case 3 the median of its maxima over 11 windows alone correlates 0.73 with the true rate (66 % within 10 %).
* New `analyse_dechirped()` and `plot_dechirped()`: the universal de-chirping. The record is converted to the time of any track, the chosen transform (Morlet wavelet, superlet-SST, multisynchrosqueezing, chirplet-SST or multitaper EHA) is computed there (the tuned scalogram, periods in kyr: where the track is right every cycle is a horizontal line, a check that does not use the truth) and mapped back to depth with the period axis in depth units at every position; with `n_sim` the same warp and transform are applied to AR1 surrogates and the significant cells are outlined. Tested on nine records (CIP Cases 2 and 3, three synthetic records, Klonk, ODP 1119, Orcadian 12/13-1, KTG-01; report `rategram_report.html`): the rate-gram plus chirp-matched tracking works with nothing assumed where precession and obliquity are above the noise and the rate changes smoothly or in a few steps (r = 0.93-0.98), and fails, visibly in the tuned scalogram, where only eccentricity is resolved (Klonk, KTG-01: alias at half or double the rate), where the noise is twice the signal, and on Orcadian; ODP 1119 is partly tracked (0.67). The failures are the records that needed a duration prior or a constraint in the earlier tests.
* `analyse_dechirped(transform = "combined")`: the geometric mean of the rank-normalised Morlet, multitaper EHA and superlet-SST maps in the time of the track. Measured on CIP Case 3 as the fraction of positions at which each astronomical line peaks at its period in the tuned map: Morlet 0.14 (weakest, steadiest at the long periods), EHA 0.21 (best at 95, 31.7 and 16.5 kyr), superlet-SST 0.19 (best at 405 kyr), MSST 0.10, chirplet-SST 0.20 (best at 18 kyr), combined 0.22 (best on average and for obliquity and 18 kyr). `plot_dechirped()` and `plot_chirp_eha()` use viridis with a quantile (histogram-equalised) colour scale by default. Benchmark report with all figures: `WaverideR_benchmark_report.pdf`.
* `spectral_moments()` now also implements the sedimentation-rate procedure of Sinnesael et al. (2018) itself (`anchor`): the moments are those of the normalised local periodogram (their equations 7-10), every spatial frequency scales with 1 / rate (equation 12), so the smoothed long-term trend of the bandwidth (or of the mean frequency, `moment`; overlap-add polynomial baseline of Verhelst 2000, `trend = "poly"`, or lowess) is anchored at one depth where an astronomical line is identified (`anchor = c(depth, period_kyr, frequency)`, or a known rate), as the paper does with the obliquity line at the top of ODP 982 and the 0.45 cycle/m precession line of U1437. The Acycle `sr_mean` path (rescaling by the first raw bandwidth value and additive re-centring on a mean rate, with a trend window that depends on the depth unit) is kept as it was, but it is not the paper's method. Two options the paper motivates but does not provide: `freq_range` (moments over a band) and `noise_floor` (median of the window periodogram subtracted), because the white-noise bias of section 2.3.3 makes the second moment follow the cutoff frequency rather than the rate: on a synthetic record (precession, obliquity, eccentricity, white noise 0.3 times the signal, rate 1-3 cm/kyr) the anchored bandwidth trend correlates 0.00 with the true rate over the whole spectrum, 0.97 over 0.1-8 cycles/m, 0.99 with the floor removed as well (91 % within 20 %); with the noise at 0.8 the band is what matters (0.70 / 0.84, the mean frequency 0.98). `plot_spectral_moments()` draws the trend and the anchored rate with the anchor marked.
* `spectral_moments()` now keeps the periodogram of every window (`spectra`: frequency x window), and `plot_spectral_moments()` draws them as an evolutive spectrum image (frequency or log2 period against depth, `yaxis`; log, quantile or linear colour scale) with the mean frequency and its bandwidth drawn on top, above the moment curves; before, only the two moment curves were shown although the spectra had been computed.
* One colour rule for all image plots (internal `.wr_colour_breaks()`): `plot_wavelet()`, `plot_superlet()`, `plot_eha_log2()`, `plot_track()`, `plot_dechirped()`, `plot_chirp_eha()`, `plot_mesa_evolutive()`, `plot_spectral_moments()` and the `flmw()` map take one argument `scale`: `"auto"` (default; logarithmic for synchrosqueezed, superlet-SST, chirplet-SST and other mostly empty maps or maps spanning more than three decades, histogram-equalised quantiles otherwise), `"log"`, `"quantile"`, `"sqrt"`, `"linear"`. Before, squeezed maps used an amplitude-linear scale under which most ridges fell in the lowest colour, and `plot_track()` drew the scalogram almost white (square root of the power on a pale palette); it now uses viridis with the power divided by its mean per period. Empty cells of de-chirped maps (outside the mapped period range) are drawn in the lowest colour instead of white. The significance styles of `plot_dechirped()` and `plot_chirp_eha()` use a diverging scale centred on the 95th percentile of the surrogates. The examples of `add_wavelet()` and `add_wavelet_avg()` drew the scalogram horizontally inside a vertical layout and overlaid an `astrochron::hilbert()` envelope of a non-centred signal (two crossing curves in the last panel); they now use `plot_horizontal = FALSE` and the wavelet amplitude of `extract_amplitude()`.
* `extract_amplitude()` replaced by the wavelet version: the analytic signal of a complex Morlet transform, either over the whole record (the wavelet equivalent of the Hilbert envelope, computed scale by scale, with AR, reflective or zero padding) or over a band around one cycle (`cycle =` a period, a tracked period curve for a varying sedimentation rate, or `"ridge"`), absolutely calibrated (a cosine of amplitude 3 returns 3 whatever `dt`, `dj`, `omega_nr`); old arguments (`pts`, `remean`, ...) are accepted and ignored with a message.
* `analyse_dechirped(n_sim > 0)`: the surrogate null distribution is pooled per period over all positions (the AR1 surrogates are stationary in depth), which is stable with a few surrogates; for synchrosqueezed transforms a ridge is compared with the ridges noise produces (the positive cells), so significance maps now exist for superlet-SST, wavelet-SST and chirplet-SST as well (before, the 95th percentile of mostly empty surrogate cells was zero).
* `plot_mesa_evolutive()`: period or frequency on a log2 axis (`log_axis = TRUE`), the common colour `scale` (log for the many decades of a maximum-entropy spectrum), viridis, and `normalise = "period"` to remove the red-noise slope. Note that the automatically chosen AR order (FPE) is low (4-6) for windows of about 100 samples, which gives one smooth peak per window; a fixed order (e.g. `order = N/3`) resolves the cycles.
* Superlet cone of influence (`analyse_superlet()`, also used by the cross superlet): the adaptive superlet uses many cycles at short periods and few at long ones, so the edge zone could be shorter at a long period than at a shorter one, which drew an hourglass-shaped cone; it is now monotonic in the period (conservative).
* `astroMCBchron()` help example rewritten: dates with 1 sigma = 80 kyr, a rate varying between 1 and 3 cm/kyr and a 405-kyr track with 2 % realisation noise; the joint model is within 21 kyr of the truth on average (95 % width 186 kyr) against 106 kyr (654 kyr) for the dates alone. The previous example (constant rate, 30-kyr dates) gave the astrochronology nothing to add.
* Function check (R CMD check, the test suite, the examples of all 221 help pages, an option sweep of the visual functions; report `WaverideR_function_check_report.pdf`). Fixed on the way: the roxygen block of `mesa_spectrum()` was attached to an internal helper, so `mesa_spectrum()` was not exported and its help page documented the wrong function; `plot_wavelet()` on a chirp-EHA result (columns at the window positions, not at the data) failed; `plot_win_fft()` without `freq_min`/`freq_max` failed; `analyse_wavelet_sst(type = "reassign")` failed on scales without reassignable coefficients; `plot_track()` failed on a chirp-matched track (no `astro_periods` attribute, no credible band) and now accepts `flags = TRUE`; `astro_anchor()` looped forever when called outside an interactive session (it waits for keyboard input) and now stops with a message; missing `@param` entries (`k`, `n_iter`, `background`) and over-long example lines in the manual.
* DESCRIPTION (CRAN review): the Description starts with what the package does instead of "Tools for", and its references are written as authors (year) <doi:...> (Torrence and Compo 1998, Moca et al. 2021, Daubechies et al. 2011, Arts et al. 2024; DOIs checked against Crossref and doi.org); `\doi{}` instead of a doi.org URL in two help pages; `.Rbuildignore` fixed (`^\.github$` and `^README\.md$` were on one line).
* Code style of the whole package (R/ and tests/): every `if`/`else`/`for`/`while`/`repeat` body is now in `{}` (about 1,470 bodies wrapped by a parser-based tool, which checked that every file still parses to the same expressions), one statement per line (no `;`), statement blocks on their own lines, 2-space block and 4-space continuation indentation (61,725 -> 64,570 lines). The apply family (`apply`, `sapply`, `lapply`, `vapply`, `Reduce`, `Filter`, `replicate`; 212 calls in R/ and 11 in the tests) is replaced by explicit `for` loops into pre-allocated objects or exact vectorised base functions (`rowSums`, `diff`, `pmax`, `max.col`); kept on purpose: `outer()` and the `parallel::mclapply()`/`parLapply()` calls inside `wr_lapply()`. The package uses no tidyverse. Checked: the objects created by the examples of all 123 help pages with examples are identical to before (tolerance 1e-12), the test suite passes.
* Graphical parameters and options are restored on exit in every function that changes them (CRAN policy): `plot.coco_adaptive()`, `plot.coco_classic()`, `plot.ecoco_adaptive()`, `plot.ecoco_classic()`, `plot.ctest()` and `plot_astro()` now save all parameters (`par(no.readonly = TRUE)`) instead of only the changed ones; `plot_win_timeOpt()` and `extract_signal(plot_residual = TRUE)` did not restore them at all; `astro_anchor()`, `delpts_tracked_period_wt()`, `dynamic_extraction()`, `ID_points()`, `tuning_pts()` and `track_chirp_hmm()` restore `options()` with `on.exit()` (before, `options(warn = -1)` stayed set when the function stopped early or was interrupted). The examples of `add_wavelet()` and `add_wavelet_avg()` reset `par()` at the end.
* `track_period_wavelet(track_peaks = FALSE)` fixed (bug present since the function was written, also in CRAN 0.5.1): the picked points were returned with depth and log2 period swapped, so the column `depth` held log2 periods, `period` held 2^depth (`Inf` for depths beyond about 1024 m) and picks at the same period row were averaged into one; `loess_auto()`, `curve2time()` and `curve2tune()` then failed ("NA/NaN/Inf in foreign function call"). Free picking now uses `locator()`: a click is snapped to the nearest depth of the wavelet and its period is taken as clicked (before, `identify()` searched every cell of the scalogram, millions of points with `dj = 1/200`); clicking a pick again removes it. `plot_horizontal` was hard-coded to `TRUE` and is now honoured (period on the x-axis, depth on the y-axis); the output is the same in both layouts. The function no longer leaves `options(warn = -1)` set. Reproduced and checked on the GR log of well GEVERIK-01 (Goeree Member, 992-1492 m); new tests with scripted clicks.
* Citation: `citation("WaverideR")` (new `inst/CITATION`), the package help page (`?WaverideR`, now also `WaverideR-package`, listed first in the package index), the DESCRIPTION and the README now ask for a double citation, the paper that introduced the package (Arts et al. 2024, Frontiers in Earth Science 12, 1357751) and the package itself (Arts 2023, CRAN).
* New `check_agemodel()`, `plot_agemodel_check()` and `as_wr_track()`: the de-chirping turned into a check of an existing age model. `as_wr_track()` builds a track from a tracked period curve (`track_period_wavelet()`), a sedimentation-rate curve or an age-depth model; `check_agemodel()` de-chirps the record with that track (Morlet, multitaper EHA, superlet-SST and chirplet-SST by default), and in every time column searches the comb of astronomical periods for the common shift that best fits the lines (a wrong rate shifts every cycle by the same factor: the harmonic-product principle of the rate-gram applied to the tuned map), which gives the sedimentation-rate correction the record asks for, per transform; the corrections are combined by a median, intervals where the combined correction exceeds `tolerance` with at least two cycles and two transforms agreeing are flagged, short flags dropped and short gaps closed (`min_length`), and a blind chirp-matched re-track is run for comparison. Each cycle is also measured on its own (`offset_single`), and `plot_agemodel_check()` shows the proxy, the tuned maps with the astronomical periods and the flagged intervals, the correction curves of all transforms and the per-cycle view. On a synthetic record whose rate steps from 5 to 8 cm/kyr under a constant-rate age model of 6.5 cm/kyr the corrections came out at 0.76 and 1.21 (truth 0.77 and 1.23) from all four transforms, with the step located within 100 kyr. Applied to seven of the user's own records (Da Silva et al. 2024 MAGS, Bisciaro Ca/Al, G14-05 GR and G17-S-01 Th, Silurian Poland GR and greyscale, Sommerodde luminosity, Altajme Ti/Al against the 405 tracks and the joint age model). The check reads the record, not the truth: a comb fit is only as good as the cycles that are expressed, and an interval flagged because the sub-eccentricity cycles sit at the wrong ratio to a well-tracked 405-kyr line says that the two do not fit together, not which of them is right.
* New `analyse_chirp_wavelet()`: the wavelet counterpart (chirped Morlet wavelets whose sweep corresponds to a rate gradient, c = -g s omega_0, the template scored on the scalogram matched to every gradient; feeds `track_chirp_hmm()`). Evaluated honestly: weaker than the EHA version (Case 3 blind: 0.57-0.77 against 0.97), because only short wavelets (precession at omega 6-16, spanning 2-4 m) can be matched by a linear chirp and they see little sweep within their own length, while the long wavelets span a whole oscillation and are capped to the plain wavelet (`chirp_max`); without the cap the steered scalogram was destroyed (a 405-kyr wavelet swept dozens of octaves). The steered wavelet scalogram (`steer`, `squeeze = TRUE`) therefore does not sharpen the lines on this record. `track_chirp_hmm()` also accepts any plain scalogram (gradient as a hidden state on the evidence of `track_astro_hmm()`): Morlet omega 10 blind 0.55 against 0.51 for the plain HMM.
* Calibration figure of what is possible: on the noise-free Case 3 series the superlet-SST and the chirplet-SST show every precession and obliquity line wiggling with the rate, so the transforms are right; the noise (red noise with lag-1 autocorrelation 0.999 plus white noise) removes that from any single transform, and the oscillation is only recovered by pooling the cycles (template) and the positions (window, continuity).
* `track_astro_hmm()`: `consistent_cycles = "partial"` (new default): cycles that leave the period range of the scalogram at one end of the sedimentation-rate range are kept, but only the changes of their evidence along the record are used (row means removed), so they still show where the rate changes without favouring the rates at which they are resolved. Cycles outside the range at every rate are still dropped (unchanged for Klonk). The weights used are returned in the attribute `cycle_weights`.
* New `plot_track()` (with `show` and `flags`).
* `track_astro_hmm()`: `consistent_cycles = TRUE` (default) uses only the cycles that stay inside the resolved period range over the whole sedimentation-rate range (before, rates that moved more cycles into the range were favoured, e.g. the ~110-kyr cycle read as obliquity in the eccentricity-only Klonk record); `track_astro_auto()` no longer drops the 405-kyr cycle at the top of the sedimentation-rate range.
* Tests on six real records (ODP 1119; Orcadian Basin 12/13-1; KTG-01; Klonk; Bisciaro; Pod Barrandovem) led to: evidence normalised by a running median along the period axis (`normalise = "median"`, removes the red-noise bias), cycles outside the resolved period range no longer favoured, `duration` argument of `track_astro_auto()` (approximate duration as a wide prior on the mean sedimentation rate; without such a prior the tracks were often a factor ~2 off, with it the durations agreed within ~15 % with the published or manual age models in 9 of 12 runs), flags from the spread of the posterior (`confidence10`) instead of the long-window estimate, and `time_median` (posterior median age) in the output.


### Even more wavelets, higher-order and superlet synchrosqueezing
* 11 further wavelet families (22 in total): exactly zero-mean Morlet (`"morlet0"`), real Morlet (`"morl"`), Hermitian hat (`"hhat"`), Poisson, gammatone, Hann (raised cosine), and the orthogonal Haar, Daubechies (`"db1"`-`"db20"`), Symlets (`"sym2"`-`"sym20"`) and Battle-Lemarie spline wavelets, computed in the frequency domain from their filters (Daubechies and Symlet filters by spectral factorisation). PyWavelets style names (`"db4"`, `"sym8"`, `"cgau3"`, `"gaus2"`, `"cmor1.5-1"`, `"shan2-1.5"`, `"fbsp2-1-1.5"`) are accepted, and `wavelet` can be a user defined function `psi_hat(omega)`.
* `analyse_wavelet_sst()`: `sst_order = 2` (second-order synchrosqueezing, Oberlin et al. 2015, 2017) corrects the instantaneous frequency for the chirp rate, so ridges follow fast changes of the sedimentation rate also with long wavelets; `type = "set"` (synchroextracting transform, Yu et al. 2017); several `omega_nr` values give adaptive synchrosqueezing (per position the wavelet with the lowest Renyi entropy, Sheu et al. 2017).
* `analyse_wavelet_sst()` also gained `n_iter` (multisynchrosqueezing, Yu et al. 2019: the frequency reassignment is repeated) and `type = "reassign"` (reassigned scalogram, Auger and Flandrin 1995: power moved in period and in depth/time).
* New `analyse_superlet_sst()`: synchrosqueezed transforms of wavelets with an increasing number of cycles (default omega_0 = 8, 12, 16; second order, 4 reassignment iterations) combined by their geometric mean (superlet principle, Moca et al. 2021). In synthetic tests (eccentricity-modulated precession and obliquity, changing sedimentation rate, AR1 noise) it put 75-90 % of the energy within 3 % of the true periods, against 25-40 % for the Morlet scalogram, and separated the 19 and 23.7 kyr precession bands best.
* New `analyse_conceft()`: ConceFT multitaper synchrosqueezing with orthonormal generalised Morse wavelets (Daubechies, Wang and Wu 2016), and `method = "multitaper"`.
* New `sst_significance()`: significance of synchrosqueezed transforms against AR noise (optionally with a changing variance) by simulation with the same settings; non-significant coefficients are removed.

### More wavelets, wavelet agnostic functions, SST defaults, FFT uncertainty
* `analyse_wavelet()` now offers 11 mother wavelets: Morlet, Paul, DOG (Mexican hat = order 2), Morse, bump, complex Morlet (`"cmor"`), analytic Gaussian derivative (`"cgau"`), Shannon, frequency B-spline (`"fbsp"`), Meyer and log-Gabor; with aliases (`"gabor"`, `"mexh"`, `"ricker"`, `"gmw"`, `"cauchy"`, ...).
* `analyse_Xwavelet()`, `analyse_wavelet_coherence()`, `analyse_Awavelet()`, `analyse_superlet()` and `analyse_Xsuperlet()` gained `wavelet` and `wavelet_order`. In the adaptive wavelet and the superlets the number of cycles is converted to the Morlet equivalent bandwidth of the chosen wavelet; the Morlet results are numerically unchanged.
* `analyse_wavelet_sst()` works with every analytic wavelet (all except DOG). Default `omega_nr = 12`, noise threshold `"riskshrink"` (universal threshold on the finest scales), amplitude linear colour levels in plots.
* New `fft_frequency_uncertainty()` (practical 1/T and Heisenberg-Gabor 1/(4 pi T) frequency uncertainty, relative uncertainty) and `fft_uncertainty()` (uncertainty of a tracked cycle and its time model, Fourier counterpart of `wavelet_uncertainty()`).

Wavelets
* `analyse_wavelet()` has new arguments `wavelet` and `wavelet_order`: besides the Morlet wavelet (default, unchanged) the Paul, derivative-of-Gaussian (DOG, e.g. Mexican hat), generalised Morse and bump wavelets can be used. The period axis is always the Fourier period. The reconstruction constant and the cone of influence follow the wavelet, so `extract_signal*()`, `wavelet_significance()` and the plotting functions work with the new wavelets. The cross-wavelet, coherence and Awavelet functions are still Morlet only.
* New `analyse_wavelet_sst()` (alias `analyze_wavelet_sst()`): wavelet synchrosqueezing transform (Daubechies et al. 2011; Thakur et al. 2013). The result has the structure of `analyse_wavelet()`, so `plot_wavelet()`, tracking and extraction work on it; ridges of cycles are much thinner and the transform stays invertible.

New functions
* `mesa_spectrum()` and `plot_mesa_spectrum()`: maximum entropy spectral analysis (MESA, Burg). Burg recursion in base R (identical coefficients to `stats::ar.burg`), AR order chosen with FPE, AIC, MDL or CAT or set by hand, one-sided power spectral density with the variance preserved, and peak table with prominence.
* `mesa_evolutive()`, `track_mesa()`, `mesa_average()` and `plot_mesa_evolutive()`: evolutive MESA with one or several window lengths (for example 30, 40 and 50 m), padding of the window at the top and bottom (`pad`) so that the spectrum spans the entire record, peak extraction at every window position, dynamic-programming tracking of a spectral ridge, average amplitude spectrum, and spectral-image plot. `mesa_spectrum()` also accepts a single `window` with padding.
* `test_tuned()`, `plot_test_k()` and `test_tuned_check()`: updated precession-constant test for tuned records (version 3). MTM peaks (optionally of several proxies pooled) are assigned to k+s3, k+s6, k+g5, k+g2 and k+g4 over a wide k range (default 19.5-23.5 kyr) instead of windows around AstroGeo22, with a pre-scan that removes intervals with little precession, ambiguity handling of peak assignments, combined k with uncertainty, an Earth-Moon distance and length-of-day conversion, and two further tests (spectral peaks, amplitude modulation against surrogates). Works better for Silurian/Devonian data.

British English
* WaverideR now uses British spelling as its primary spelling: `analyse_wavelet()`, `analyse_superlet()`, `analyse_Awavelet()`, `analyse_Xwavelet()`, `analyse_Xsuperlet()`, `analyse_wavelet_coherence()` and `wr_coloured_noise()`. The American names (`analyze_*`, `wr_colored_noise()`) remain available as identical aliases that share the same help page. Arguments with American spelling also accept British synonyms (`colour_brewer`, `colour_grid`, `standardise`, `centre_limits`, `normalise`); the original argument names keep working.

Removed (overlap with 'astrochron')
* `berger_periods()` removed (overlaps with `astrochron::bergerPeriods()`); `compare_published()` no longer contains the Berger et al. row and the data file `Berger1989Periods.txt` was removed.
* `Hilbert_transform()` is no longer exported; use `astrochron::hilbert()`. A padded internal helper is still used by `extract_amplitude()`.

Parallel computing
* foreach/doSNOW replaced by an internal base-R `parallel` helper. Sequential by default; `run_multicore = TRUE` plus the new `n_cores` argument switch parallel on. Cores are capped (physical cores - 1, max 4, tasks, free RAM, CRAN limit) and work is done in batches with incremental combining so memory stays bounded. Results are reproducible with `set.seed()` and independent of the number of cores.
* Forked workers on Linux/macOS, PSOCK workers on Windows and in RStudio (WaverideR must be installed for PSOCK workers).

Plot windows
* New `wr_open_device()`: no more Windows-only code. Interactive tracking (`identify()`) works on Windows, macOS (XQuartz/Quartz) and Linux (X11) and gives a clear message when no interactive device is available.
* `dev_new` now defaults to TRUE only in the RStudio GUI (separate window needed for raster plots); override with `options(waverider.dev_new = ...)`.

Integrated functions
* All acycle2R functions (astrochronology, Earth-Moon, COCO/eCOCO, etc.), with data in `inst/extdata`.
* `astroMCBchron()` and `modBchron()` with print/summary/predict/plot methods and a memory guard.

Dependencies
* Removed: foreach, doSNOW, DescTools, Hmisc, Matrix, matrixStats, reshape2, truncnorm, colorednoise, trapezoid, fANCOVA, DecomposeR, rlist, RColorBrewer, colorRamps, viridis (base-R replacements). scico and magick are now optional (Suggests).

Fixes
* `track_ecoco()`: new interactive tracker for eCOCO results, the counterpart of `track_period()`. Plot the eCOCO map (`pcoco`, `rho`, ...), click ridge points (click again to deselect) and get the tracked sedimentation rate against depth. Points can also be given as a table (`points =`) for scripted, reproducible tracking, and the result can be interpolated to all windows.
* `analyze_superlet()` was checked against a direct convolution implementation of Moca et al. (2021) and corrected: (1) the default settings and `order_scaling = "none"` multiplied the order by `c1` (and `o = NULL` used order `c1` instead of 1); (2) `order_scaling = "power"` used the wrong operator precedence and `"quadratic"` was documented but missing; (3) records with a non-zero mean or a trend (typical for stratigraphic proxies) gave power at the longest periods and large edge effects, so a linear trend (or mean) is now removed and the record standardised (`detrend`, `standardize`; use `detrend = "none", standardize = FALSE` for the old behaviour); (4) amplitudes are computed from sums of logarithms (no underflow for high orders) and phases as circular means; (5) a phase sign convention that was opposite to `analyze_wavelet()`; (6) new output `coi` (cone of influence) and `order`, and a warning when `lowerPeriod` is below the Nyquist period. (7) the additive cycle rule (`mult = FALSE`) is now checked against Moca et al. (2021): the i-th wavelet has `c1 + i - 1` cycles (it was `c1 + i`, so order 1 differed from the multiplicative case); additive results change accordingly.
* `analyze_Xsuperlet()`: `order_scaling` and `order_alpha` were not passed on to `analyze_superlet()`; the default `upperPeriod`/`lowerPeriod` values were swapped; the cross spectrum was divided by the period although superlet amplitudes are not scale dependent; the cross phase was taken from a complex geometric mean, whose phase is not defined (it wraps), and is now the difference of the circular-mean phases. Records are now trimmed to their common, complete part.
* New `wavelet_significance()`: significance levels for wavelet and superlet scalograms. `method = "analytic"` is the instantaneous Torrence and Compo (1998) red-noise test (local chi-square levels and the global spectrum with its degrees of freedom; wavelets only; about 4-5 percent false positives in simulations with Morlet omega 6, 8 and 10). `method = "surrogate"` runs noise surrogates through the same transform (wavelets and superlets) with white, AR(1), signal-corrected AR(1) (candidate periods removed before the noise is estimated) or power-law backgrounds, gives the scale-wise 95 percent level outside the cone of influence and the global spectrum level (after the Weesp analysis), and optionally p-values for candidate cycles. Surrogates run on several cores with results independent of the number of cores. `plot_wavelet()` and `plot_superlet()` have the new argument `add_signif` to outline the significant cells and draw the global level.
* New `wavelet_band_significance()`: significance of the scale-averaged wavelet power (Torrence and Compo 1998, section 4d) over a period band such as the obliquity or precession band, against a white or AR(1) background, instantaneous and only at points where the whole band lies outside the cone of influence. `method = "exact"` (default) derives the null distribution from the eigenvalues of the band-averaged quadratic form, so it is valid for any `omega_nr`, `dj`, band and background (checked against Monte Carlo and 4.6-5.1 percent false positives for omega 4, 6, 10 and 12). `method = "torrence"` is the published chi-square form; its constant (decorrelation length 0.60) only holds for omega 6 and makes the test too conservative for larger omega (about 2 percent false positives at omega 10), so it is rescaled to 0.6 * 6 / omega_nr by default. Includes `print` and `plot` methods.
* The scalogram plots (`plot_wavelet()`, `plot_Awavelet()`, `plot_eha_log2()`, `plot_superlet()`, `plot_Xwavelet()`, `plot_Xsuperlet()`, `plot_wavelet_coherence()`) now share one set of panel functions (`R/plot_wavelet_panels.R`) instead of 8 to 16 hand-copied layout branches each: the seven files shrank from about 17,700 to about 2,000 lines (plus 400 lines of shared panels), and the colour palette code that was copied into 17 functions is one helper. The output was checked pixel by pixel against the previous version for several hundred option combinations (data, average, MTM panels, horizontal/vertical, direction, period limits, lines, points, ablines, p-values, significance contours, arrows); it is identical except for the bug fixes below.
* Plot bug fixes found while doing this: `color_brewer = "viridis"` and `"colorRamps"` failed in several plot functions (a stray space in the generated function name); the `heat.colors`, `terrain.colors`, `topo.colors` and `cm.colors` palettes failed (they have no start/end arguments); `analyze_wavelet(pval = TRUE)` stopped with an error; `plot_eha_log2(plot_opt = "Power")` used a non-existing field for the average; the vertical `plot_eha_log2()` layout without extra panels drew no period axis, frame or user lines; the vertical `plot_Xwavelet()` layout with data panel and no average failed; the cross spectrum and superlet plots with data and average panels drew a stray y axis in the global spectrum panel; p-value overlays redrew axes. `plot_eha_log2()` and `analyze_wavelet()` lost their large commented-out old copies.
* `analyze_superlet()`: the cone of influence followed the integer order of the widest wavelet and therefore had steps; it now follows the fractional order and is smooth. The wavelets are made per frequency instead of all being kept in memory, which lowers the memory use for long records, and `run_multicore`/`n_cores` compute the frequencies on several cores (identical results). The object also stores `mult`, `detrend` and `standardize`.
* `flmw()` accepts superlet objects.
* `plot_superlet()`, `plot_Xsuperlet()` and `track_period()` now shade the cone of influence of superlet objects (`coi.1`, `coi.2`), in both orientations, as for wavelets.
* `flmw()`: classic method rewritten with vectorised code, about 50-60 times faster with identical output (checked against the previous version, including `run_random = TRUE`). The old code stopped with "missing value where TRUE/FALSE needed" for some spectra; this no longer happens. New arguments `method` ("classic", "weighted", "template"), `tolerance` and `whiten`; results gain `score_mat` and `best_sedrate`.
* New `flmw_path()`: continuity along depth (Viterbi path) and posterior probabilities with credible intervals, optional sedimentation rate bounds and hiatuses. In a synthetic test (6 cycles, noise equal to the signal, sedimentation rate changing 3-6 cm/kyr) the root mean square error fell from about 2.0 cm/kyr (classic, best score per depth) to about 0.15-0.2 cm/kyr (template or weighted score plus path).
* `model_red_noise_wt()` computes the Morlet wavelets once and transforms all scales with one `mvfft()` call (about 2 to 4 times faster, identical output). `sum_power_sedrate()` is rewritten without per-position loops over data frames: band limits and local maxima are computed once and window maxima come from a sparse table, about 20 times faster with identical output (checked on several settings).
* `wavelet_significance()`: the critical level of the global spectrum was `NA` (red dashed line stopped) at periods without any point outside the cone of influence, while the global spectrum itself was drawn there. The analytic level now uses the instantaneous level (2 degrees of freedom) there and the surrogate level uses all points, so the line is continuous over the whole period range; `gws` (observed) is still `NA` there.
* Found with a plot gallery of all examples and layouts: `plot_eha_log2()` failed when `lowerPeriod`/`upperPeriod` were not given (it looked for a non-existing `$Period`; they now default to the range of the object). The phase arrows of `plot_Xwavelet()`, `plot_Xsuperlet()` and `plot_wavelet_coherence()` had fixed lengths in depth and period units (`T_unit = 15`, `P_unit = 0.5`), which is far too large for short records and distorts the angle; the default is now automatic (`T_unit = NULL`, `P_unit = NULL`: 0.3 inch whatever the record length, with the screen angle equal to the phase angle). Colour key labels no longer repeat (0.001, 0.001, 0.000) for small powers, and the p-value axis of the global spectrum no longer collides with the colour key.
* `analyze_superlet()`: the left and right cone-of-influence polygons overlapped (darker band, crossing lines) at long periods, where the edge zones are wider than half the record. The cone of influence is now one polygon, whose sides meet in the middle of the record like the dome of the wavelet plots (`coi` itself is unchanged).
* `plot_astro_anchor()`: the example used `time_dir = FALSE` for anchor points in which time increases with depth, so all tie lines crossed in one point; it now uses `time_dir = TRUE` (set `time_dir = FALSE` only for sections where time decreases with depth).
* New `astro_anchor_interactive()`: a more flexible version of `astro_anchor()`. Proxy record and astronomical solution are shown together with their peaks; tie points are set by clicking a peak in each record (in either order, in any place of the record), removed by clicking a tied peak again, and the windows can be zoomed and moved. *Tune and inspect* opens a second window at any moment with the record tuned to the current tie points, a cycle extracted from it, the astronomical solution (with the correlation) and the sedimentation rate between the tie points; tie points in the wrong order are flagged. Tie lines connect the points in both windows; partners outside the other window are shown as stubs with their position, and *Sync* moves one window to the partners of the tie points in the other. The result has the same columns as `astro_anchor()`. Existing tie points can be passed in `ties`, and `locator_fun` allows scripted use.
* New `astro_anchor_app()`: a Shiny app (runs locally) for tying a proxy record to an astronomical solution. Proxy record, astronomical solution and the record tuned to the current tie points (with an extracted cycle, the astronomical solution and the sedimentation rate) are shown in one window, joined by numbered tie lines, and update while tie points are made, selected, moved, untied (double click) and undone. Astronomical tie points always sit on a maximum/minimum; proxy tie points can be placed on the raw record or on a cycle extracted with a Taner filter (`flow`, `fhigh`, `roll = 10^x`), need not lie on peaks, and can be moved with fine/coarse buttons, by typing a position or an offset from the nearest peak, by clicking or by dragging. The astronomical target can be a provided solution (or a list of them), a metronome (cosine with chosen period, phase and interval) or a Laskar solution (La2004 ETP/eccentricity/obliquity/precession, La2010a-d, La2011, insolation; downloaded on request) over a chosen time interval. Both records have window sliders, number boxes, mouse-wheel zoom and Shift+drag panning, and both axes can be flipped and the time direction reversed (reference frame). Needs the suggested package `shiny`.
* `add_MTM = TRUE` in `plot_wavelet()` ("object 'MTM_res_1' not found") is covered by a regression test.
* `retrack_wt_MC()`: duplicated sequential branch removed.
* Random draws from truncated normals use an inverse-CDF method and are not bit-identical to earlier versions for a given seed.

# Changes in version 0.5.1

removed issues resulting in notes
Resolved issues
Possibly misspelled words in DESCRIPTION:
Found the following (possibly) invalid URLs:
Examples with CPU (user + system) or elapsed time > 5s

# Changes in version 0.5.0

maintained code
added functions

# Changes in version 0.4.0

maintained code
added functions

# Changes in version 0.3.3

fixed an error in: track_period_wavelet


# Changes in version 0.3.2

added functions:

add_wavelet
adds a wavelet plot to a larger plot

add_wavelet_avg
adds the average wavelet power to a larger plot

curve2time_unc
converts a tracked period curve to time taking the uncertainty of the curve 
into account

curve2time_unc_anchor
converts a tracked period curve to time by taking the uncertainty of the curve 
and anchored radiometric date into account 

dur_gaps
calculates the duration off a hiatus based on the missing part of an astronomical
cycle

geo_col
generate standard color codes for the Geological Time Scale

geo_loc
generates ages for the boundaries of a geochronological subdivision

geo_mid
generate the mean age of a geological subdivision

lithlog_disc
discriticizes lithologs from a block pattern to a continous signal which
allows it to analyzed using spectral analysis techniques

minimal_tuning
Creates an age model using the minimal tuning technique.
This means that the distance between 2 peaks of an extracted cycle is 
set to duration of the interpreted astronomical cycle

retrack_wt_MC
Re-track cycles using a Monte-Carlo simulation
Then input is an already tracked period of an astronomical cycle and the wavelet
in which it was tracked. the already tracked periods are combined to create a 
pseudo-curve which is then re-tracked in one of the wavelets (randomly selected)
the result is a bunch of retracked curves 

sedrate2tune
converts a data set using the sedimentation rate to a tuned age-model

updated functions: 

Hilbert_transform
added side padding to avoid edge effects

plot_wavelet
added additional color plotting option 

plot_win_fft
added additional color plotting option 

track_period_wavelet
added additional color plotting option 
added the option to select any point within the plot

since submission:
## CRAN package check 
Kurt Hornik writes:

Still open ...

Please fix before 2023-09-10 to safely retain your package on CRAN.

Dear maintainer,

R-devel will soon find

Found the following (possibly) invalid DOIs:
  DOI: 1016/j.epsl.2018.02.010
    From: man/mag.Rd
    Status: 404
    Message: Not Found
  DOI: <doi:10.1190/1.1441328>
    From: man/analyze_wavelet.Rd
          man/plot_wavelet.Rd
    Status: 404
    Message: Not Found
  DOI: <doi:10.1190/1.1441329>
    From: man/analyze_wavelet.Rd
    Status: 404
    Message: Not Found
  DOI: doi.org/10.1051/0004-6361/201117504
    From: man/astro_anchor.Rd
    Status: 404
    Message: Not Found

for you package: previously DOIs in Rd files were missed in the CRAN
submission checks.

Please fix.

When doing so, please note that for the valid DOIs you do things like

win_fft.Rd:time series analysis <\doi{doi:10.1016/j.earscirev.2018.11.015}>

where you enclose the \doi{} inside <...>: please don't do this (and
drop the extra doi:).  One needs to do e.g.

  <doi:10.1016/j.earscirev.2018.11.015>

 in the package DESCRIPTION Description, but should do

  \doi{10.1016/j.earscirev.2018.11.015}

in Rd.

Please fix before 2023-08-23.

Best
-k

The (possibly) invalid links were fixed
The way links should be written in the documentation and documentation
is also adjusted to be inline with CRAN guidelines


# Changes in version 0.3.1

Examples contained dontrun{} which were replaced with donttest{} when examples are executed in >5 seconds
otherwise dontrun{} is removed completely. Searched all files using CTR:+SHIFT+F to see if dontrun{} was still
present and dontrun{} was not found after replacement

Checked and can now confirm that only 1 core is used by default in functions which can run
using multiple cores by default the parameter to run multiple cores is set to (run_multicore = FALSE).
In examples run_multicore is set to; (run_multicore = FALSE).

Information messages were still written to the console which is replaced with the if(verbose)cat(..) option
for interactive scripts where printing text to the console is important text is added to the examples and arguments of text indicating that functionality is reduced if verbose is set to FALSE meaning that only basic functionality is
maintained

package names, software names and API are now written in quotes e.g. '----'.

some plotting functions opened new windows for new plots, these plot would not show up in Rmarkdown documents.
An option was added to turn off the open plots in new windows. 

added support for more color palette options.
added support to switch the x and y-axis when plotting.


# CRAN feedback to version 0.3.0

Please always write package names, software names and API (application 
programming interface) names in single quotes in title and description. 
e.g: --> 'WaverideR'
Please note that package names are case sensitive.

You still have example wrapped within dontrun{} which I think is not 
necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace 
dontrun{} with donttest{} or explain why dontrun{} is necessary for 
these examples.

You still write information messages to the console that cannot be 
easily suppressed. It is more R like to generate objects that can be 
used to extract the information a user is interested in, and then 
print() that object.
Instead of print()/cat() rather use message()/warning()  or 
if(verbose)cat(..) (or maybe stop()) if you really have to write text to 
the console.
(except for print, summary, interactive functions)

Please confirm that you do not use more than 2 cores in your examples, 
vignettes, etc.


# Changes in version 0.3.0

General changes: 

* Added References both in the description file of the package and the description texts of individual functions.

* Changed print() statements to if(verbose)cat(..) to allow suppression of messages in the console
An exception was made for the astro_anchor function where cat() is used instead because because "feedback" needs
to entered into the console

* spell-checked code and removed misspelled words 

* Removed dontrun from examples and replaced with donttest when examples are executed in >5 seconds
otherwise dontrun is removed and extra code is added. Many examples still have donttest added to them as the examples of functions use real-world examples of already published records. Using real-world example data-sets is preferred over synthetic data sets as  

* Previous iterations of functions often contained code which changed the user's options. by adding oldpar <- par(no.readonly = TRUE) followed by  on.exit(par(oldpar)) to functions this issue is avoided. 
This however makes images/plots uneditable after plotting. To keep images/plot editable the parameter "keep_editable" is added to functions which generate images/plots, If keep_editable is set to TRUE plot/images remain editable but by default keep_editable is set to FALSE to avoid changing the user's options.

* Functions which  have the option to run on  multiple cores have the parameter "run_multicore" added to them. On default run_multicore is set to FALSE to avoid running the process on multiple cores overloading the user's computer. In all examples run_multicore is set to FALSE meaning that only 1 core is used when a function is executed. 

Function specific changes:

* added astrosignal_example data set
* added depth_rank_example data set 
* added extract_amplitude function 
* added flmw function function
* added lithlog_disc function
* added plot_sed_model function
* added plot_win_fft function
* added win_fft  function
* major revision of the max_detect and min_detect functions which reduced the possibility of generating erroneous output 



# WaverideR 0.2.0

* Added a `NEWS.md` file to track changes to the package.
