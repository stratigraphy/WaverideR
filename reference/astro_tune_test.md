# Test a minimal tuning: multi-band rate tracking with surrogates run through the whole chain

Estimates the sedimentation-rate history of one or more stratigraphic
records from all astronomical cycles at once, turns it into a
(minimally) tuned age model with uncertainty, and tests the tuned
records for astronomical amplitude and frequency modulation against
surrogates that keep the power spectrum of the data. The rate tracker
here is a fast version of
[`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
(Morlet power, no synchrosqueezing), so that every surrogate can be
tracked, tuned and tested in exactly the same way as the data; this is
what calibrates the gain that any tuning gives to any series. For the
best age model of a record use
[`track_astro_auto`](https://stratigraphy.github.io/WaverideR/reference/track_astro_auto.md),
[`track_astro_hmm`](https://stratigraphy.github.io/WaverideR/reference/track_astro_hmm.md)
or
[`track_astro_iterative`](https://stratigraphy.github.io/WaverideR/reference/track_astro_iterative.md);
this function tests whether such a (minimal) tuning is supported by
astronomical structure beyond the power spectrum, and which rate range
(interpretation) is supported. Tests of an already tuned record are in
[`test_tuned`](https://stratigraphy.github.io/WaverideR/reference/test_tuned.md).

## Usage

``` r
astro_tune_test(
  data = NULL,
  sedmin = 2,
  sedmax = 25,
  numsed = 50,
  age_ma = NULL,
  targetE = NULL,
  targetP = NULL,
  targetO = NULL,
  weights = c(0.5, 1.5, 1, 2),
  sigma10 = 0.05,
  lambda = 0.3,
  step = 0.5,
  omega0 = 6,
  n_sur = 200,
  n_paths = 200,
  detrend = TRUE,
  run_multicore = FALSE,
  n_cores = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Evenly spaced records (matrix or data frame): depth in metres in
  column 1, and one or more proxies measured on the same depths in the
  further columns (no missing values). Every proxy is detrended (see
  `detrend`) and standardised; the column names are used as record
  names. Required, no default (`Default=NULL` gives an error).

- sedmin, sedmax:

  Lowest and highest sedimentation rate of the rate states of the hidden
  Markov model, in cm/kyr (numeric, \> 0, `sedmin` \< `sedmax`).
  Comparing the marginal log-likelihood of different ranges compares
  rate interpretations. `Default=2` and `Default=25`.

- numsed:

  Number of log-spaced rate states between `sedmin` and `sedmax`
  (integer, \>= 2). `Default=50`.

- age_ma:

  Age of the record in Ma (numeric) used to compute the target periods
  from AstroGeo22 (see
  [`astro_sedrate_test`](https://stratigraphy.github.io/WaverideR/reference/astro_sedrate_test.md));
  `NULL` uses present-day periods. `Default=NULL`.

- targetE, targetP, targetO:

  Your own eccentricity, precession and obliquity target periods in kyr
  (numeric vectors); the longest value of `targetE` is taken as the
  405-kyr term. They override the age-dependent values; `NULL` keeps
  those. `Default=NULL`.

- weights:

  Weights of the four term groups in the rate likelihood (numeric vector
  of length 4, \>= 0, in the order 405-kyr eccentricity, short
  eccentricity, obliquity, precession); each group weight is shared
  equally among its terms. `Default=c(0.5, 1.5, 1, 2)`.

- sigma10:

  Standard deviation of the change in log sedimentation rate per 10 m of
  depth (numeric, \> 0): the continuity prior of the rate path; smaller
  values give smoother rate paths. `Default=0.05`.

- lambda:

  Tempering factor applied to the rate log-likelihood, to account for
  the dependence of neighbouring depths (numeric, between 0 and 1; 1 =
  no tempering). `Default=0.3`.

- step:

  Depth step in metres between the states of the hidden Markov model
  (numeric, \> 0; at least one sample interval). `Default=0.5`.

- omega0:

  Central angular frequency of the Morlet wavelet used for the wavelet
  power (numeric, \> 0); larger values give a better frequency and a
  poorer depth resolution. `Default=6`.

- n_sur:

  Number of phase-randomised surrogates run through the whole chain to
  compute the p-values (integer, \> 0). `Default=200`.

- n_paths:

  Number of rate paths sampled from the posterior (forward-filtering
  backward-sampling) for the uncertainty of the age model and the rate
  (integer, \> 0). `Default=200`.

- detrend:

  Remove a linear trend from every record before the analysis (logical,
  TRUE or FALSE). `Default=TRUE`.

- run_multicore:

  Process the surrogates in parallel on several cores (logical, TRUE or
  FALSE). `Default=FALSE`.

- n_cores:

  Number of cores used for the surrogates (integer); `NULL` uses the
  option `waverider.n_cores` or, with `run_multicore = TRUE`, all but
  one detected core (at most `getOption("waverider.max_cores", 4)`).
  `Default=NULL`.

- verbose:

  Print progress messages (logical, TRUE or FALSE). `Default=FALSE`.

## Value

A list of class `astro_tune_test` with

- age_model:

  depth, time (kyr, median of the sampled paths), its 2.5 and 97.5%
  quantiles, and the time of the most probable (Viterbi) path, which is
  used for the tests

- rate:

  depth, most probable rate, posterior mean rate and its 2.5 and 97.5%
  quantiles

- tests:

  per record: r_prec, r_eam, r_efm, S_am, S_all and their probabilities

- summary:

  mean rate and duration (median of the sampled paths) with its 95%
  interval, the duration of the Viterbi path, marginal log-likelihood
  and its probability, mean S_am and S_all and their probabilities

- surrogates:

  the surrogate statistics

- settings:

  the settings and target periods

## Details

**1. Rate likelihood from all cycles.** The Morlet wavelet power of each
record is computed (in depth) at the frequencies that every trial
sedimentation rate predicts for the 405-kyr eccentricity, short
eccentricity, obliquity and precession terms. Power is divided by a
robust background (median power over depth, smoothed with a running
median over frequency), and the weighted sum of the log power ratios
over all terms gives the log-likelihood of every rate at every depth.
Combining all terms (rather than following one cycle) is the multi-order
approach of tacholess order tracking in machine diagnostics (Leclere et
al., 2016) and of multiband wavelet age modelling (Duesing et al.,
2021).

**2. Smooth rate path (hidden Markov model).** The log rate follows a
random walk with standard deviation `sigma10` per 10 m (continuity
prior, as the sedimentation-rate transition model of HMM-Match; Lin et
al., 2014). The Viterbi algorithm gives the most probable rate path, the
forward-backward algorithm the rate posterior at every depth, and
forward-filtering backward-sampling gives `n_paths` sampled rate paths,
which give the age-model uncertainty. `lambda` tempers the likelihood
(neighbouring depths are not independent). The marginal log-likelihood
of the model is reported for comparing rate ranges (e.g. two competing
interpretations) through `sedmin` and `sedmax`.

**3. Tuning.** Time is the integral of 1 / rate along the most probable
path. Because the path is smooth (a few degrees of freedom per 405-kyr
cycle), this is a minimal tuning: it does not impose precession-scale
frequency modulation (Zeeden et al., 2015).

**4. Modulation tests on the tuned records.** (a) precession amplitude
versus the record's own eccentricity band: wide precession filter (65%
bandwidth), Hilbert amplitude, low-pass (\< 1/75 kyr), Spearman
correlation (Zeeden et al., 2015); (b) short-eccentricity amplitude
versus the record's 405-kyr band (wide filter, low-pass \< 1/250 kyr);
(c) short-eccentricity frequency modulation (instantaneous frequency)
versus the 405-kyr band (Laurin et al., 2016). For an
eccentricity-forced record (a) and (b) have the sign of the proxy
response and (c) the opposite sign; `S_am = |a + b| / 2` and
`S_all = |a + b - c| / 3` use this sign consistency. Wide filters and
low-pass envelopes avoid the amplitude modulation that narrow filters
create from frequency modulation (Huybers and Aharonson, 2010).

**5. Significance.** Phase-randomised surrogates of the depth records
(Ebisuzaki, 1997; with the same random phases for all records, so the
correlation between records is kept; Prichard and Theiler, 1994) are run
through the whole chain (rate likelihood, path, tuning and tests). The
surrogates have exactly the power spectrum of the data, including the
power of any cycles, so the tests ask whether the data hold astronomical
structure beyond their power spectrum: a smooth rate path that aligns
power at all astronomical frequency ratios at the same depths
(`p_loglik`), and the modulation patterns above.

## References

Duesing, W. et al. (2021). Multiband wavelet age modeling for a ~293 m
(~600 kyr) sediment core from Chew Bahir Basin, southern Ethiopian Rift.
Frontiers in Earth Science 9, 594047 \<doi:10.3389/feart.2021.594047\>

Ebisuzaki, W. (1997). A method to estimate the statistical significance
of a correlation when the data are serially correlated. Journal of
Climate 10, 2147-2153
\<doi:10.1175/1520-0442(1997)010\<2147:AMTETS\>2.0.CO;2\>

Huybers, P. and Aharonson, O. (2010). Orbital tuning, eccentricity, and
the frequency modulation of climatic precession. Paleoceanography 25,
PA4228 \<doi:10.1029/2010PA001952\>

Laurin, J., Meyers, S.R., Galeotti, S. and Lanci, L. (2016). Frequency
modulation reveals the phasing of orbital eccentricity during Cretaceous
Oceanic Anoxic Event II and the Eocene hyperthermals. Earth and
Planetary Science Letters 442, 143-156
\<doi:10.1016/j.epsl.2016.02.047\>

Leclere, Q., Andre, H. and Antoni, J. (2016). A multi-order
probabilistic approach for instantaneous angular speed tracking:
debriefing of the CMMNO'14 diagnosis contest. Mechanical Systems and
Signal Processing 81, 375-386 \<doi:10.1016/j.ymssp.2016.02.053\>

Lin, L., Khider, D., Lisiecki, L.E. and Lawrence, C.E. (2014).
Probabilistic sequence alignment of stratigraphic records.
Paleoceanography 29, 976-989 \<doi:10.1002/2014PA002713\>

Prichard, D. and Theiler, J. (1994). Generating surrogate data for time
series with several simultaneously measured variables. Physical Review
Letters 73, 951-954 \<doi:10.1103/PhysRevLett.73.951\>

Zeeden, C., Meyers, S.R., Lourens, L.J. and Hilgen, F.J. (2015). Testing
astronomically tuned age models. Paleoceanography 30, 369-383
\<doi:10.1002/2014PA002762\>

## Examples

``` r
# \donttest{
set.seed(1)
res <- astro_tune_test(data = mag[1:600, ], sedmin = 0.3, sedmax = 2, numsed = 30,
    n_sur = 10, n_paths = 20)
res$summary
#>   mean_rate duration_kyr duration_viterbi duration_lo duration_hi     loglik
#> 1 0.4442067     6742.356         6742.356    6315.403    6742.356 -0.4605271
#>    p_loglik       S_am    p_S_am      S_all   p_S_all
#> 1 0.4545455 0.01822007 0.9090909 0.01169109 0.9090909
# }
```
