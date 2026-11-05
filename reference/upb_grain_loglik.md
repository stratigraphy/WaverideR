# Eruption-age likelihoods from single-zircon U-Pb dates (for astroMCBchron)

`upb_grain_loglik()` turns single-zircon dates of one or more ash beds
into one eruption-age log-likelihood function per bed, for the
`date_loglik` argument of
[`astroMCBchron`](https://stratigraphy.github.io/WaverideR/reference/astroMCBchron.md).
Instead of a weighted mean, every grain is modelled as a draw from a
mixture of

- a crystallisation interval: the grain crystallised uniformly between
  the eruption age *te* and *te* + *D* (antecrystic growth before
  eruption), with the duration *D* marginalised over a uniform prior
  between 0 and `D_max`;

- inherited or recycled grains (probability `p_inherit`): uniform over
  `W_inherit` before the eruption;

- grains with residual Pb loss (probability `p_pbloss`): uniform over
  `W_pbloss` after the eruption,

each convolved with the analytical uncertainty of the grain. The
likelihood of *te* is computed on a grid, so no grain has to be rejected
beforehand. With `D_max` = 0 and no mixture components the result is the
likelihood of a weighted mean; with a large `D_max` the eruption age
approaches the youngest concordant grains (as in the Bayesian
eruption-age estimator of Keller et al., 2018). Which regime applies
depends on the ratio of the crystallisation interval to the analytical
uncertainty of the grains, which is why the choice is left to the data
rather than fixed beforehand.

The functions return a log-likelihood (not a posterior): no prior is
placed on *te*, so the stratigraphic order and the astrochronology in
`astroMCBchron` act as the prior. `astroMCBchron` normalises each
function numerically and combines it with its outlier mixture.

## Usage

``` r
upb_grain_loglik(
  grains,
  D_max = 0.6,
  p_inherit = 0.1,
  p_pbloss = 0.05,
  W_inherit = 20,
  W_pbloss = 3,
  step = 0.002,
  nD = 61
)
```

## Arguments

- grains:

  Single-grain dates (data frame) with the columns `id` (name of the
  bed, character or factor), `age` (single-grain date, numeric, in Ma or
  any other age unit used consistently with `D_max`, `W_inherit`,
  `W_pbloss` and `step`) and either `sd` (analytical 1 sigma
  uncertainty) or `sd2` (analytical 2 sigma uncertainty, halved). Rows
  with missing ages or uncertainties are skipped. Required, no default.

- D_max:

  Maximum duration of the crystallisation interval before the eruption
  (numeric, \>= 0, in age units); the duration *D* is marginalised over
  a uniform prior between 0 and `D_max`. 0 gives the likelihood of a
  weighted mean, large values move the eruption age towards the youngest
  grains. `Default=0.6` (Ma).

- p_inherit:

  Probability that a grain is inherited or recycled (numeric, \>= 0;
  `p_inherit + p_pbloss` must be \< 1); such grains are modelled as
  uniform over `W_inherit` older than the eruption. `Default=0.1`.

- p_pbloss:

  Probability that a grain suffered residual Pb loss (numeric, \>= 0);
  such grains are modelled as uniform over `W_pbloss` younger than the
  eruption. `Default=0.05`.

- W_inherit:

  Width of the inherited-grain component, i.e. how much older than the
  eruption inherited grains can be (numeric, \> 0, in age units).
  `Default=20` (Ma).

- W_pbloss:

  Width of the Pb-loss component, i.e. how much younger than the
  eruption grains with Pb loss can appear (numeric, \> 0, in age units).
  `Default=3` (Ma).

- step:

  Grid step of the eruption age on which the likelihood is computed
  (numeric, \> 0, in age units); smaller steps are more precise and
  slower. `Default=0.002` (2 kyr for ages in Ma).

- nD:

  Number of grid values of the crystallisation duration *D* between 0
  and `D_max` used for the marginalisation (integer, \>= 1).
  `Default=61`.

## Value

A list with `loglik` (named list of functions, one per bed, of the
eruption age in age units) and `summary` (data.frame per bed: number of
grains, maximum- likelihood eruption age, mean and 2 sd of the
normalised likelihood, youngest grain and weighted mean with MSWD for
comparison).

## References

Keller, C.B., Schoene, B., Samperton, K.M., 2018. A stochastic sampling
approach to zircon eruption age interpretation. Geochemical Perspectives
Letters 8, 31-35.

## Examples

``` r
set.seed(2)
g <- data.frame(id = rep(c("ashA", "ashB"), each = 8),
                age = c(430 + runif(8, 0, 0.1) + rnorm(8, 0, 0.04),
                        431 + runif(8, 0, 0.1) + rnorm(8, 0, 0.04)),
                sd = 0.04)
g$age[3] <- 433.2                                     # one inherited grain
gl <- upb_grain_loglik(g)
gl$summary
#>     id n    te_ml  te_mean    te_2sd youngest youngest_2sd   wm_all wm_all_2sd
#> 1 ashA 8 430.0012 429.9858 0.1104326 430.0072         0.08 430.4641 0.02828427
#> 2 ashB 8 431.0032 430.9885 0.1037042 430.9652         0.08 431.0666 0.02828427
#>     mswd_all
#> 1 765.789909
#> 2   2.186294
```
