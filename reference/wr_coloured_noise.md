# Colour noise and palette helpers

`wr_coloured_noise` generates AR(1) red noise (replaces
`colorednoise::colored_noise`); `wr_brewer_palettes` lists the available
ColorBrewer palette names and their maximum number of colours.

## Usage

``` r
wr_coloured_noise(timesteps, mean, sd, phi)

wr_brewer_palettes()

wr_colored_noise(timesteps, mean, sd, phi)
```

## Arguments

- timesteps:

  Number of values to generate (integer, \>= 1; with 1 a single normal
  value is returned). Required, no default.

- mean:

  Mean of the noise (numeric, in the unit of the series to be
  simulated). Required, no default.

- sd:

  Standard deviation of the noise (numeric, \> 0); the marginal standard
  deviation of the AR(1) process. Required, no default.

- phi:

  Lag-1 autocorrelation of the AR(1) process (numeric, between 0 and 1,
  below 1; 0 gives white noise, values close to 1 give strongly red
  noise). Required, no default.

## Value

`wr_coloured_noise`: numeric vector. `wr_brewer_palettes`: data frame
with the palettes.

## Examples

``` r
set.seed(1)
x <- wr_coloured_noise(100, mean = 0, sd = 1, phi = 0.8)
head(wr_brewer_palettes())
#>      maxcolors
#> BrBG        11
#> PiYG        11
#> PRGn        11
#> PuOr        11
#> RdBu        11
#> RdGy        11
```
