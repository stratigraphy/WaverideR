# Discriticizes lithologs

Discriticizes lithologs to allow further time-series analysis first the
Greatest common divisor/highest common factor is calculated which is
then used to discriticize the litholog to an evenly sampled data series.
The function is designed to place the boundary at the original depth
level of the bed boundaries. The Greatest common divisor/highest common
factor can be a very small number as such the discriticized data set can
be large which impacts computational performance later on therefore a
linear interpolation option is added to downscale the data to allow for
computational efficiency later on. This is made to discriticize
lithologs created using the 'StratigrapheR' package. as such the same
data format for input is used. eg. column 1 is bottom of the bed, column
2 is top of bed, column is depth rank/proxy value

## Usage

``` r
lithlog_disc(
  litholog = NULL,
  subset_fact = 10,
  lin_interp = FALSE,
  dt = NULL,
  genplot = FALSE,
  x_lab = "rank",
  y_lab = "depth (m)",
  keep_editable = FALSE
)
```

## Arguments

- litholog:

  Litholog as a three-column matrix or data frame (the bed format of
  'StratigrapheR'), one row per bed ordered from the base of the section
  upwards: column 1 the base (bottom) depth of the bed, column 2 the top
  depth of the bed (larger than column 1, e.g. in metres) and column 3
  the depth rank or proxy value of the bed. Required (`Default=NULL`).

- subset_fact:

  Divisor of the greatest common divisor of the bed thicknesses
  (numeric, \> 0): the litholog is discretised with a step of (greatest
  common divisor) / `subset_fact`, so larger values give a finer, larger
  series. `Default=10`.

- lin_interp:

  Linearly interpolate the discretised series to the regular step `dt`,
  e.g. to reduce a very finely sampled series (logical, TRUE or FALSE).
  `Default=FALSE`.

- dt:

  Sampling interval of the interpolated series, in the depth unit of
  `litholog` (numeric, \> 0); required when `lin_interp = TRUE`,
  otherwise ignored. `Default=NULL`.

- genplot:

  Plot the discretised series, value against depth (logical, TRUE or
  FALSE). `Default=FALSE`.

- x_lab:

  Label of the x-axis (rank or proxy value) of the plot (character).
  `Default="rank"`.

- y_lab:

  Label of the y-axis (depth) of the plot (character).
  `Default="depth (m)"`.

- keep_editable:

  Keep the graphical parameters (`par`) of the plot after the function
  returns, so that features can be added to it (logical, TRUE or FALSE);
  `FALSE` restores the previous settings. `Default=FALSE`.

## Value

Returns a matrix with 2 columns, the first column is depth the second
columns is the depth/rank proxy If genplot is `Default=TRUE` then a plot
of the discriticizes time series is plotted

## References

Wouters, S., Da Silva, A.-C., Boulvain, F., and Devleeschouwer, X..
2021. StratigrapheR: Concepts for Litholog Generation in R. The R
Journal. \<doi:10.32614/RJ-2021-039\>

## Examples

``` r
# Convert depth rank record to a discrete proxy record to allow for further
# analysis in which discrete time series are needed
depth_rank_example_disc <- lithlog_disc(litholog = depth_rank_example,
           subset_fact = 10,
           genplot = FALSE,
           x_lab = "rank",
           y_lab = "depth (m)",
           keep_editable=FALSE)

```
