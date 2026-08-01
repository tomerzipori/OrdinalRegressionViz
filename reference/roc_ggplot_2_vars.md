# ROC curve panel for one pair of 2-level factors (frequentist)

Internal workhorse of
[`roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/roc_plot.md):
turns emmeans cumulative- and exceedance-probability grids into a single
ROC panel.

## Usage

``` r
roc_ggplot_2_vars(
  ems_sensitivity,
  ems_specificity,
  var_group,
  CI = 0.95,
  palette_groups = 2,
  palette_thresholds = 2,
  group_labels = NULL,
  empirical = NULL,
  ttl = ""
)
```

## Arguments

- ems_sensitivity:

  `emmGrid` of cumulative probabilities (`mode = "cum.prob"`) at the
  first (signal) level of the SDT variable.

- ems_specificity:

  `emmGrid` of exceedance probabilities (`mode = "exc.prob"`) at the
  second (noise) level of the SDT variable.

- var_group:

  Name (string) of a 2-level factor covariate; one ROC curve is drawn
  per level.

- CI:

  Confidence level of the intervals (between 0 and 1).

- palette_groups:

  Integer index of a qualitative palette passed to
  [`ggplot2::scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
  for the confidence bands.

- palette_thresholds:

  Integer index of a diverging palette passed to
  [`ggplot2::scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
  for the threshold points.

- group_labels:

  Optional character vector of length 2 with legend labels for the
  levels of `var_group`. Defaults to its title-cased levels.

- empirical:

  Optional data frame of observed hit/false-alarm rates (from
  `empirical_roc_points()`) to overlay.

- ttl:

  Plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
