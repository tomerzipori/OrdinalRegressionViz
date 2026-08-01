# ROC curve panel for one pair of 2-level factors (Bayesian)

Internal workhorse of
[`bayesian_roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_roc_plot.md):
turns a grid of posterior category-probability rvars into a single ROC
panel.

## Usage

``` r
bayesian_roc_ggplot_2_vars(
  grid,
  var_signal,
  var_group,
  CI = 0.95,
  centrality = "mean",
  palette_thresholds = 7,
  palette_curves = "viridis",
  group_labels = NULL,
  empirical = NULL,
  ttl = ""
)
```

## Arguments

- grid:

  Data frame with the columns `var_signal`, `var_group` and one
  [posterior::rvar](https://mc-stan.org/posterior/reference/rvar.html)
  column per response category.

- var_signal:

  Name (string) of the 2-level factor playing the role of the SDT
  signal, e.g. old/new or true/fake.

- var_group:

  Name (string) of a 2-level factor covariate, e.g. pre/post. One panel
  is drawn per level.

- CI:

  Width of the credible intervals (between 0 and 1).

- centrality:

  Posterior point summary drawn as the curve: `"mean"` (default) or
  `"median"`.

- palette_thresholds:

  Integer index of a diverging palette passed to
  [`ggplot2::scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
  for the threshold points.

- palette_curves:

  Name of a viridis option (passed to
  [`ggplot2::scale_fill_viridis_d()`](https://ggplot2.tidyverse.org/reference/scale_viridis.html))
  used for the credible bands.

- group_labels:

  Optional character vector of length 2 with panel subtitles. Defaults
  to the title-cased levels of `var_group`.

- empirical:

  Optional data frame of observed hit/false-alarm rates (from
  `empirical_roc_points()`) to overlay.

- ttl:

  Plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
