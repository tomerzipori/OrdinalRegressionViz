# ROC curve of a frequentist ordinal regression model

Draws a receiver operating characteristic (ROC) curve — hit rate against
false-alarm rate across the response thresholds — from a cumulative link
model fitted with
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) or
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html). One
curve is drawn per level of `var_group`, with confidence bands computed
from
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
cumulative-probability estimates; when `var_facet` is given,
side-by-side ROC panels are drawn for its two levels.

## Usage

``` r
roc_plot(
  model,
  var_signal,
  var_group,
  var_facet = NULL,
  CI = 0.95,
  palette_groups = 2,
  palette_thresholds = 2,
  group_labels = NULL,
  facet_labels = NULL,
  show_empirical = FALSE,
  ttl = ""
)
```

## Arguments

- model:

  A `clm` or `clmm` ordinal regression model.

- var_signal:

  Name (string) of the 2-level factor playing the role of the SDT
  signal, e.g. old/new or true/fake.

- var_group:

  Name (string) of a 2-level factor covariate; one ROC curve is drawn
  per level.

- var_facet:

  Optional name (string) of a second 2-level factor covariate; one panel
  is drawn per level.

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

- facet_labels:

  Optional character vector of length 2 with panel titles for the levels
  of `var_facet`. Defaults to its title-cased levels.

- show_empirical:

  If `TRUE`, the observed (empirical) hit and false-alarm rates are
  overlaid as crosses — a quick visual check of how well the
  model-implied ROC matches the data.

- ttl:

  Plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object (or a patchwork of two panels when `var_facet` is given).

## See also

[`bayesian_roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_roc_plot.md)
for [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
models, and
[`SDT_distributions_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/SDT_distributions_plot.md)
for the corresponding latent distributions.

## Examples

``` r
fit <- ordinal::clm(value ~ target * time, data = sdt_ratings)
roc_plot(fit, var_signal = "target", var_group = "time")
```
