# ROC curve of a Bayesian ordinal regression model

Draws a receiver operating characteristic (ROC) curve — hit rate against
false-alarm rate across the response thresholds — from a Bayesian
cumulative ("ordinal") regression model fitted with
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html). One
curve is drawn per level of `var_group`, with curvewise credible bands;
when `var_facet` is given, side-by-side ROC panels are drawn for its two
levels.

## Usage

``` r
bayesian_roc_plot(
  b_model,
  var_signal,
  var_group,
  var_facet = NULL,
  CI = 0.95,
  centrality = c("mean", "median"),
  palette_thresholds = 7,
  palette_curves = "viridis",
  group_labels = NULL,
  facet_labels = NULL,
  show_empirical = FALSE,
  ttl = ""
)
```

## Arguments

- b_model:

  A `brmsfit` cumulative regression model.

- var_signal:

  Name (string) of the 2-level factor playing the role of the SDT
  signal, e.g. old/new or true/fake.

- var_group:

  Name (string) of a 2-level factor covariate, e.g. pre/post. One panel
  is drawn per level.

- var_facet:

  Optional name (string) of a second 2-level factor covariate; panels
  are duplicated for each of its levels.

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

## Details

Population-level ("fixed effects only") predicted probabilities are used
(`re_formula = NA` in
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)).

## See also

[`roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/roc_plot.md)
for [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) models.

## Examples

``` r
if (FALSE) { # \dontrun{
bayesian_roc_plot(
  b_fit,
  var_signal = "target", var_group = "time"
)
} # }
```
