# Latent SDT distributions of a Bayesian ordinal regression model

Plots the latent ("perceived signal") distributions implied by a
Bayesian cumulative regression model, together with the posterior
distributions of the response thresholds (criteria), in the style of
signal detection theory (SDT). One panel is drawn for every level of
`var_group`, optionally split by a third variable (`var_facet`).

## Usage

``` r
bayesian_SDT_distribution_plot(
  b_model,
  var_signal,
  var_group,
  var_facet = NULL,
  signal_labels = NULL,
  group_labels = NULL,
  facet_labels = NULL,
  palette = 9,
  alpha = 0.6,
  plot_range = c(-11, 11),
  ci = 0.9,
  ndraws = NULL,
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

- signal_labels:

  Optional character vector of length 2 with legend labels for the two
  latent distributions (reference level first). Defaults to the
  title-cased levels of `var_signal`.

- group_labels:

  Optional character vector of length 2 with panel subtitles. Defaults
  to the title-cased levels of `var_group`.

- facet_labels:

  Optional character vector of length 2 with panel titles for the levels
  of `var_facet`. Defaults to its title-cased levels.

- palette:

  Integer index of a sequential palette passed to
  [`ggplot2::scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
  for the threshold slabs.

- alpha:

  Opacity of the threshold slabs.

- plot_range:

  Numeric vector of length 2: x-axis limits of the latent scale.

- ci:

  Width of the curvewise credible band around each density curve
  (between 0 and 1).

- ndraws:

  Optional number of posterior draws used to compute the density bands.
  The default (`NULL`) uses all draws; a value such as `500` makes the
  plot considerably faster with little visual change.

- ttl:

  Plot title.

## Value

A
[patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
object combining the panels; it can be printed, modified with `&`, or
saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The model is expected to be a
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) fit
with a cumulative family. The latent distribution follows from the
model's link function: normal for `probit` (the classic SDT display),
logistic for `logit`, Cauchy for `cauchit`, and Gumbel for `cloglog`.
When the discrimination parameter `disc` is modeled (unequal-variance
SDT), the posterior of each distribution's scale is used; otherwise all
latent distributions have unit scale. Models without some interaction
terms are supported: terms the model does not contain are simply treated
as zero.

## See also

[`SDT_distributions_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/SDT_distributions_plot.md)
for [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) models,
and
[`bayesian_roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_roc_plot.md)
for the corresponding ROC curves.

## Examples

``` r
if (FALSE) { # \dontrun{
# b_fit <- brms::brm(
#   brms::bf(value ~ target * time, disc ~ target * time),
#   family = brms::cumulative("probit"),
#   data = sdt_ratings
# )
bayesian_SDT_distribution_plot(
  b_fit,
  var_signal = "target", var_group = "time",
  plot_range = c(-4, 6), ndraws = 500
)
} # }
```
