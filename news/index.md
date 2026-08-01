# Changelog

## OrdinalRegressionViz 0.2.0

### Breaking changes

- `var_signal` and `var_group` are now required arguments in all plot
  functions. They previously defaulted to the study-specific variable
  names `"target"` and `"time"`.
- The `filename`, `path`, `width`, and `height` arguments were removed
  from all plot functions. The functions now only return plot objects;
  save them with
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
- The `response` / `response_scale` arguments were removed; the response
  variable and its scale are now derived from the model.
- [`ordinal_model_ppd_check()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/ordinal_model_ppd_check.md):
  the `vars` argument (a count of columns) was replaced by `group_vars`,
  a character vector of model variable names. Positional assumptions
  about the columns of the model data were removed.
- [`bayesian_SDT_distribution_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_SDT_distribution_plot.md):
  the `plot_range` and label arguments were reorganized. Legend labels
  now derive from the levels of `var_signal` instead of the hardcoded
  `"True tweets"` / `"Fake tweets"`, and threshold labels derive from
  the model’s response scale.
- brms and ordinal moved from Imports to Suggests; their availability is
  checked at runtime.
- The joke startup message was removed.

### New features

- New
  [`sdt_indices()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/sdt_indices.md):
  tidy signal detection indices (d’, sigma ratio, criteria, AUC) per
  design cell, for `clm`/`clmm` (point estimates) and `brmsfit` models
  (posterior medians with credible intervals).
- New `show_empirical` argument in
  [`roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/roc_plot.md)
  and
  [`bayesian_roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_roc_plot.md):
  overlays the observed hit and false-alarm rates on the model-implied
  ROC as a visual goodness-of-fit check.
- All latent-distribution and ROC functions accept any 2-level factors,
  and work with models with or without interaction terms (terms absent
  from the model are treated as zero).
- Link functions beyond probit are supported: logit (logistic), cauchit
  (Cauchy), and cloglog/loglog (Gumbel) latent distributions.
- Equal-variance Bayesian models (no `disc` part) are now supported by
  [`bayesian_SDT_distribution_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_SDT_distribution_plot.md).
- New arguments: `signal_labels`, `group_labels`, and `facet_labels` in
  all plot functions; `CI` in
  [`roc_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/roc_plot.md);
  `ci` and `ndraws` (for speed) in
  [`bayesian_SDT_distribution_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_SDT_distribution_plot.md);
  `author_var`, `study_var`, and `effect_label` in
  [`bayesian_forest()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_forest.md).
- New packaged dataset `sdt_ratings`: simulated perceived-truth ratings
  in a 2x2x2 signal-detection design, plus a getting-started vignette
  ([`vignette("visualizing-ordinal-models")`](https://tomerzipori.github.io/OrdinalRegressionViz/articles/visualizing-ordinal-models.md)).

### Bug fixes

- ROC curves of `clm` models now correctly distinguish groups by
  linetype; previously a constant string was mapped to the linetype
  aesthetic.
- [`rvar_cumsum()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/rvar_cumsum.md)
  no longer corrupts length-1 inputs.
- Fixed hardcoded publication-specific coefficient names (`b_targetfake`
  and similar) that failed on any other dataset.
- Threshold labels no longer assume a 13-point (or 7-point) response
  scale.
- Removed [`library()`](https://rdrr.io/r/base/library.html) calls
  inside package functions and undeclared dependencies (tidyverse,
  ggpubr, parameters, glue, insight, purrr).

### Infrastructure

- Added a testthat suite.
- Added GitHub Actions workflows for R CMD check and pkgdown.
- Completed roxygen2 documentation for all exported and internal
  functions.
