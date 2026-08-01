# Forest plot for a Bayesian random-effects meta-analysis

Draws a forest plot — posterior densities and intervals per study, plus
the pooled effect — for a Bayesian random-effects meta-analysis fitted
with [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
The model is expected to contain nested grouping terms of the form
`(1 | author/study)` (defaults: `Author`, `Study`); densities are
colored by the higher-level (author) grouping.

## Usage

``` r
bayesian_forest(
  b_meta_analysis,
  author_var = "Author",
  study_var = "Study",
  palette = 11,
  effect_label = "Standardized Mean Difference",
  ttl = "",
  hjust_ttl = 1
)
```

## Arguments

- b_meta_analysis:

  A `brmsfit` random-effects meta-analysis with nested author/study
  grouping terms.

- author_var:

  Name (string) of the higher-level grouping variable.

- study_var:

  Name (string) of the study-level grouping variable nested in
  `author_var`.

- palette:

  Integer: number of colors drawn from the Spectral palette (see
  [`RColorBrewer::brewer.pal()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html))
  before interpolation.

- effect_label:

  Label of the x axis (the summary measure).

- ttl:

  Plot title.

- hjust_ttl:

  Horizontal adjustment of the title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

The layout is adapted from Harrer, Cuijpers, Furukawa, & Ebert (2021),
*Doing Meta-Analysis with R: A Hands-On Guide*
(<https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/>).

## Examples

``` r
if (FALSE) { # \dontrun{
# b_meta <- brms::brm(
#   yi | se(sei) ~ 1 + (1 | Author / Study),
#   data = meta_data
# )
bayesian_forest(b_meta)
} # }
```
