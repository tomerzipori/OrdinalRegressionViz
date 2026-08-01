# Grouped posterior predictive check for a Bayesian ordinal model

Draws a grouped bar plot comparing the observed response distribution
with posterior predictive draws (via
[`bayesplot::ppc_bars_grouped()`](https://mc-stan.org/bayesplot/reference/PPC-discrete.html)),
with one facet per combination of the grouping variables.

## Usage

``` r
ordinal_model_ppd_check(b_model, group_vars, ndraws = 40, ttl = "")
```

## Arguments

- b_model:

  A `brmsfit` cumulative regression model.

- group_vars:

  Character vector with the names of the model variables whose
  combinations define the facets (e.g. `c("target", "time")`).

- ndraws:

  Number of posterior predictive draws to use.

- ttl:

  Plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
ordinal_model_ppd_check(b_fit, group_vars = c("target", "time"))
} # }
```
