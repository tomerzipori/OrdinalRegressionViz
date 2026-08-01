# Signal detection indices of an ordinal regression model

Computes the signal detection theory (SDT) quantities implied by a
cumulative regression model, per design cell: the sensitivity index
(`d_prime`), the scale ratio of the signal and noise distributions
(`sigma_ratio`), the response thresholds (`criterion ...`), and the area
under the implied ROC curve (`auc`).

## Usage

``` r
sdt_indices(
  model,
  var_signal,
  var_group = NULL,
  var_facet = NULL,
  ci = 0.95,
  ndraws = NULL,
  ...
)
```

## Arguments

- model:

  A cumulative regression model:
  [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html),
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html), or a
  `brmsfit` from
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- var_signal:

  Name (string) of the 2-level factor playing the role of the SDT
  signal. Its first level is the noise distribution.

- var_group, var_facet:

  Optional names (strings) of additional 2-level factor covariates;
  indices are computed for every combination of their levels.

- ci:

  Width of the credible intervals (`brmsfit` only).

- ndraws:

  Optional number of posterior draws to use (`brmsfit` only); `NULL`
  uses all draws.

- ...:

  Passed on to methods.

## Value

A data frame with one row per index and design cell: the levels of
`var_group`/`var_facet` (when supplied), `index`, `estimate`, `lower`,
and `upper`.

## Details

Definitions follow the classic SDT conventions, on the latent scale of
the model's link function (normal for probit):

- `d_prime` is the location difference between the signal (second level
  of `var_signal`) and noise (first level) distributions, in units of
  the noise distribution's scale.

- `sigma_ratio` is the signal scale divided by the noise scale (1 unless
  the model includes a `disc` part;
  [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html) scale
  models are not supported).

- `criterion <k|k+1>` are the response thresholds of the cell, in raw
  latent units (matching
  [`SDT_distributions_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/SDT_distributions_plot.md)).

- `auc` is the probability that a random signal draw exceeds a random
  noise draw. It has a closed form under probit and is computed by
  numeric integration for other links.

For `brmsfit` models the full posterior of every index is used, and the
returned `estimate` is the posterior median with a `ci` central credible
interval. For `clm`/`clmm` models the indices are point estimates
(`lower`/`upper` are `NA`); use the Bayesian interface if you need
uncertainty for the derived quantities.

## See also

[`SDT_distributions_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/SDT_distributions_plot.md)
and
[`bayesian_SDT_distribution_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/bayesian_SDT_distribution_plot.md)
for the corresponding latent-distribution displays.

## Examples

``` r
fit <- ordinal::clm(value ~ target * time, data = sdt_ratings, link = "probit")
sdt_indices(fit, var_signal = "target", var_group = "time")
#>    time         index   estimate lower upper
#> 1   pre       d_prime  1.1929128    NA    NA
#> 2   pre   sigma_ratio  1.0000000    NA    NA
#> 3   pre           auc  0.8005302    NA    NA
#> 4   pre criterion 1|2 -0.9950796    NA    NA
#> 5   pre criterion 2|3 -0.2661625    NA    NA
#> 6   pre criterion 3|4  0.3947510    NA    NA
#> 7   pre criterion 4|5  1.0447956    NA    NA
#> 8   pre criterion 5|6  1.6681026    NA    NA
#> 9  post       d_prime  1.3120230    NA    NA
#> 10 post   sigma_ratio  1.0000000    NA    NA
#> 11 post           auc  0.8232289    NA    NA
#> 12 post criterion 1|2 -0.9050570    NA    NA
#> 13 post criterion 2|3 -0.1761400    NA    NA
#> 14 post criterion 3|4  0.4847736    NA    NA
#> 15 post criterion 4|5  1.1348182    NA    NA
#> 16 post criterion 5|6  1.7581252    NA    NA
```
