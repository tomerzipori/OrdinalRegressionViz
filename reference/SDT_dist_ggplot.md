# Single latent-distribution panel (frequentist)

Internal workhorse of
[`SDT_distributions_plot()`](https://tomerzipori.github.io/OrdinalRegressionViz/reference/SDT_distributions_plot.md):
draws the reference and signal latent distributions of one design cell
together with its response thresholds.

## Usage

``` r
SDT_dist_ggplot(
  ref_mean = 0,
  group_mean = 1,
  thresholds,
  dens = link_density("probit"),
  signal_labels = c("Noise", "Signal"),
  palette = 9,
  alpha = 0.7,
  plot_limits = c(-3, 5),
  ttl = "",
  x_label = ""
)
```

## Arguments

- ref_mean:

  Location of the reference ("noise") distribution.

- group_mean:

  Location of the signal distribution.

- thresholds:

  Named numeric vector of response thresholds.

- dens:

  Density function of the latent distribution, see `link_density()`.

- signal_labels:

  Character vector of length 2 with the legend labels (reference first).

- palette:

  Integer index of a diverging palette passed to
  [`ggplot2::scale_color_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
  for the threshold lines.

- alpha:

  Opacity of the threshold lines.

- plot_limits:

  Numeric vector of length 2: x-axis limits of the latent scale.

- ttl:

  Plot title.

- x_label:

  Label of the x (latent) axis.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
