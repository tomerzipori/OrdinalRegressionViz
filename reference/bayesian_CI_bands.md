# Credible/confidence band polygon for ROC curves

Builds the polygon that connects the lower and upper interval bounds of
a ROC curve so it can be drawn with
[`ggplot2::geom_polygon()`](https://ggplot2.tidyverse.org/reference/geom_polygon.html).

## Usage

``` r
bayesian_CI_bands(data, x, xmin, xmax, y, ymin, ymax, group = NULL)
```

## Arguments

- data:

  Data frame holding the curve and its interval bounds.

- x, xmin, xmax:

  Names (strings) of the x coordinate and its bounds.

- y, ymin, ymax:

  Names (strings) of the y coordinate and its bounds.

- group:

  Optional character vector of grouping column names.

## Value

A data frame of polygon vertices with columns `x` and `y` (plus any
grouping columns).
