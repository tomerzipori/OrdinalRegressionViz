% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_SDT_distribution_plot.R
\name{bayesian_SDT_distribution_plot}
\alias{bayesian_SDT_distribution_plot}
\title{Plot of latent normal distributions for two or three 2-level categorical variables Bayesian ordinal probit regression models}
\usage{
bayesian_SDT_distribution_plot(
  b_model,
  var_signal = "target",
  var_group = "time",
  var_facet = NULL,
  response = "value",
  response_scale = c(1:7),
  palette = 9,
  plot_range = c(-11, 11),
  ttl = "",
  filename = NULL,
  path = getwd(),
  width = 2450,
  height = 1446
)
}
\arguments{
\item{b_model}{Bayesian ordinal probit regression model - a brms object}

\item{var_signal}{variable analogous to the classic SDT signal, e.g. old/new}

\item{var_group}{variable name of the co-variate, e.g. experimental group}

\item{var_facet}{variable to facet by}

\item{response}{variable name of the response}

\item{response_scale}{integer vector of response scale levels, e.g 1 to 7}

\item{palette}{integer representing a divergent palette in the scale_color_brewer function}

\item{plot_range}{2-integer vector representing x-axis limits of the plot}

\item{ttl}{plot's title}

\item{filename}{if not NULL, the name of the png file to save the plot as}

\item{path}{folder path to save the plot in}

\item{width}{width of the saved plot - in pixels}

\item{height}{height of the saved plot - in pixels}
}
\value{
a ggplot plot object
}
\description{
This function makes a latent distributions plot for two or three 2-level categorical variables Bayesian ordinal regression models.
}
