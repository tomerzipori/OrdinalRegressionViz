% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bayesian_distributions_plot_2_vars}
\alias{bayesian_distributions_plot_2_vars}
\title{Plot of latent normal distributions for two 2-level categorical variables Bayesian ordinal probit regression models}
\usage{
bayesian_distributions_plot_2_vars(
  b_model,
  var_signal = "target",
  var_group = "time",
  response = "value",
  group1_ttl = "",
  group2_ttl = "",
  response_scale = c(1:7),
  palette = 9,
  plot_range = c(-11, 11),
  ttl = ""
)
}
\arguments{
\item{b_model}{Bayesian ordinal probit regression model - a brms object}

\item{var_signal}{variable analogous to the classic SDT signal, e.g. old/new}

\item{var_group}{variable name of the co-variate, e.g. experimental group}

\item{response}{variable name of the response}

\item{response_scale}{integer vector of response scale levels, e.g 1 to 7}

\item{palette}{integer representing a divergent palette in the scale_color_brewer function}

\item{plot_range}{2-integer vector representing x-axis limits of the plot}

\item{ttl}{plot's title}
}
\value{
a ggplot2::ggplot plot object
}
\description{
This function makes a latent distributions plot for two 2-level categorical variables Bayesian ordinal regression models.
}
\keyword{internal}
