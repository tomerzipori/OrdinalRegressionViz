% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bayesian_roc_ggplot_2_vars}
\alias{bayesian_roc_ggplot_2_vars}
\title{ROC Curve for two 2-level categorical predictor Bayesian ordinal probit regression models}
\usage{
bayesian_roc_ggplot_2_vars(
  grid,
  var_signal = "target",
  var_group = "time",
  response = "value",
  CI = 0.95,
  centrality = "mean",
  palette = 7,
  ttl = ""
)
}
\arguments{
\item{grid}{data frame of conditional probabilities posteriors}

\item{var_signal}{variable analogous to the classic SDT signal, e.g. old/new}

\item{var_group}{variable name of the co-variate, e.g. experimental group}

\item{response}{variable name of the response}

\item{CI}{credible interval width between 0 and 1}

\item{centrality}{centrality measure for posterior, either "mean" or "median"}

\item{palette}{integer representing a divergent palette in the scale_color_brewer function}

\item{ttl}{plot's title}
}
\value{
a ggplot2::ggplot plot object
}
\description{
This function makes a ROC curve for two categorical variables Bayesian ordinal probit regression models.
}
\keyword{internal}
