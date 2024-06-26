% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian_forest.R
\name{bayesian_forest}
\alias{bayesian_forest}
\title{Forest plot for Bayesian random effects meta-analysis}
\usage{
bayesian_forest(
  b_meta_analysis,
  palette = 11,
  ttl = "",
  hjust_ttl = 1,
  filename = NULL,
  path = getwd(),
  width = 2450,
  height = 1446
)
}
\arguments{
\item{b_meta_analysis}{Bayesian random effects meta-analysis - a brms object}

\item{palette}{integer representing a palette in the scale_color_brewer function}

\item{ttl}{plot's title}

\item{hjust_ttl}{constant offset of plot title to keep in in the center}

\item{filename}{if not NULL, the name of the png file to save the plot as}

\item{path}{folder path to save the plot in}

\item{width}{width of the saved plot - in pixels}

\item{height}{height of the saved plot - in pixels}
}
\value{
a ggplot plot object
}
\description{
This function makes a forest plot for Bayesian random effects meta-analyses. Code is heavily inspired by Harrer et al. online guide - https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/.
Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021). Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London: Chapman & Hall/CRC Press. ISBN 978-0-367-61007-4.
}
\details{
Important - this function requires a 'Study' variable for individual studies and an 'Author' clustering variable of studies.
Important - In the final plot studies are colored by their author.
}
