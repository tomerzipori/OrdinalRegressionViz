#' Plot of latent normal distributions for two or three 2-level categorical variables Bayesian ordinal probit regression models
#'
#' This function makes a latent distributions plot for two or three 2-level categorical variables Bayesian ordinal regression models.
#' @param b_model Bayesian ordinal probit regression model - a brms object
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param var_facet variable to facet by
#' @param response variable name of the response
#' @param response_scale integer vector of response scale levels, e.g 1 to 7
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param plot_range 2-integer vector representing x-axis limits of the plot
#' @param ttl plot's title
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
bayesian_SDT_distribution_plot <- function(b_model,
                                           var_signal = "target",
                                           var_group = "time",
                                           var_facet = NULL,
                                           response = "value",
                                           response_scale = c(1:7),
                                           palette = 9,
                                           alpha = 0.6,
                                           plot_range = c(-11, 11),
                                           ttl = "",
                                           filename = NULL,
                                           path = getwd(),
                                           width = 2450,
                                           height = 1446) {

  if (is.null(var_facet)) {
    out_plot <- bayesian_distributions_plot_2_vars(b_model = b_model,
                                                   var_signal = var_signal,
                                                   var_group = var_group,
                                                   response = response,
                                                   response_scale = response_scale,
                                                   palette = palette,
                                                   alpha = alpha,
                                                   plot_range = plot_range,
                                                   ttl = ttl)
  }

  if (!is.null(var_facet)) {
    out_plot <- bayesian_distributions_plot_3_vars(b_model = b_model,
                                                   var_signal = var_signal,
                                                   var_group = var_group,
                                                   var_facet = var_facet,
                                                   response = response,
                                                   response_scale = response_scale,
                                                   palette = palette,
                                                   alpha = alpha,
                                                   plot_range = plot_range,
                                                   ttl = ttl)
  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

