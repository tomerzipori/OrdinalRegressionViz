#' ROC Curve for two or three 2-level categorical predictor probit ordinal regression models
#'
#' This function makes a ROC curve for two or three categorical variables ordinal regression models.
#' @param model ordinal regression model, as returned by clm or clmm functions
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param var_facet variable to facet by
#' @param response variable name of the response
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param ttl plot's title
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
roc_plot <- function(model,
                     var_signal = "target",
                     var_group = "time",
                     var_facet = NULL,
                     response = "value",
                     palette = 2,
                     ttl = "",
                     filename = NULL,
                     path = getwd(),
                     width = 2450,
                     height = 1446) {

  signal_levels <- levels(model$model[var_signal][1,])

  if (is.null(var_facet)) {

    f <- formula(paste0("~ cut | ", var_group, " + ", var_signal))
    ems_sens <- emmeans::emmeans(model, f, at = list(target = signal_levels[1]), mode = "cum.prob")
    ems_spec <- emmeans::emmeans(model, f, at = list(target = signal_levels[2]), mode = "exc.prob")
    out_plot <- roc_ggplot_2_vars(ems_sensitivity = ems_sens,
                                  ems_specificity = ems_spec,
                                  var_signal = var_signal,
                                  var_group = var_group,
                                  response = response,
                                  palette = palette,
                                  ttl = ttl)

  } else if (!is.null(var_facet)) {

    facet_levels <- levels(model$model[var_facet][1,])
    f <- formula(paste0("~ cut | ", var_group, " + ", var_facet, " + ", var_signal))

    ems_sens_1 <- emmeans::emmeans(model, f, at = list(target = signal_levels[1], condition = facet_levels[1]), mode = "cum.prob")
    ems_spec_1 <- emmeans::emmeans(model, f, at = list(target = signal_levels[2], condition = facet_levels[1]), mode = "exc.prob")

    plot1 <- roc_ggplot_2_vars(ems_sensitivity = ems_sens_1,
                               ems_specificity = ems_spec_1,
                               var_signal = var_signal,
                               var_group = var_group,
                               response = response,
                               palette = palette,
                               ttl = stringr::str_to_title(facet_levels[1]))

    ems_sens_2 <- emmeans::emmeans(model, f, at = list(target = signal_levels[1], condition = facet_levels[2]), mode = "cum.prob")
    ems_spec_2 <- emmeans::emmeans(model, f, at = list(target = signal_levels[2], condition = facet_levels[2]), mode = "exc.prob")

    plot2 <- roc_ggplot_2_vars(ems_sensitivity = ems_sens_2,
                               ems_specificity = ems_spec_2,
                               var_signal = var_signal,
                               var_group = var_group,
                               response = response,
                               palette = palette,
                               ttl = stringr::str_to_title(facet_levels[2]))

    out_plot <- (plot1 + plot2) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(title = ttl, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = 0.5)))

  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

