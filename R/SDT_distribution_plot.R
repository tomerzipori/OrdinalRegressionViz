#' Latent normal distributions for two or three 2-level categorical predictor probit ordinal regression models
#'
#' This function makes a Latent normal distributions plot for two or three categorical variables ordinal regression models.
#' @param model ordinal regression model, as returned by clm or clmm functions
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param var_facet variable to facet by
#' @param response variable name of the response
#' @param plot_limits 2-integer vector representing x-axis limits of the plot
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param alpha opacity of lines representing thresholds
#' @param ttl plot's title
#' @param group1_ttl label for the first group subplot
#' @param group2_ttl label for the second group subplot
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
SDT_distributions_plot <- function(model,
                                   var_signal = "target",
                                   var_group = "time",
                                   var_facet = NULL,
                                   response = "value",
                                   plot_limits = c(-3, 5),
                                   palette = 9,
                                   alpha = 0.8,
                                   ttl = "",
                                   group1_ttl = "",
                                   group2_ttl = "",
                                   filename = NULL,
                                   path = getwd(),
                                   width = 2450,
                                   height = 1446) {

  model_data <- model$model
  model_coef <- names(coef(model))

  thresholds <- coef(model)[as.numeric(min(model_data[response][,1])):as.numeric(max(model_data[response][,1]))-1]
  b_group <- coef(model)[paste0(var_group, levels(model_data[var_group][,1])[2])]
  b_signal <- coef(model)[paste0(var_signal, levels(model_data[var_signal][,1])[2])]
  b_groupXsignal <- ifelse(is.null(var_facet),
                           coef(model)[model_coef[str_detect(model_coef, var_signal) & str_detect(model_coef, var_group)]],
                           coef(model)[model_coef[str_detect(model_coef, var_signal) & str_detect(model_coef, var_group) & str_detect(model_coef, var_facet, negate = T)]])

  if (!is.null(var_facet)) {

    b_facet <- coef(model)[paste0(var_facet, levels(model_data[var_facet][,1])[2])]
    b_groupXfacet <- coef(model)[model_coef[str_detect(model_coef, var_group) & str_detect(model_coef, var_facet) & str_detect(model_coef, var_signal, negate = T)]]
    b_signalXfacet <- coef(model)[model_coef[str_detect(model_coef, var_group, negate = T) & str_detect(model_coef, var_facet) & str_detect(model_coef, var_signal)]]
    b_groupXsignalXfacet <- coef(model)[model_coef[str_detect(model_coef, var_group) & str_detect(model_coef, var_facet) & str_detect(model_coef, var_signal)]]

  }

  # distribution means
  ## facet1
  ### group1
  thresholds_facet1_group1 <- thresholds
  mean_facet1_group1_noise <- 0
  mean_facet1_group1_signal <- b_signal

  facet1_group1_plot <- SDT_dist_ggplot(ref_mean = mean_facet1_group1_noise, group_mean = mean_facet1_group1_signal, thresholds = thresholds_facet1_group1,
                                        palette = palette, alpha = alpha, plot_limits = plot_limits, ttl = group1_ttl, x_label = "Obs. signal")

  ### group2
  thresholds_facet1_group2 <- thresholds - b_group
  mean_facet1_group2_noise <- 0
  mean_facet1_group2_signal <- b_signal + b_groupXsignal

  facet1_group2_plot <- SDT_dist_ggplot(ref_mean = mean_facet1_group2_noise, group_mean = mean_facet1_group2_signal, thresholds = thresholds_facet1_group2,
                                        palette = palette, alpha = alpha, plot_limits = plot_limits, ttl = group2_ttl, x_label = "Obs. signal")

  if (!is.null(var_facet)) {

    ## facet2
    ### group1
    thresholds_facet2_group1 <- thresholds - b_facet
    mean_facet2_group1_noise <- 0
    mean_facet2_group1_signal <- b_signal + b_signalXfacet

    facet2_group1_plot <- SDT_dist_ggplot(ref_mean = mean_facet2_group1_noise, group_mean = mean_facet2_group1_signal, thresholds = thresholds_facet2_group1,
                                          palette = palette, alpha = alpha, plot_limits = plot_limits, ttl = group1_ttl, x_label = "Obs. signal")

    ### group2
    thresholds_facet2_group2 <- thresholds - b_facet - b_group - b_groupXfacet
    mean_facet2_group2_noise <- 0
    mean_facet2_group2_signal <- b_signal + b_groupXsignal + b_signalXfacet + b_groupXsignalXfacet

    facet2_group2_plot <- SDT_dist_ggplot(ref_mean = mean_facet2_group2_noise, group_mean = mean_facet2_group2_signal, thresholds = thresholds_facet2_group2,
                                          palette = palette, alpha = alpha, plot_limits = plot_limits, x_label = "Obs. signal")

    row_label_1 <- wrap_elements(panel = ggpubr::text_grob(str_to_title(levels(model_data[var_facet][,1])[1]), face = "bold", family = "serif"))
    row_label_2 <- wrap_elements(panel = ggpubr::text_grob(str_to_title(levels(model_data[var_facet][,1])[2]), face = "bold", family = "serif"))

    out_plot <- (row_label_1 / (facet1_group1_plot | facet1_group2_plot) / row_label_2 / (facet2_group1_plot | facet2_group2_plot)) +
      plot_layout(guides = "collect", heights = c(.17,1,.17,1), nrow = 4) +
      plot_annotation(title = ttl, theme = theme(plot.title = element_text(size = 19, hjust = 0.5, family = "serif")))

  } else if (is.null(var_facet)) {

    out_plot <- (facet1_group1_plot / facet1_group2_plot) +
      plot_layout(guides = "collect") +
      plot_annotation(title = ttl, theme = theme(plot.title = element_text(family = "serif", size = 20, hjust = .5)))

  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

