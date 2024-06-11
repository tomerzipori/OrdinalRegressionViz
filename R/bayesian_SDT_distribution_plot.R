#' Plot of latent normal distributions for two 2-level categorical variables Bayesian ordinal probit regression models
#'
#' This function makes a latent distributions plot for two 2-level categorical variables Bayesian ordinal regression models.
#' @param b_model Bayesian ordinal probit regression model - a brms object
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param response variable name of the response
#' @param response_scale integer vector of response scale levels, e.g 1 to 7
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param plot_range 2-integer vector representing x-axis limits of the plot
#' @param ttl plot's title
#' @return a ggplot plot object
#' @keywords internal
bayesian_distributions_plot_2_vars <- function(b_model,
                                               var_signal = "target",
                                               var_group = "time",
                                               response = "value",
                                               response_scale = c(1:7),
                                               palette = 9,
                                               plot_range = c(-11, 11),
                                               ttl = "") {

  coefs <- get_variables(b_model)
  coefs <- coefs[str_detect(coefs, "b_")]
  coefs <- coefs[str_detect(coefs, paste0(var_signal, "|", var_group))]
  coef_signal <- coefs[str_detect(coefs, var_group, negate = T)]
  coef_group <- coefs[str_detect(coefs, var_signal, negate = T)]
  coef_interaction <- coefs[str_detect(coefs, ":")]


  criteria <- gather_rvars(b_model, b_Intercept[Response], !!sym(coef_group[1]))
  criteria_control <- criteria[response_scale[-max(response_scale)],]
  criteria_feedback <- criteria_control |>
    mutate(.value = .value - criteria$.value[criteria$.variable == coef_group[1]])

  plot_breaks <- ceiling((max(plot_range) - min(plot_range))/0.05) # 0.05 is completely arbitrary

  signal_dist_pre_true <- spread_draws(b_model, b_disc_Intercept) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  signal_dist_pre_fake <- spread_draws(b_model, !!sym(coef_signal[1]), b_disc_Intercept, !!sym(coef_signal[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!sym(coef_signal[2]))) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!sym(coef_signal[1]), sd = b_disc_Intercept * b_disc_targetfake)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  signal_dist_post_true <- spread_draws(b_model, b_disc_Intercept, !!sym(coef_group[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!sym(coef_group[2]))) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_timepost)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  signal_dist_post_fake <- spread_draws(b_model, !!sym(coef_signal[1]), !!sym(coef_interaction[1]),
                                        b_disc_Intercept, !!sym(coef_group[2]), !!sym(coef_signal[2]), !!sym(coef_interaction[2])) |>
    rename(b_timepostXtargetfake = !!sym(coef_interaction[1]),
           b_disc_timepostXtargetfake = !!sym(coef_interaction[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(b_disc_timepost),
           b_disc_targetfake = 1/exp(b_disc_targetfake),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = b_targetfake + b_timepostXtargetfake,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_timepostXtargetfake)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  plot1 <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_pre_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_pre_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_pre_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_pre_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = 9,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = ttl, subtitle = "Pre-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, family = "serif", hjust = .5),
          plot.subtitle = element_text(size = 15, family = "serif"))

  plot2 <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_post_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_post_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_post_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_post_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = 9,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = "Post-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size = 15, family = "serif"))

  return((plot1 / plot2) + plot_layout(guides = "collect"))
}


#' Plot of latent normal distributions for three 2-level categorical variables Bayesian ordinal probit regression models
#'
#' This function makes a latent distributions plot for three 2-level categorical variables Bayesian ordinal regression models.
#' @param b_model Bayesian ordinal probit regression model - a brms object
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param var_facet variable to facet by
#' @param response variable name of the response
#' @param response_scale integer vector of response scale levels, e.g 1 to 7
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param plot_range 2-integer vector representing x-axis limits of the plot
#' @param ttl plot's title
#' @return a ggplot plot object
#' @keywords internal
bayesian_distributions_plot_3_vars <- function(b_model,
                                               var_signal = "target",
                                               var_group = "time",
                                               var_facet = "condition",
                                               response = "value",
                                               response_scale = c(1:7),
                                               palette = 9,
                                               plot_range = c(-11, 11),
                                               ttl = "") {


  model_data <- b_model$data
  coefs <- get_variables(b_model)
  coefs <- coefs[str_detect(coefs, "b_")]
  coefs <- coefs[str_detect(coefs, paste0(var_signal, "|", var_group, "|", var_facet))]
  coef_signal <- coefs[str_detect(coefs, paste0(var_group, "|", var_facet), negate = T)]
  coef_group <- coefs[str_detect(coefs, paste0(var_signal, "|", var_facet), negate = T)]
  coef_facet <- coefs[str_detect(coefs, paste0(var_signal, "|", var_group), negate = T)]
  coef_interactions <- coefs[str_detect(coefs, ":")]
  coef_signalXgroup <- coef_interactions[str_detect(coef_interactions, var_facet, negate = T)]
  coef_signalXfacet <- coef_interactions[str_detect(coef_interactions, var_group, negate = T)]
  coef_groupXfacet <- coef_interactions[str_detect(coef_interactions, var_signal, negate = T)]
  coef_signalXgroupXfacet <- coefs[str_detect(coefs, var_signal) & str_detect(coefs, var_group) & str_detect(coefs, var_facet)]

  library(tidyverse)
  library(tidybayes)
  library(patchwork)

  # calculating criteria
  criteria <- gather_rvars(b_model, b_Intercept[Response], !!sym(coef_group[1]), !!sym(coef_facet[1]), !!sym(coef_groupXfacet[1]))
  criteria_control_pre <- criteria[response_scale[-max(response_scale)],]
  criteria_control_post <- criteria_control_pre |>
    mutate(.value = .value - criteria$.value[criteria$.variable == coef_group[1]])
  criteria_feedback_pre <- criteria_control_pre |>
    mutate(.value = .value - criteria$.value[criteria$.variable == coef_facet[1]])
  criteria_feedback_post <- criteria_control_pre |>
    mutate(.value = .value - criteria$.value[criteria$.variable == coef_facet[1]] - criteria$.value[criteria$.variable == coef_group[1]] - criteria$.value[criteria$.variable == coef_groupXfacet[1]])

  plot_breaks <- ceiling((max(plot_range) - min(plot_range))/0.05) # 0.35 is completely arbitrary

  # Control
  ## Pre
  ### true
  signal_dist_control_pre_true <- spread_draws(b_model, b_disc_Intercept) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_control_pre_fake <- spread_draws(b_model, !!sym(coef_signal[1]), b_disc_Intercept, !!sym(coef_signal[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!sym(coef_signal[2]))) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!sym(coef_signal[1]), sd = b_disc_Intercept * b_disc_targetfake)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ## Post
  ### true
  signal_dist_control_post_true <- spread_draws(b_model, b_disc_Intercept, !!sym(coef_group[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!sym(coef_group[2]))) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_timepost)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_control_post_fake <- spread_draws(b_model, !!sym(coef_signal[1]), !!sym(coef_signalXgroup[1]), b_disc_Intercept, !!sym(coef_group[2]), !!sym(coef_signal[2]), !!sym(coef_signalXgroup[2])) |>
    rename(b_timepostXtargetfake = !!sym(coef_signalXgroup[1]),
           b_disc_timepostXtargetfake = !!sym(coef_signalXgroup[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!sym(coef_group[2])),
           b_disc_targetfake = 1/exp(!!sym(coef_signal[2])),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!sym(coef_signal[1]) + b_timepostXtargetfake,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_timepostXtargetfake)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  # Feedback
  ## Pre
  ### true
  signal_dist_feedback_pre_true <- spread_draws(b_model, b_disc_Intercept, !!sym(coef_facet[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_conditionfeedback = 1/exp(!!sym(coef_facet[2]))) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_conditionfeedback)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_feedback_pre_fake <- spread_draws(b_model, !!sym(coef_signal[1]), !!sym(coef_signalXfacet[1]), b_disc_Intercept, !!sym(coef_signal[2]), !!sym(coef_facet[2]), !!sym(coef_signalXfacet[2])) |>
    rename(b_targetfakeXconditionfeedback = !!sym(coef_signalXfacet[1]),
           b_disc_targetfakeXconditionfeedback = !!sym(coef_signalXfacet[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!sym(coef_signal[2])),
           b_disc_conditionfeedback = 1/exp(!!sym(coef_facet[2])),
           b_disc_targetfakeXconditionfeedback = 1/exp(b_disc_targetfakeXconditionfeedback)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!sym(coef_signal[1]) + b_targetfakeXconditionfeedback,
                sd = b_disc_Intercept * b_disc_targetfake * b_disc_conditionfeedback * b_disc_targetfakeXconditionfeedback)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ## Post
  ### true
  signal_dist_feedback_post_true <- spread_draws(b_model, b_disc_Intercept, !!sym(coef_group[2]), !!sym(coef_facet[2]), !!sym(coef_groupXfacet[2])) |>
    rename(b_disc_timepostXconditionfeedback = !!sym(coef_groupXfacet[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!sym(coef_group[2])),
           b_disc_conditionfeedback = 1/exp(!!sym(coef_facet[2])),
           b_disc_timepostXconditionfeedback = 1/exp(b_disc_timepostXconditionfeedback)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_conditionfeedback * b_disc_timepostXconditionfeedback)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_feedback_post_fake <- spread_draws(b_model, !!sym(coef_signal[1]), !!sym(coef_signalXgroup[1]),
                                                 !!sym(coef_signalXfacet[1]),
                                                 !!sym(coef_signalXgroupXfacet[1]),
                                                 b_disc_Intercept, !!sym(coef_group[2]), !!sym(coef_signal[2]), !!sym(coef_facet[2]),
                                                 !!sym(coef_signalXgroup[2]), !!sym(coef_groupXfacet[2]), !!sym(coef_signalXfacet[2]),
                                                 !!sym(coef_signalXgroupXfacet[2])) |>
    rename(b_timepostXtargetfake = !!sym(coef_signalXgroup[1]),
           b_targetfakeXconditionfeedback = !!sym(coef_signalXfacet[1]),
           b_timepostXtargetfakeXconditionfeedback = !!sym(coef_signalXgroupXfacet[1]),
           b_disc_timepostXtargetfake = !!sym(coef_signalXgroup[2]),
           b_disc_timepostXconditionfeedback = !!sym(coef_groupXfacet[2]),
           b_disc_targetfakeXconditionfeedback = !!sym(coef_signalXfacet[2]),
           b_disc_timepostXtargetfakeXconditionfeedback = !!sym(coef_signalXgroupXfacet[2])) |>
    mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!sym(coef_group[2])),
           b_disc_targetfake = 1/exp(!!sym(coef_signal[2])),
           b_disc_conditionfeedback = 1/exp(!!sym(coef_facet[2])),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake),
           b_disc_timepostXconditionfeedback = 1/exp(b_disc_timepostXconditionfeedback),
           b_disc_targetfakeXconditionfeedback = 1/exp(b_disc_targetfakeXconditionfeedback),
           b_disc_timepostXtargetfakeXconditionfeedback = 1/exp(b_disc_timepostXtargetfakeXconditionfeedback)) |>
    group_by(.draw) |>
    reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!sym(coef_signal[1]) + b_timepostXtargetfake + b_targetfakeXconditionfeedback + b_timepostXtargetfakeXconditionfeedback,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_conditionfeedback * b_disc_timepostXtargetfake * b_disc_timepostXconditionfeedback * b_disc_targetfakeXconditionfeedback * b_disc_timepostXtargetfakeXconditionfeedback)
    ) |>
    ungroup() |>
    curve_interval(.along = x, .width = 0.9)

  # Control plot
  plot1_control <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_pre_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_control_pre_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_pre_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_control_pre_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control_pre) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = str_to_title(levels(model_data[,var_facet])[1]), subtitle = "Pre-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size = 10, family = "serif"),
          plot.title = element_text(size = 16, family = "serif", hjust = .5))

  plot2_control <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_post_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_control_post_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_post_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_control_post_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control_post) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = "Post-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size = 10, family = "serif"))

  # Feedback plot
  plot1_feedback <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_pre_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_feedback_pre_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_pre_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_feedback_pre_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback_pre) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = str_to_title(levels(model_data[,var_facet])[2]), subtitle = "Pre-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size = 10, family = "serif"),
          plot.title = element_text(size = 16, family = "serif", hjust = .5))

  plot2_feedback <- ggplot() +
    # Noise
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_post_true,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "True tweets"), data = signal_dist_feedback_post_true) +
    # Noise + Signal
    geom_ribbon(aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_post_fake,
                fill = "grey", alpha = 0.4) +
    geom_line(aes(x, d, linetype = "Fake tweets"), data = signal_dist_feedback_post_fake) +
    # Thresholds
    stat_slab(aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback_post) +
    # Theme and scales
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = "Post-intervention") +
    scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(size = 10, family = "serif"))

  plot_out <- ((plot1_control / plot2_control) | (plot1_feedback / plot2_feedback)) +
    plot_layout(guides = "collect") +
    plot_annotation(title = ttl, theme = theme(plot.title = element_text(size = 20, family = "serif", hjust = .5)))

  return(plot_out)
}


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
                                                   plot_range = plot_range,
                                                   ttl = ttl)
  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

