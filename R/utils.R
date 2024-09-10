####### Helper functions for Bayesian ROC plots #######

#' Cumulative sum for posterior::rvar objects
#'
#' This function is a simple loop to calculate a cumulative sum of rvar objects
#' @param rvars a posterior rvar object
#' @return numeric vector
#' @keywords internal
rvar_cumsum <- function(rvars) {
  out <- rvars
  for (i in c(2:length(rvars))) {
    out[i] <- out[i-1] + rvars[i]
  }
  return(out)
}


#' Calculating credible interval bands for ROC curves
#'
#' This is a helper function for the main Bayesian ROC curve plotting function, calculating the CI bands around the curve
#' @param xmin lower bound of CI on the x axis
#' @param xmax upper bound of CI on the x axis
#' @param ymin lower bound of CI on the y axis
#' @param ymax upper bound of CI on the y axis
#' @return data frame
#' @keywords internal
bayesian_CI_bands <- function(data,
                              x, xmin, xmax,
                              y, ymin, ymax,
                              group = NULL) {
  bands_low <- data |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(group))) |>
    dplyr::reframe(
      x = c(dplyr::pick(dplyr::all_of(xmin))[[1]], dplyr::pick(dplyr::all_of(x))[[1]]),
      y = c(dplyr::pick(dplyr::all_of(y))[[1]], dplyr::pick(dplyr::all_of(ymin))[[1]]),
    ) |>
    dplyr::mutate(
      x = dplyr::coalesce(x, y),
      y = dplyr::coalesce(y, x)
    ) |>
    dplyr::arrange(x) |>
    dplyr::distinct()


  bands_high <- data |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(group))) |>
    dplyr::reframe(
      x = c(dplyr::pick(dplyr::all_of(xmax))[[1]], dplyr::pick(dplyr::all_of(x))[[1]]),
      y = c(dplyr::pick(dplyr::all_of(y))[[1]], dplyr::pick(dplyr::all_of(ymax))[[1]]),
    ) |>
    dplyr::mutate(
      x = dplyr::coalesce(x, y),
      y = dplyr::coalesce(y, x)
    ) |>
    dplyr::arrange(desc(x)) |>
    dplyr::distinct()

  bands <- dplyr::bind_rows(bands_low, bands_high) |>
    dplyr::ungroup()

  return(bands)
}


#' ROC Curve for two 2-level categorical predictor Bayesian ordinal probit regression models
#'
#' This function makes a ROC curve for two categorical variables Bayesian ordinal probit regression models.
#' @param grid data frame of conditional probabilities posteriors
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param response variable name of the response
#' @param CI credible interval width between 0 and 1
#' @param centrality centrality measure for posterior, either "mean" or "median"
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param ttl plot's title
#' @return a ggplot2::ggplot plot object
#' @keywords internal
bayesian_roc_ggplot_2_vars <- function(grid,
                                       var_signal = "target",
                                       var_group = "time",
                                       response = "value",
                                       CI = 0.95,
                                       centrality = "mean",
                                       palette = 7,
                                       ttl = "") {

  grid_long <- grid |>
    tidyr::pivot_longer(cols = names(grid)[-c(1:2)],
                 names_to = "cut",
                 values_to = "prob") |>
    dplyr::mutate(cut = as.numeric(cut)) |>
    tidyr::pivot_wider(values_from = prob, names_from = !!dplyr::sym(var_signal)) |>
    dplyr::group_by(!!dplyr::sym(var_group)) |>
    dplyr::mutate(
      Sensitivity = dplyr::lag(rvar_cumsum(!!dplyr::sym(levels(grid[var_signal][1,])[1])), default = 0),
      Specificity = rev(rvar_cumsum(rev(!!dplyr::sym(levels(grid[var_signal][1,])[2])))),
      Threshold = paste0(dplyr::lag(cut), "|", cut)
    )|>
    dplyr::rows_append(data.frame(Sensitivity = 1, Specificity = 0)) |>
    dplyr::mutate(
      Threshold = ifelse(mean(Sensitivity) %in% c(0, 1), NA, Threshold)
    )

  grid_long_Sensitivity <- grid_long |>
    dplyr::group_by(Threshold, !!dplyr::sym(var_group)) |>
    tidyr::drop_na() |>
    dplyr::select(Threshold, !!dplyr::sym(var_group), Sensitivity) |>
    dplyr::mutate(Sensitivity_low = quantile(Sensitivity, probs = (1-CI)/2),
           Sensitivity_high = quantile(Sensitivity, probs = (1-CI)/2 + CI),
           Sensitivity = dplyr::case_when(centrality == "mean" ~ mean(Sensitivity),
                                          centrality == "median" ~ median(Sensitivity)))

  grid_long_Specificity <- grid_long |>
    dplyr::group_by(Threshold, !!dplyr::sym(var_group)) |>
    tidyr::drop_na() |>
    dplyr::select(Threshold, !!dplyr::sym(var_group), Specificity) |>
    dplyr::mutate(Specificity_low = quantile(Specificity, probs = (1-CI)/2),
           Specificity_high = quantile(Specificity, probs = (1-CI)/2 + CI),
           Specificity = dplyr::case_when(centrality == "mean" ~ mean(Specificity),
                                   centrality == "median" ~ median(Specificity)))

  edges <- data.frame(
    time = rep(levels(unique(grid[var_group])[1,]), each = length(levels(unique(grid[var_group])[1,]))),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )

  names(edges)[1] <- var_group

  roc_data_grid <- grid_long_Sensitivity |>
    dplyr::full_join(grid_long_Specificity, by = dplyr::join_by(Threshold, !!dplyr::sym(var_group))) |>
    dplyr::rows_append(edges) |>
    dplyr::mutate(FAR = 1 - Specificity,
           FAR_low = 1 - Specificity_low,
           FAR_high = 1 - Specificity_high) |>
    dplyr::arrange(!!dplyr::sym(var_group), Sensitivity)

  bands <- roc_data_grid |>
    bayesian_CI_bands("FAR", "FAR_low", "FAR_high",
                      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
                      group = c(var_group))

  ggplot2::ggplot(roc_data_grid, ggplot2::aes(FAR, Sensitivity)) +
    ggplot2::geom_polygon(ggplot2::aes(x, y, fill = !!dplyr::sym(var_group)), data = bands, alpha = 0.4) +
    ggplot2::scale_fill_manual(values = c("grey80", "grey44"), name = stringr::str_to_title(var_group), labels = stringr::str_to_title(levels(grid[var_group][1,]))) +
    ggnewscale::new_scale_fill() +

    ggplot2::geom_path(ggplot2::aes(linetype = !!dplyr::sym(var_group)), linewidth = 0.8, show.legend = F) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid")) +

    ggplot2::geom_linerange(ggplot2::aes(xmin = FAR_low, xmax = FAR_high), color = "grey40") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = Sensitivity_low, ymax = Sensitivity_high), color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(fill = ordered(Threshold)), shape = 21, size = 3,
               data = \(d) tidyr::drop_na(d, Threshold)) +

    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::expand_limits(x = c(0,1), y = c(0,1)) +
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, fill = var_group, x = "False Alarm Rate", y = "Hit Rate", title = ttl) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 19, family = "serif", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(size = 15, family = "serif", hjust = 0.5))

}


####### Helper functions for Bayesian latent distributions plots #######

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
#' @return a ggplot2::ggplot plot object
#' @keywords internal
bayesian_distributions_plot_2_vars <- function(b_model,
                                               var_signal = "target",
                                               var_group = "time",
                                               response = "value",
                                               group1_ttl = "",
                                               group2_ttl = "",
                                               response_scale = c(1:7),
                                               palette = 9,
                                               plot_range = c(-11, 11),
                                               ttl = "") {

  model_data <- b_model$data

  coefs <- tidybayes::get_variables(b_model)
  coefs <- coefs[stringr::str_detect(coefs, "b_")]
  coefs <- coefs[stringr::str_detect(coefs, paste0(var_signal, "|", var_group))]
  coef_signal <- coefs[stringr::str_detect(coefs, var_group, negate = T)]
  coef_group <- coefs[stringr::str_detect(coefs, var_signal, negate = T)]
  coef_interaction <- coefs[stringr::str_detect(coefs, ":")]


  criteria <- tidybayes::gather_rvars(b_model, b_Intercept[Response], !!dplyr::sym(coef_group[1]))
  criteria_control <- criteria[response_scale[-max(response_scale)],]
  criteria_feedback <- criteria_control |>
    dplyr::mutate(.value = .value - criteria$.value[criteria$.variable == coef_group[1]])

  plot_breaks <- ceiling((max(plot_range) - min(plot_range))/0.05) # 0.05 is completely arbitrary

  signal_dist_pre_true <- tidybayes::spread_draws(b_model, b_disc_Intercept) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  signal_dist_pre_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), b_disc_Intercept, !!dplyr::sym(coef_signal[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2]))) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!dplyr::sym(coef_signal[1]), sd = b_disc_Intercept * b_disc_targetfake)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  signal_dist_post_true <- tidybayes::spread_draws(b_model, b_disc_Intercept, !!dplyr::sym(coef_group[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2]))) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_timepost)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  signal_dist_post_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), !!dplyr::sym(coef_interaction[1]),
                                        b_disc_Intercept, !!dplyr::sym(coef_group[2]), !!dplyr::sym(coef_signal[2]), !!dplyr::sym(coef_interaction[2])) |>
    dplyr::rename(b_timepostXtargetfake = !!dplyr::sym(coef_interaction[1]),
           b_disc_timepostXtargetfake = !!dplyr::sym(coef_interaction[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2])),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2])),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = b_targetfake + b_timepostXtargetfake,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_timepostXtargetfake)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  plot1 <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_pre_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_pre_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_pre_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_pre_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = ttl, subtitle = ifelse(group1_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[1]), group1_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = .5),
                   plot.subtitle = ggplot2::element_text(size = 15, family = "serif"))

  plot2 <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_post_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_post_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_post_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_post_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = ifelse(group2_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[2]), group2_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 15, family = "serif"))

  return((plot1 / plot2) + patchwork::plot_layout(guides = "collect"))
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
#' @return a ggplot2::ggplot plot object
#' @keywords internal
bayesian_distributions_plot_3_vars <- function(b_model,
                                               var_signal = "target",
                                               var_group = "time",
                                               var_facet = "condition",
                                               response = "value",
                                               group1_ttl = "",
                                               group2_ttl = "",
                                               response_scale = c(1:7),
                                               palette = 9,
                                               plot_range = c(-11, 11),
                                               ttl = "") {


  model_data <- b_model$data
  coefs <- tidybayes::get_variables(b_model)
  coefs <- coefs[stringr::str_detect(coefs, "b_")]
  coefs <- coefs[stringr::str_detect(coefs, paste0(var_signal, "|", var_group, "|", var_facet))]
  coef_signal <- coefs[stringr::str_detect(coefs, paste0(var_group, "|", var_facet), negate = T)]
  coef_group <- coefs[stringr::str_detect(coefs, paste0(var_signal, "|", var_facet), negate = T)]
  coef_facet <- coefs[stringr::str_detect(coefs, paste0(var_signal, "|", var_group), negate = T)]
  coef_interactions <- coefs[stringr::str_detect(coefs, ":")]
  coef_signalXgroup <- coef_interactions[stringr::str_detect(coef_interactions, var_facet, negate = T)]
  coef_signalXfacet <- coef_interactions[stringr::str_detect(coef_interactions, var_group, negate = T)]
  coef_groupXfacet <- coef_interactions[stringr::str_detect(coef_interactions, var_signal, negate = T)]
  coef_signalXgroupXfacet <- coefs[stringr::str_detect(coefs, var_signal) & stringr::str_detect(coefs, var_group) & stringr::str_detect(coefs, var_facet)]

  library(tidyverse)
  library(tidybayes)
  library(patchwork)

  # calculating criteria
  criteria <- tidybayes::gather_rvars(b_model, b_Intercept[Response], !!dplyr::sym(coef_group[1]), !!dplyr::sym(coef_facet[1]), !!dplyr::sym(coef_groupXfacet[1]))
  criteria_control_pre <- criteria[response_scale[-max(response_scale)],]
  criteria_control_post <- criteria_control_pre |>
    dplyr::mutate(.value = .value - criteria$.value[criteria$.variable == coef_group[1]])
  criteria_feedback_pre <- criteria_control_pre |>
    dplyr::mutate(.value = .value - criteria$.value[criteria$.variable == coef_facet[1]])
  criteria_feedback_post <- criteria_control_pre |>
    dplyr::mutate(.value = .value - criteria$.value[criteria$.variable == coef_facet[1]] - criteria$.value[criteria$.variable == coef_group[1]] - criteria$.value[criteria$.variable == coef_groupXfacet[1]])

  plot_breaks <- ceiling((max(plot_range) - min(plot_range))/0.05) # 0.35 is completely arbitrary

  # Control
  ## Pre
  ### true
  signal_dist_control_pre_true <- tidybayes::spread_draws(b_model, b_disc_Intercept) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_control_pre_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), b_disc_Intercept, !!dplyr::sym(coef_signal[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2]))) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!dplyr::sym(coef_signal[1]), sd = b_disc_Intercept * b_disc_targetfake)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ## Post
  ### true
  signal_dist_control_post_true <- tidybayes::spread_draws(b_model, b_disc_Intercept, !!dplyr::sym(coef_group[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2]))) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_timepost)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_control_post_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), !!dplyr::sym(coef_signalXgroup[1]), b_disc_Intercept, !!dplyr::sym(coef_group[2]), !!dplyr::sym(coef_signal[2]), !!dplyr::sym(coef_signalXgroup[2])) |>
    dplyr::rename(b_timepostXtargetfake = !!dplyr::sym(coef_signalXgroup[1]),
           b_disc_timepostXtargetfake = !!dplyr::sym(coef_signalXgroup[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2])),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2])),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!dplyr::sym(coef_signal[1]) + b_timepostXtargetfake,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_timepostXtargetfake)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  # Feedback
  ## Pre
  ### true
  signal_dist_feedback_pre_true <- tidybayes::spread_draws(b_model, b_disc_Intercept, !!dplyr::sym(coef_facet[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_conditionfeedback = 1/exp(!!dplyr::sym(coef_facet[2]))) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0, sd = b_disc_Intercept * b_disc_conditionfeedback)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_feedback_pre_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), !!dplyr::sym(coef_signalXfacet[1]), b_disc_Intercept, !!dplyr::sym(coef_signal[2]), !!dplyr::sym(coef_facet[2]), !!dplyr::sym(coef_signalXfacet[2])) |>
    dplyr::rename(b_targetfakeXconditionfeedback = !!dplyr::sym(coef_signalXfacet[1]),
           b_disc_targetfakeXconditionfeedback = !!dplyr::sym(coef_signalXfacet[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2])),
           b_disc_conditionfeedback = 1/exp(!!dplyr::sym(coef_facet[2])),
           b_disc_targetfakeXconditionfeedback = 1/exp(b_disc_targetfakeXconditionfeedback)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!dplyr::sym(coef_signal[1]) + b_targetfakeXconditionfeedback,
                sd = b_disc_Intercept * b_disc_targetfake * b_disc_conditionfeedback * b_disc_targetfakeXconditionfeedback)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ## Post
  ### true
  signal_dist_feedback_post_true <- tidybayes::spread_draws(b_model, b_disc_Intercept, !!dplyr::sym(coef_group[2]), !!dplyr::sym(coef_facet[2]), !!dplyr::sym(coef_groupXfacet[2])) |>
    dplyr::rename(b_disc_timepostXconditionfeedback = !!dplyr::sym(coef_groupXfacet[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2])),
           b_disc_conditionfeedback = 1/exp(!!dplyr::sym(coef_facet[2])),
           b_disc_timepostXconditionfeedback = 1/exp(b_disc_timepostXconditionfeedback)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = 0,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_conditionfeedback * b_disc_timepostXconditionfeedback)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  ### fake
  signal_dist_feedback_post_fake <- tidybayes::spread_draws(b_model, !!dplyr::sym(coef_signal[1]), !!dplyr::sym(coef_signalXgroup[1]),
                                                 !!dplyr::sym(coef_signalXfacet[1]),
                                                 !!dplyr::sym(coef_signalXgroupXfacet[1]),
                                                 b_disc_Intercept, !!dplyr::sym(coef_group[2]), !!dplyr::sym(coef_signal[2]), !!dplyr::sym(coef_facet[2]),
                                                 !!dplyr::sym(coef_signalXgroup[2]), !!dplyr::sym(coef_groupXfacet[2]), !!dplyr::sym(coef_signalXfacet[2]),
                                                 !!dplyr::sym(coef_signalXgroupXfacet[2])) |>
    dplyr::rename(b_timepostXtargetfake = !!dplyr::sym(coef_signalXgroup[1]),
           b_targetfakeXconditionfeedback = !!dplyr::sym(coef_signalXfacet[1]),
           b_timepostXtargetfakeXconditionfeedback = !!dplyr::sym(coef_signalXgroupXfacet[1]),
           b_disc_timepostXtargetfake = !!dplyr::sym(coef_signalXgroup[2]),
           b_disc_timepostXconditionfeedback = !!dplyr::sym(coef_groupXfacet[2]),
           b_disc_targetfakeXconditionfeedback = !!dplyr::sym(coef_signalXfacet[2]),
           b_disc_timepostXtargetfakeXconditionfeedback = !!dplyr::sym(coef_signalXgroupXfacet[2])) |>
    dplyr::mutate(b_disc_Intercept = 1/exp(b_disc_Intercept),
           b_disc_timepost = 1/exp(!!dplyr::sym(coef_group[2])),
           b_disc_targetfake = 1/exp(!!dplyr::sym(coef_signal[2])),
           b_disc_conditionfeedback = 1/exp(!!dplyr::sym(coef_facet[2])),
           b_disc_timepostXtargetfake = 1/exp(b_disc_timepostXtargetfake),
           b_disc_timepostXconditionfeedback = 1/exp(b_disc_timepostXconditionfeedback),
           b_disc_targetfakeXconditionfeedback = 1/exp(b_disc_targetfakeXconditionfeedback),
           b_disc_timepostXtargetfakeXconditionfeedback = 1/exp(b_disc_timepostXtargetfakeXconditionfeedback)) |>
    dplyr::group_by(.draw) |>
    dplyr::reframe(
      x = seq(plot_range[1], plot_range[2], length = plot_breaks),
      d = dnorm(x, mean = !!dplyr::sym(coef_signal[1]) + b_timepostXtargetfake + b_targetfakeXconditionfeedback + b_timepostXtargetfakeXconditionfeedback,
                sd = b_disc_Intercept * b_disc_timepost * b_disc_targetfake * b_disc_conditionfeedback * b_disc_timepostXtargetfake * b_disc_timepostXconditionfeedback * b_disc_targetfakeXconditionfeedback * b_disc_timepostXtargetfakeXconditionfeedback)
    ) |>
    dplyr::ungroup() |>
    tidybayes::curve_interval(.along = x, .width = 0.9)

  # Control plot
  plot1_control <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_pre_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_control_pre_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_pre_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_control_pre_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control_pre) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = stringr::str_to_title(levels(model_data[,var_facet])[1]), subtitle = ifelse(group1_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[1]), group1_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 10, family = "serif"),
                   plot.title = ggplot2::element_text(size = 16, family = "serif", hjust = .5))

  plot2_control <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_post_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_control_post_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_control_post_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_control_post_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_control_post) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = ifelse(group2_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[2]), group2_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 10, family = "serif"))

  # Feedback plot
  plot1_feedback <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_pre_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_feedback_pre_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_pre_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_feedback_pre_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback_pre) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:12, " | ", 2:13),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = stringr::str_to_title(levels(model_data[,var_facet])[2]), subtitle = ifelse(group1_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[1]), group1_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 10, family = "serif"),
                   plot.title = ggplot2::element_text(size = 16, family = "serif", hjust = .5))

  plot2_feedback <- ggplot2::ggplot() +
    # Noise
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_post_true,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "True tweets"), data = signal_dist_feedback_post_true) +
    # Noise + Signal
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
                data = signal_dist_feedback_post_fake,
                fill = "grey", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(x, d, linetype = "Fake tweets"), data = signal_dist_feedback_post_fake) +
    # Thresholds
    tidybayes::stat_slab(ggplot2::aes(xdist = .value, fill = ordered(Response)),
              color = "gray", alpha = 0.6, key_glyph = "polygon",
              data = criteria_feedback_post) +
    # Theme and scales
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      labels = paste0(1:6, " | ", 2:7),
                      na.translate = FALSE) +
    ggplot2::labs(color = NULL, linetype = NULL, x = "", y = NULL, title = NULL, subtitle = ifelse(group2_ttl == "", stringr::str_to_title(levels(model_data[var_group][,1])[2]), group2_ttl)) +
    ggplot2::scale_x_continuous(limits = plot_range, breaks = seq(plot_range[1], plot_range[2]), labels = seq(plot_range[1], plot_range[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 10, family = "serif"))

  plot_out <- ((plot1_control / plot2_control) | (plot1_feedback / plot2_feedback)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = ttl, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = .5)))

  return(plot_out)
}


####### Helper function for Frequentist ROC plots #######

#' ROC Curve for 2-categorical variable models
#'
#' This function makes a ROC curve for two categorical variables ordinal regression models.
#' @param ems_sensitivity emmeans object containing conditional probabilities for the ordinal regression model + sensitivity calculation for each level
#' @param ems_sensitivity emmeans object containing conditional probabilities for the ordinal regression model + specificity calculation for each level
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param response variable name of the response
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param ttl plot's title
#' @return a ggplot2::ggplot plot object
#' @keywords internal
roc_ggplot_2_vars <- function(ems_sensitivity,
                              ems_specificity,
                              var_signal = "target",
                              var_group = "time",
                              response = "value",
                              palette = 2,
                              ttl = "") {

  levels <- unique(data.frame(ems_sensitivity)[var_group])

  df_Sensitivity <- ems_sensitivity |>
    data.frame() |>
    dplyr::rename(Sensitivity = cumprob,
           Sensitivity_low = asymp.LCL,
           Sensitivity_high = asymp.UCL) |>
    dplyr::select(cut, !!dplyr::sym(var_group), dplyr::starts_with("Sensitivity"))

  df_Specificity <- ems_specificity |>
    data.frame() |>
    dplyr::rename(Specificity = exc.prob,
           Specificity_low = asymp.LCL,
           Specificity_high = asymp.UCL) |>
    dplyr::select(cut, !!dplyr::sym(var_group), dplyr::starts_with("Specificity"))

  edges <- data.frame(
    group = rep(c(levels[1,], levels[2,]), each = 2),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )

  names(edges)[1] <- var_group

  roc_data <- df_Sensitivity |>
    dplyr::full_join(df_Specificity, by = dplyr::join_by(cut, !!dplyr::sym(var_group))) |>
    dplyr::rows_append(edges) |>
    dplyr::mutate(FAR = 1 - Specificity,
           FAR_low = 1 - Specificity_low,
           FAR_high = 1 - Specificity_high) |>
    dplyr::arrange(!!dplyr::sym(var_group), Sensitivity)

  bands <- roc_data |>
    bayesian_CI_bands("FAR", "FAR_low", "FAR_high",
                      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
                      group = c(var_group))


  out_plot <- ggplot2::ggplot(roc_data, ggplot2::aes(FAR, Sensitivity)) +
    ggplot2::geom_polygon(ggplot2::aes(x, y, fill = !!dplyr::sym(var_group)), data = bands,
                 alpha = 0.4) +
    ggplot2::scale_fill_brewer(type = "qual", palette = palette, name = stringr::str_to_title(var_group), labels = stringr::str_to_title(levels[,1])) +
    ggnewscale::new_scale_fill() +

    ggplot2::geom_path(ggplot2::aes(linetype = var_group), linewidth = 1, show.legend = F) +

    ggplot2::geom_linerange(ggplot2::aes(xmin = FAR_low, xmax = FAR_high), color = "grey40") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = Sensitivity_low, ymax = Sensitivity_high), color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(fill = ordered(cut)), shape = 21, size = 3,
               data = \(d) tidyr::drop_na(d, cut)) +

    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::expand_limits(x = c(0,1), y = c(0,1)) +
    ggplot2::scale_fill_brewer("Threshold", type = "div", palette = palette,
                      na.translate = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(color = NULL, linetype = NULL, fill = NULL, y = "HR", title = ttl, subtitle = "95% Confidence Interval") +
    ggplot2::coord_fixed() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(size = 14, family = "serif", hjust = 0.5))

  out_plot

}


####### Helper function for Frequentist latent distributions plots #######

#' Plot of latent normal distributions for two 2-level categorical variables ordinal probit regression models
#'
#' This function makes a latent distributions plot for two 2-level categorical variables ordinal regression models.
#' @param ref_mean mean of the normal distribution representing "noise"
#' @param group_mean mean of the normal distribution representing "noise + signal"
#' @param thresholds numeric vector representing thresholds of the ordinal model
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param alpha opacity of lines representing thresholds
#' @param plot_limits 2-integer vector representing x-axis limits of the plot
#' @param ttl plot's title
#' @param x_label label of the x axis
#' @return a ggplot2::ggplot plot object
#' @keywords internal
SDT_dist_ggplot <- function(ref_mean = 0,
                            group_mean = 1,
                            thresholds,
                            palette = 9,
                            alpha = 0.7,
                            plot_limits = c(-3, 5),
                            ttl = "",
                            x_label = "") {

  ggplot2::ggplot() +

    # Group1 + Noise
    ggplot2::stat_function(ggplot2::aes(linetype = "True"), fun = dnorm,
                  args = list(mean = ref_mean, sd = 1),
                  linewidth = 1) +
    # Group1 + Noise + Signal
    ggplot2::stat_function(ggplot2::aes(linetype = "Fake"), fun = dnorm,
                  args = list(mean = group_mean, sd = 1),
                  linewidth = 1) +
    # Thresholds
    ggplot2::geom_vline(ggplot2::aes(xintercept = thresholds, color = names(thresholds)),
               linewidth = 1.5, alpha = alpha) +
    ggplot2::scale_color_brewer("Threshold", type = "div", palette = palette,
                       labels = stringr::str_replace(names(thresholds), pattern = "\\|", replacement = " | ")) +
    ggplot2::labs(y = NULL, linetype = NULL, x = x_label, title = ttl) +
    ggplot2::expand_limits(x = plot_limits, y = 0.45) +
    ggplot2::scale_x_continuous(breaks = seq(plot_limits[1], plot_limits[2], 1), labels = seq(plot_limits[1], plot_limits[2], 1)) +
    ggplot2::theme_classic()
}

