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
    group_by(pick(any_of(group))) |>
    reframe(
      x = c(pick(all_of(xmin))[[1]], pick(all_of(x))[[1]]),
      y = c(pick(all_of(y))[[1]], pick(all_of(ymin))[[1]]),
    ) |>
    mutate(
      x = coalesce(x, y),
      y = coalesce(y, x)
    ) |>
    arrange(x) |>
    distinct()


  bands_high <- data |>
    group_by(pick(any_of(group))) |>
    reframe(
      x = c(pick(all_of(xmax))[[1]], pick(all_of(x))[[1]]),
      y = c(pick(all_of(y))[[1]], pick(all_of(ymax))[[1]]),
    ) |>
    mutate(
      x = coalesce(x, y),
      y = coalesce(y, x)
    ) |>
    arrange(desc(x)) |>
    distinct()

  bands <- bind_rows(bands_low, bands_high) |>
    ungroup()

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
#' @return a ggplot plot object
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
    pivot_longer(cols = names(grid)[-c(1:2)],
                 names_to = "cut",
                 values_to = "prob") |>
    mutate(cut = as.numeric(cut)) |>
    pivot_wider(values_from = prob, names_from = !!sym(var_signal)) |>
    group_by(!!sym(var_group)) |>
    mutate(
      Sensitivity = lag(rvar_cumsum(!!sym(levels(grid[var_signal][1,])[1])), default = 0),
      Specificity = rev(rvar_cumsum(rev(!!sym(levels(grid[var_signal][1,])[2])))),
      Threshold = paste0(lag(cut), "|", cut)
    )|>
    rows_append(data.frame(Sensitivity = 1, Specificity = 0)) |>
    mutate(
      Threshold = ifelse(mean(Sensitivity) %in% c(0, 1), NA, Threshold)
    )

  grid_long_Sensitivity <- grid_long |>
    group_by(Threshold, !!sym(var_group)) |>
    drop_na() |>
    select(Threshold, !!sym(var_group), Sensitivity) |>
    mutate(Sensitivity_low = quantile(Sensitivity, probs = (1-CI)/2),
           Sensitivity_high = quantile(Sensitivity, probs = (1-CI)/2 + CI),
           Sensitivity = case_when(centrality == "mean" ~ mean(Sensitivity),
                                   centrality == "median" ~ median(Sensitivity)))

  grid_long_Specificity <- grid_long |>
    group_by(Threshold, !!sym(var_group)) |>
    drop_na() |>
    select(Threshold, !!sym(var_group), Specificity) |>
    mutate(Specificity_low = quantile(Specificity, probs = (1-CI)/2),
           Specificity_high = quantile(Specificity, probs = (1-CI)/2 + CI),
           Specificity = case_when(centrality == "mean" ~ mean(Specificity),
                                   centrality == "median" ~ median(Specificity)))

  edges <- data.frame(
    time = rep(levels(unique(grid[var_group])[1,]), each = length(levels(unique(grid[var_group])[1,]))),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )

  roc_data_grid <- grid_long_Sensitivity |>
    full_join(grid_long_Specificity, by = join_by(Threshold, !!sym(var_group))) |>
    rows_append(edges) |>
    mutate(FAR = 1 - Specificity,
           FAR_low = 1 - Specificity_low,
           FAR_high = 1 - Specificity_high) |>
    arrange(!!sym(var_group), Sensitivity)

  bands <- roc_data_grid |>
    bayesian_CI_bands("FAR", "FAR_low", "FAR_high",
                      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
                      group = c(var_group))

  ggplot(roc_data_grid, aes(FAR, Sensitivity)) +
    geom_polygon(aes(x, y, fill = !!sym(var_group)), data = bands,
                 alpha = 0.4) +
    scale_fill_brewer(type = "qual", palette = 2) +
    ggnewscale::new_scale_fill() +

    geom_path(aes(linetype = !!sym(var_group)), linewidth = 0.8, show.legend = F) +
    scale_linetype_manual(values = c("dashed", "solid")) +

    geom_linerange(aes(xmin = FAR_low, xmax = FAR_high), color = "grey40") +
    geom_linerange(aes(ymin = Sensitivity_low, ymax = Sensitivity_high), color = "grey40") +
    geom_point(aes(fill = ordered(Threshold)), shape = 21, size = 3,
               data = \(d) drop_na(d, Threshold)) +

    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    expand_limits(x = c(0,1), y = c(0,1)) +
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      na.translate = FALSE) +
    labs(color = NULL, fill = var_group, y = "HR", title = ttl, subtitle = paste0(as.character(100*CI), "% Credible Interval")) +
    coord_fixed() +
    theme_classic() +
    theme(plot.title = element_text(size = 19, family = "serif", hjust = 0.5),
          plot.subtitle = element_text(size = 15, family = "serif", hjust = 0.5))

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
#' @return a ggplot plot object
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
    rename(Sensitivity = cumprob,
           Sensitivity_low = asymp.LCL,
           Sensitivity_high = asymp.UCL) |>
    select(cut, !!sym(var_group), starts_with("Sensitivity"))

  df_Specificity <- ems_specificity |>
    data.frame() |>
    rename(Specificity = exc.prob,
           Specificity_low = asymp.LCL,
           Specificity_high = asymp.UCL) |>
    select(cut, !!sym(var_group), starts_with("Specificity"))

  edges <- data.frame(
    group = rep(c(levels[1,], levels[2,]), each = 2),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )

  names(edges)[1] <- var_group

  roc_data <- df_Sensitivity |>
    full_join(df_Specificity, by = join_by(cut, !!sym(var_group))) |>
    rows_append(edges) |>
    mutate(FAR = 1 - Specificity,
           FAR_low = 1 - Specificity_low,
           FAR_high = 1 - Specificity_high) |>
    arrange(!!sym(var_group), Sensitivity)

  bands <- roc_data |>
    bayesian_CI_bands("FAR", "FAR_low", "FAR_high",
                      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
                      group = c(var_group))


  out_plot <- ggplot(roc_data, aes(FAR, Sensitivity)) +
    geom_polygon(aes(x, y, fill = !!sym(var_group)), data = bands,
                 alpha = 0.4) +
    scale_fill_brewer(type = "qual", palette = palette, name = str_to_title(var_group)) +
    ggnewscale::new_scale_fill() +

    geom_path(aes(linetype = var_group), linewidth = 1, show.legend = F) +

    geom_linerange(aes(xmin = FAR_low, xmax = FAR_high), color = "grey40") +
    geom_linerange(aes(ymin = Sensitivity_low, ymax = Sensitivity_high), color = "grey40") +
    geom_point(aes(fill = ordered(cut)), shape = 21, size = 3,
               data = \(d) drop_na(d, cut)) +

    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    expand_limits(x = c(0,1), y = c(0,1)) +
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      na.translate = FALSE) +
    theme_classic() +
    labs(color = NULL, linetype = NULL, fill = NULL, y = "HR", title = ttl, subtitle = "95% Confidence Interval") +
    coord_fixed() +
    theme(plot.title = element_text(size = 20, family = "serif", hjust = 0.5),
          plot.subtitle = element_text(size = 14, family = "serif", hjust = 0.5))

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
#' @return a ggplot plot object
#' @keywords internal
SDT_dist_ggplot <- function(ref_mean = 0,
                            group_mean = 1,
                            thresholds,
                            palette = 9,
                            alpha = 0.7,
                            plot_limits = c(-3, 5),
                            ttl = "",
                            x_label = "") {

  ggplot() +

    # Group1 + Noise
    stat_function(aes(linetype = "True"), fun = dnorm,
                  args = list(mean = ref_mean, sd = 1),
                  linewidth = 1) +
    # Group1 + Noise + Signal
    stat_function(aes(linetype = "Fake"), fun = dnorm,
                  args = list(mean = group_mean, sd = 1),
                  linewidth = 1) +
    # Thresholds
    geom_vline(aes(xintercept = thresholds, color = names(thresholds)),
               linewidth = 1.5, alpha = alpha) +
    scale_color_brewer("Threshold", type = "div", palette = palette,
                       labels = str_replace(names(thresholds), pattern = "\\|", replacement = " | ")) +
    labs(y = NULL, linetype = NULL, x = x_label, title = ttl) +
    expand_limits(x = plot_limits, y = 0.45) +
    scale_x_continuous(breaks = seq(plot_limits[1], plot_limits[2], 1), labels = seq(plot_limits[1], plot_limits[2], 1)) +
    theme_classic()
}

