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

    ems_sens_1 <- emmeans(model, f, at = list(target = signal_levels[1], condition = facet_levels[1]), mode = "cum.prob")
    ems_spec_1 <- emmeans(model, f, at = list(target = signal_levels[2], condition = facet_levels[1]), mode = "exc.prob")

    plot1 <- roc_ggplot_2_vars(ems_sensitivity = ems_sens_1,
                               ems_specificity = ems_spec_1,
                               var_signal = var_signal,
                               var_group = var_group,
                               response = response,
                               palette = palette,
                               ttl = str_to_title(facet_levels[1]))

    ems_sens_2 <- emmeans(model, f, at = list(target = signal_levels[1], condition = facet_levels[2]), mode = "cum.prob")
    ems_spec_2 <- emmeans(model, f, at = list(target = signal_levels[2], condition = facet_levels[2]), mode = "exc.prob")

    plot2 <- roc_ggplot_2_vars(ems_sensitivity = ems_sens_2,
                               ems_specificity = ems_spec_2,
                               var_signal = var_signal,
                               var_group = var_group,
                               response = response,
                               palette = palette,
                               ttl = str_to_title(facet_levels[2]))

    out_plot <- (plot1 + plot2) +
      plot_layout(guides = "collect") +
      plot_annotation(title = ttl, theme = theme(plot.title = element_text(size = 20, family = "serif", hjust = 0.5)))

  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

