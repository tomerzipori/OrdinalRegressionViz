#' ROC curve of a frequentist ordinal regression model
#'
#' Draws a receiver operating characteristic (ROC) curve — hit rate against
#' false-alarm rate across the response thresholds — from a cumulative link
#' model fitted with [ordinal::clm()] or [ordinal::clmm()]. One curve is
#' drawn per level of `var_group`, with confidence bands computed from
#' [emmeans::emmeans()] cumulative-probability estimates; when `var_facet`
#' is given, side-by-side ROC panels are drawn for its two levels.
#'
#' @param model A `clm` or `clmm` ordinal regression model.
#' @param var_signal Name (string) of the 2-level factor playing the role of
#'   the SDT signal, e.g. old/new or true/fake.
#' @param var_group Name (string) of a 2-level factor covariate; one ROC
#'   curve is drawn per level.
#' @param var_facet Optional name (string) of a second 2-level factor
#'   covariate; one panel is drawn per level.
#' @param CI Confidence level of the intervals (between 0 and 1).
#' @param palette_groups Integer index of a qualitative palette passed to
#'   [ggplot2::scale_fill_brewer()] for the confidence bands.
#' @param palette_thresholds Integer index of a diverging palette passed to
#'   [ggplot2::scale_fill_brewer()] for the threshold points.
#' @param group_labels Optional character vector of length 2 with legend
#'   labels for the levels of `var_group`. Defaults to its title-cased
#'   levels.
#' @param facet_labels Optional character vector of length 2 with panel
#'   titles for the levels of `var_facet`. Defaults to its title-cased
#'   levels.
#' @param show_empirical If `TRUE`, the observed (empirical) hit and
#'   false-alarm rates are overlaid as crosses — a quick visual check of how
#'   well the model-implied ROC matches the data.
#' @param ttl Plot title.
#' @return A [ggplot2::ggplot] object (or a patchwork of two panels when
#'   `var_facet` is given).
#' @seealso [bayesian_roc_plot()] for [brms::brm()] models, and
#'   [SDT_distributions_plot()] for the corresponding latent distributions.
#' @examplesIf requireNamespace("ordinal", quietly = TRUE)
#' fit <- ordinal::clm(value ~ target * time, data = sdt_ratings)
#' roc_plot(fit, var_signal = "target", var_group = "time")
#' @export
roc_plot <- function(model,
                     var_signal,
                     var_group,
                     var_facet = NULL,
                     CI = 0.95,
                     palette_groups = 2,
                     palette_thresholds = 2,
                     group_labels = NULL,
                     facet_labels = NULL,
                     show_empirical = FALSE,
                     ttl = "") {
  check_clm(model)
  check_prob(CI, "CI")

  model_data <- model$model
  sig <- check_factor_var(model_data, var_signal, "var_signal")
  grp <- check_factor_var(model_data, var_group, "var_group")
  if (!is.null(var_facet)) {
    check_factor_var(model_data, var_facet, "var_facet")
  }
  group_labels <- check_labels2(group_labels, levels(grp), "group_labels")
  signal_levels <- levels(sig)

  if (is.null(var_facet)) {
    f <- stats::formula(paste0("~ cut | ", var_group, " + ", var_signal))
    ems_sens <- emmeans::emmeans(
      model, f,
      at = stats::setNames(list(signal_levels[1]), var_signal),
      mode = "cum.prob", level = CI
    )
    ems_spec <- emmeans::emmeans(
      model, f,
      at = stats::setNames(list(signal_levels[2]), var_signal),
      mode = "exc.prob", level = CI
    )
    out_plot <- roc_ggplot_2_vars(
      ems_sensitivity = ems_sens,
      ems_specificity = ems_spec,
      var_group = var_group,
      CI = CI,
      palette_groups = palette_groups,
      palette_thresholds = palette_thresholds,
      group_labels = group_labels,
      empirical = if (show_empirical) {
        empirical_roc_points(model_data, response_name(model_data), var_signal, var_group)
      },
      ttl = ttl
    )
    return(out_plot)
  }

  facet_levels <- levels(model_data[[var_facet]])
  facet_labels <- check_labels2(facet_labels, facet_levels, "facet_labels")
  f <- stats::formula(
    paste0("~ cut | ", var_group, " + ", var_facet, " + ", var_signal)
  )

  panels <- lapply(1:2, function(i) {
    at_sens <- stats::setNames(
      list(signal_levels[1], facet_levels[i]), c(var_signal, var_facet)
    )
    at_spec <- stats::setNames(
      list(signal_levels[2], facet_levels[i]), c(var_signal, var_facet)
    )
    roc_ggplot_2_vars(
      ems_sensitivity = emmeans::emmeans(model, f, at = at_sens, mode = "cum.prob", level = CI),
      ems_specificity = emmeans::emmeans(model, f, at = at_spec, mode = "exc.prob", level = CI),
      var_group = var_group,
      CI = CI,
      palette_groups = palette_groups,
      palette_thresholds = palette_thresholds,
      group_labels = group_labels,
      empirical = if (show_empirical) {
        facet_data <- model_data[model_data[[var_facet]] == facet_levels[i], , drop = FALSE]
        empirical_roc_points(facet_data, response_name(model_data), var_signal, var_group)
      },
      ttl = facet_labels[i]
    )
  })

  (panels[[1]] + panels[[2]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = ttl,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = 0.5)
      )
    )
}


#' ROC curve panel for one pair of 2-level factors (frequentist)
#'
#' Internal workhorse of [roc_plot()]: turns emmeans cumulative- and
#' exceedance-probability grids into a single ROC panel.
#'
#' @param ems_sensitivity `emmGrid` of cumulative probabilities
#'   (`mode = "cum.prob"`) at the first (signal) level of the SDT variable.
#' @param ems_specificity `emmGrid` of exceedance probabilities
#'   (`mode = "exc.prob"`) at the second (noise) level of the SDT variable.
#' @param empirical Optional data frame of observed hit/false-alarm rates
#'   (from `empirical_roc_points()`) to overlay.
#' @inheritParams roc_plot
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
roc_ggplot_2_vars <- function(ems_sensitivity,
                              ems_specificity,
                              var_group,
                              CI = 0.95,
                              palette_groups = 2,
                              palette_thresholds = 2,
                              group_labels = NULL,
                              empirical = NULL,
                              ttl = "") {
  df_sens_raw <- as.data.frame(ems_sensitivity)
  group_levels <- levels(df_sens_raw[[var_group]])
  if (is.null(group_labels)) {
    group_labels <- stringr::str_to_title(group_levels)
  }

  # Column names differ across emmeans versions ("cumprob" vs "cum.prob",
  # asymptotic vs t-based intervals).
  df_Sensitivity <- df_sens_raw |>
    dplyr::rename(dplyr::any_of(c(
      Sensitivity = "cumprob", Sensitivity = "cum.prob",
      Sensitivity_low = "asymp.LCL", Sensitivity_low = "lower.CL",
      Sensitivity_high = "asymp.UCL", Sensitivity_high = "upper.CL"
    ))) |>
    dplyr::select(cut, !!dplyr::sym(var_group), dplyr::starts_with("Sensitivity"))

  df_Specificity <- as.data.frame(ems_specificity) |>
    dplyr::rename(dplyr::any_of(c(
      Specificity = "exc.prob",
      Specificity_low = "asymp.LCL", Specificity_low = "lower.CL",
      Specificity_high = "asymp.UCL", Specificity_high = "upper.CL"
    ))) |>
    dplyr::select(cut, !!dplyr::sym(var_group), dplyr::starts_with("Specificity"))

  edges <- data.frame(
    group = factor(rep(group_levels, each = 2), levels = group_levels),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )
  names(edges)[1] <- var_group

  roc_data <- df_Sensitivity |>
    dplyr::full_join(
      df_Specificity,
      by = dplyr::join_by(cut, !!dplyr::sym(var_group))
    ) |>
    dplyr::rows_append(edges) |>
    dplyr::mutate(
      FAR = 1 - Specificity,
      FAR_low = 1 - Specificity_low,
      FAR_high = 1 - Specificity_high
    ) |>
    dplyr::arrange(!!dplyr::sym(var_group), Sensitivity)

  bands <- roc_data |>
    bayesian_CI_bands(
      "FAR", "FAR_low", "FAR_high",
      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
      group = c(var_group)
    )

  out_plot <- ggplot2::ggplot(roc_data, ggplot2::aes(FAR, Sensitivity)) +
    ggplot2::geom_polygon(
      ggplot2::aes(x, y, fill = !!dplyr::sym(var_group)),
      data = bands, alpha = 0.4
    ) +
    ggplot2::scale_fill_brewer(
      type = "qual", palette = palette_groups,
      name = stringr::str_to_title(var_group),
      labels = group_labels
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_path(
      ggplot2::aes(linetype = !!dplyr::sym(var_group)),
      linewidth = 1, show.legend = FALSE
    ) +
    ggplot2::scale_linetype_manual(values = c("dashed", "solid")) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = FAR_low, xmax = FAR_high),
      color = "grey40",
      data = \(d) tidyr::drop_na(d, FAR_low, FAR_high)
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = Sensitivity_low, ymax = Sensitivity_high),
      color = "grey40",
      data = \(d) tidyr::drop_na(d, Sensitivity_low, Sensitivity_high)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = ordered(cut)),
      shape = 21, size = 3,
      data = \(d) tidyr::drop_na(d, cut)
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +
    ggplot2::scale_fill_brewer("Threshold",
      type = "div", palette = palette_thresholds,
      na.translate = FALSE
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      color = NULL, linetype = NULL,
      x = "False Alarm Rate", y = "Hit Rate",
      title = ttl,
      subtitle = sprintf("%g%% Confidence Interval", 100 * CI)
    ) +
    ggplot2::coord_fixed() +
    serif_titles(subtitle_size = 14)

  if (!is.null(empirical)) {
    out_plot <- out_plot +
      ggplot2::geom_point(data = empirical, shape = 4, size = 2.2, stroke = 0.9)
  }
  out_plot
}
