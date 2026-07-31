#' ROC curve of a Bayesian ordinal regression model
#'
#' Draws a receiver operating characteristic (ROC) curve — hit rate against
#' false-alarm rate across the response thresholds — from a Bayesian
#' cumulative ("ordinal") regression model fitted with [brms::brm()]. One
#' curve is drawn per level of `var_group`, with curvewise credible bands;
#' when `var_facet` is given, side-by-side ROC panels are drawn for its two
#' levels.
#'
#' Population-level ("fixed effects only") predicted probabilities are used
#' (`re_formula = NA` in [brms::posterior_epred()]).
#'
#' @inheritParams bayesian_SDT_distribution_plot
#' @param CI Width of the credible intervals (between 0 and 1).
#' @param centrality Posterior point summary drawn as the curve: `"mean"`
#'   (default) or `"median"`.
#' @param palette_thresholds Integer index of a diverging palette passed to
#'   [ggplot2::scale_fill_brewer()] for the threshold points.
#' @param palette_curves Name of a viridis option (passed to
#'   [ggplot2::scale_fill_viridis_d()]) used for the credible bands.
#' @return A [ggplot2::ggplot] object (or a patchwork of two panels when
#'   `var_facet` is given).
#' @seealso [roc_plot()] for `ordinal::clm()` models.
#' @examples
#' \dontrun{
#' bayesian_roc_plot(
#'   b_fit,
#'   var_signal = "target", var_group = "time"
#' )
#' }
#' @export
bayesian_roc_plot <- function(b_model,
                              var_signal,
                              var_group,
                              var_facet = NULL,
                              CI = 0.95,
                              centrality = c("mean", "median"),
                              palette_thresholds = 7,
                              palette_curves = "viridis",
                              group_labels = NULL,
                              facet_labels = NULL,
                              ttl = "") {
  check_brmsfit(b_model)
  check_cumulative(b_model)
  check_prob(CI, "CI")
  centrality <- match.arg(centrality)

  model_data <- b_model$data
  check_factor_var(model_data, var_signal, "var_signal")
  grp <- check_factor_var(model_data, var_group, "var_group")
  if (!is.null(var_facet)) {
    check_factor_var(model_data, var_facet, "var_facet")
  }
  group_labels <- check_labels2(group_labels, levels(grp), "group_labels")

  grid <- unique(model_data[, c(var_signal, var_group, var_facet), drop = FALSE])
  # Population-level predictions
  probs_rvar <- posterior::rvar(
    brms::posterior_epred(b_model, newdata = grid, re_formula = NA)
  )
  grid <- cbind(grid, as.data.frame(probs_rvar))

  if (is.null(var_facet)) {
    out_plot <- bayesian_roc_ggplot_2_vars(
      grid,
      var_signal = var_signal, var_group = var_group,
      CI = CI, centrality = centrality,
      palette_thresholds = palette_thresholds,
      palette_curves = palette_curves,
      group_labels = group_labels, ttl = ttl
    )
    return(out_plot)
  }

  facet_levels <- levels(model_data[[var_facet]])
  facet_labels <- check_labels2(facet_labels, facet_levels, "facet_labels")

  panels <- lapply(1:2, function(i) {
    grid_i <- grid[grid[[var_facet]] == facet_levels[i], , drop = FALSE]
    grid_i[[var_facet]] <- NULL
    bayesian_roc_ggplot_2_vars(
      grid_i,
      var_signal = var_signal, var_group = var_group,
      CI = CI, centrality = centrality,
      palette_thresholds = palette_thresholds,
      palette_curves = palette_curves,
      group_labels = group_labels, ttl = facet_labels[i]
    )
  })

  (panels[[1]] + panels[[2]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = ttl,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20, hjust = 0.5, family = "serif")
      )
    )
}


#' ROC curve panel for one pair of 2-level factors (Bayesian)
#'
#' Internal workhorse of [bayesian_roc_plot()]: turns a grid of posterior
#' category-probability rvars into a single ROC panel.
#'
#' @param grid Data frame with the columns `var_signal`, `var_group` and one
#'   [posterior::rvar] column per response category.
#' @inheritParams bayesian_roc_plot
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
bayesian_roc_ggplot_2_vars <- function(grid,
                                       var_signal,
                                       var_group,
                                       CI = 0.95,
                                       centrality = "mean",
                                       palette_thresholds = 7,
                                       palette_curves = "viridis",
                                       group_labels = NULL,
                                       ttl = "") {
  signal_levels <- levels(grid[[var_signal]])
  group_levels <- levels(grid[[var_group]])
  if (is.null(group_labels)) {
    group_labels <- stringr::str_to_title(group_levels)
  }
  prob_cols <- setdiff(names(grid), c(var_signal, var_group))

  grid_long <- grid |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(prob_cols),
      names_to = "cut",
      values_to = "prob"
    ) |>
    dplyr::mutate(cut = as.numeric(cut)) |>
    tidyr::pivot_wider(values_from = prob, names_from = !!dplyr::sym(var_signal)) |>
    dplyr::group_by(!!dplyr::sym(var_group)) |>
    dplyr::mutate(
      Sensitivity = dplyr::lag(rvar_cumsum(!!dplyr::sym(signal_levels[1])), default = 0),
      Specificity = rev(rvar_cumsum(rev(!!dplyr::sym(signal_levels[2])))),
      Threshold = ifelse(
        is.na(dplyr::lag(cut)), NA_character_,
        paste0(dplyr::lag(cut), "|", cut)
      )
    )

  cuts <- sort(unique(grid_long$cut))
  threshold_levels <- paste0(utils::head(cuts, -1), "|", utils::tail(cuts, -1))

  grid_long_Sensitivity <- grid_long |>
    dplyr::group_by(Threshold, !!dplyr::sym(var_group)) |>
    tidyr::drop_na() |>
    dplyr::select(Threshold, !!dplyr::sym(var_group), Sensitivity) |>
    dplyr::mutate(
      Sensitivity_low = stats::quantile(Sensitivity, probs = (1 - CI) / 2),
      Sensitivity_high = stats::quantile(Sensitivity, probs = (1 - CI) / 2 + CI),
      Sensitivity = dplyr::case_when(
        centrality == "mean" ~ mean(Sensitivity),
        centrality == "median" ~ stats::median(Sensitivity)
      )
    )

  grid_long_Specificity <- grid_long |>
    dplyr::group_by(Threshold, !!dplyr::sym(var_group)) |>
    tidyr::drop_na() |>
    dplyr::select(Threshold, !!dplyr::sym(var_group), Specificity) |>
    dplyr::mutate(
      Specificity_low = stats::quantile(Specificity, probs = (1 - CI) / 2),
      Specificity_high = stats::quantile(Specificity, probs = (1 - CI) / 2 + CI),
      Specificity = dplyr::case_when(
        centrality == "mean" ~ mean(Specificity),
        centrality == "median" ~ stats::median(Specificity)
      )
    )

  edges <- data.frame(
    group = factor(rep(group_levels, each = 2), levels = group_levels),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )
  names(edges)[1] <- var_group

  roc_data_grid <- grid_long_Sensitivity |>
    dplyr::full_join(
      grid_long_Specificity,
      by = dplyr::join_by(Threshold, !!dplyr::sym(var_group))
    ) |>
    dplyr::rows_append(edges) |>
    dplyr::mutate(
      FAR = 1 - Specificity,
      FAR_low = 1 - Specificity_low,
      FAR_high = 1 - Specificity_high
    ) |>
    dplyr::arrange(!!dplyr::sym(var_group), Sensitivity)

  bands <- roc_data_grid |>
    bayesian_CI_bands(
      "FAR", "FAR_low", "FAR_high",
      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
      group = c(var_group)
    )

  ggplot2::ggplot(roc_data_grid, ggplot2::aes(FAR, Sensitivity)) +
    ggplot2::geom_polygon(
      ggplot2::aes(x, y, fill = !!dplyr::sym(var_group)),
      data = bands, alpha = 0.4
    ) +
    ggplot2::scale_fill_viridis_d(
      option = palette_curves,
      name = stringr::str_to_title(var_group),
      labels = group_labels
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_path(
      ggplot2::aes(linetype = !!dplyr::sym(var_group)),
      linewidth = 0.8, show.legend = FALSE
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
      ggplot2::aes(fill = ordered(Threshold, levels = threshold_levels)),
      shape = 21, size = 3,
      data = \(d) tidyr::drop_na(d, Threshold)
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +
    ggplot2::scale_fill_brewer("Threshold",
      type = "div", palette = palette_thresholds,
      na.translate = FALSE
    ) +
    ggplot2::labs(
      color = NULL, x = "False Alarm Rate", y = "Hit Rate", title = ttl
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    serif_titles(title_size = 19)
}
