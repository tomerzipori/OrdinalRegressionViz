#' Latent SDT distributions of a frequentist ordinal regression model
#'
#' Plots the latent ("perceived signal") distributions implied by a
#' cumulative link model fitted with [ordinal::clm()] or [ordinal::clmm()],
#' together with the estimated response thresholds (criteria), in the style
#' of signal detection theory (SDT). One panel is drawn for every level of
#' `var_group`, optionally split by a third variable (`var_facet`).
#'
#' The latent distribution follows from the model's link function: normal
#' for `probit` (the classic SDT display), logistic for `logit`, Cauchy for
#' `cauchit`, and Gumbel for `cloglog`/`loglog`. All latent distributions
#' have unit scale (scale/disc effects are not supported). Models without
#' some interaction terms are supported: terms the model does not contain
#' are simply treated as zero.
#'
#' @inheritParams roc_plot
#' @param plot_limits Numeric vector of length 2: x-axis limits of the
#'   latent scale.
#' @param palette Integer index of a diverging palette passed to
#'   [ggplot2::scale_color_brewer()] for the threshold lines.
#' @param alpha Opacity of the threshold lines.
#' @param signal_labels Optional character vector of length 2 with legend
#'   labels for the two latent distributions (reference level first).
#'   Defaults to the title-cased levels of `var_signal`.
#' @param group_labels Optional character vector of length 2 with panel
#'   titles for the levels of `var_group`. Defaults to its title-cased
#'   levels.
#' @param x_label Label of the x (latent) axis.
#' @return A [patchwork][patchwork::patchwork-package] object combining the
#'   panels.
#' @seealso [bayesian_SDT_distribution_plot()] for [brms::brm()] models, and
#'   [roc_plot()] for the corresponding ROC curves.
#' @examplesIf requireNamespace("ordinal", quietly = TRUE)
#' fit <- ordinal::clm(value ~ target * time, data = sdt_ratings, link = "probit")
#' SDT_distributions_plot(fit, var_signal = "target", var_group = "time")
#' @export
SDT_distributions_plot <- function(model,
                                   var_signal,
                                   var_group,
                                   var_facet = NULL,
                                   plot_limits = c(-3, 5),
                                   palette = 9,
                                   alpha = 0.8,
                                   signal_labels = NULL,
                                   group_labels = NULL,
                                   facet_labels = NULL,
                                   ttl = "",
                                   x_label = "Latent signal") {
  check_clm(model)

  model_data <- model$model
  sig <- check_factor_var(model_data, var_signal, "var_signal")
  grp <- check_factor_var(model_data, var_group, "var_group")
  if (!is.null(var_facet)) {
    check_factor_var(model_data, var_facet, "var_facet")
  }
  signal_labels <- check_labels2(signal_labels, levels(sig), "signal_labels")
  group_labels <- check_labels2(group_labels, levels(grp), "group_labels")

  dens <- link_density(model$link)

  thresholds <- model$alpha
  coefs <- stats::coef(model)
  beta_names <- setdiff(names(coefs), names(thresholds))

  term_sig <- paste0(var_signal, levels(sig)[2])
  term_grp <- paste0(var_group, levels(grp)[2])
  term_facet <- if (is.null(var_facet)) {
    NULL
  } else {
    paste0(var_facet, levels(model_data[[var_facet]])[2])
  }
  match_term(beta_names, term_sig, required = TRUE)

  build_panel <- function(group_on, facet_on, panel_ttl) {
    on_terms <- c(
      if (group_on) term_grp,
      if (facet_on) term_facet
    )
    specs <- cell_specs(term_sig, on_terms)
    shift <- sum_matched_terms(coefs, beta_names, specs$shifts)
    signal_mean <- sum_matched_terms(coefs, beta_names, specs$mean_signal)

    SDT_dist_ggplot(
      ref_mean = 0,
      group_mean = signal_mean,
      thresholds = thresholds - shift,
      dens = dens,
      signal_labels = signal_labels,
      palette = palette,
      alpha = alpha,
      plot_limits = plot_limits,
      ttl = panel_ttl,
      x_label = x_label
    )
  }

  if (is.null(var_facet)) {
    out_plot <- (
      build_panel(FALSE, FALSE, panel_ttl = group_labels[1]) /
        build_panel(TRUE, FALSE, panel_ttl = group_labels[2])
    ) +
      patchwork::plot_layout(guides = "collect")
  } else {
    facet_labels <- check_labels2(
      facet_labels, levels(model_data[[var_facet]]), "facet_labels"
    )
    row_label <- function(lbl) {
      patchwork::wrap_elements(
        panel = grid::textGrob(lbl, gp = grid::gpar(fontface = "bold", fontfamily = "serif"))
      )
    }
    out_plot <- (
      row_label(facet_labels[1]) /
        (build_panel(FALSE, FALSE, group_labels[1]) | build_panel(TRUE, FALSE, group_labels[2])) /
        row_label(facet_labels[2]) /
        (build_panel(FALSE, TRUE, group_labels[1]) | build_panel(TRUE, TRUE, group_labels[2]))
    ) +
      patchwork::plot_layout(guides = "collect", heights = c(.17, 1, .17, 1), nrow = 4)
  }

  out_plot +
    patchwork::plot_annotation(
      title = ttl,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = 0.5)
      )
    )
}


#' Single latent-distribution panel (frequentist)
#'
#' Internal workhorse of [SDT_distributions_plot()]: draws the reference and
#' signal latent distributions of one design cell together with its response
#' thresholds.
#'
#' @param ref_mean Location of the reference ("noise") distribution.
#' @param group_mean Location of the signal distribution.
#' @param thresholds Named numeric vector of response thresholds.
#' @param dens Density function of the latent distribution, see
#'   `link_density()`.
#' @param signal_labels Character vector of length 2 with the legend labels
#'   (reference first).
#' @inheritParams SDT_distributions_plot
#' @return A [ggplot2::ggplot] object.
#' @keywords internal
SDT_dist_ggplot <- function(ref_mean = 0,
                            group_mean = 1,
                            thresholds,
                            dens = link_density("probit"),
                            signal_labels = c("Noise", "Signal"),
                            palette = 9,
                            alpha = 0.7,
                            plot_limits = c(-3, 5),
                            ttl = "",
                            x_label = "") {
  threshold_names <- factor(names(thresholds), levels = names(thresholds))

  ggplot2::ggplot() +
    # Reference ("noise") distribution
    ggplot2::stat_function(
      ggplot2::aes(linetype = signal_labels[1]),
      fun = dens,
      args = list(location = ref_mean, scale = 1),
      linewidth = 1
    ) +
    # Signal distribution
    ggplot2::stat_function(
      ggplot2::aes(linetype = signal_labels[2]),
      fun = dens,
      args = list(location = group_mean, scale = 1),
      linewidth = 1
    ) +
    # Thresholds
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = thresholds, color = threshold_names),
      linewidth = 1.5, alpha = alpha
    ) +
    ggplot2::scale_linetype_manual(
      breaks = signal_labels,
      values = c("solid", "dashed")
    ) +
    ggplot2::scale_color_brewer("Threshold",
      type = "div", palette = palette,
      labels = stringr::str_replace(names(thresholds), pattern = "\\|", replacement = " | ")
    ) +
    ggplot2::labs(y = NULL, linetype = NULL, x = x_label, title = ttl) +
    ggplot2::expand_limits(x = plot_limits, y = 0.45) +
    ggplot2::scale_x_continuous(
      breaks = seq(plot_limits[1], plot_limits[2], 1),
      labels = seq(plot_limits[1], plot_limits[2], 1)
    ) +
    ggplot2::theme_classic()
}
