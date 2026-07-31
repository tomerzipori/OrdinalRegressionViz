#' Latent SDT distributions of a Bayesian ordinal regression model
#'
#' Plots the latent ("perceived signal") distributions implied by a Bayesian
#' cumulative regression model, together with the posterior distributions of
#' the response thresholds (criteria), in the style of signal detection
#' theory (SDT). One panel is drawn for every level of `var_group`,
#' optionally split by a third variable (`var_facet`).
#'
#' The model is expected to be a [brms::brm()] fit with a cumulative family.
#' The latent distribution follows from the model's link function: normal
#' for `probit` (the classic SDT display), logistic for `logit`, Cauchy for
#' `cauchit`, and Gumbel for `cloglog`. When the discrimination parameter
#' `disc` is modeled (unequal-variance SDT), the posterior of each
#' distribution's scale is used; otherwise all latent distributions have
#' unit scale. Models without some interaction terms are supported: terms
#' the model does not contain are simply treated as zero.
#'
#' @param b_model A `brmsfit` cumulative regression model.
#' @param var_signal Name (string) of the 2-level factor playing the role of
#'   the SDT signal, e.g. old/new or true/fake.
#' @param var_group Name (string) of a 2-level factor covariate, e.g.
#'   pre/post. One panel is drawn per level.
#' @param var_facet Optional name (string) of a second 2-level factor
#'   covariate; panels are duplicated for each of its levels.
#' @param signal_labels Optional character vector of length 2 with legend
#'   labels for the two latent distributions (reference level first).
#'   Defaults to the title-cased levels of `var_signal`.
#' @param group_labels Optional character vector of length 2 with panel
#'   subtitles. Defaults to the title-cased levels of `var_group`.
#' @param facet_labels Optional character vector of length 2 with panel
#'   titles for the levels of `var_facet`. Defaults to its title-cased
#'   levels.
#' @param palette Integer index of a sequential palette passed to
#'   [ggplot2::scale_fill_brewer()] for the threshold slabs.
#' @param alpha Opacity of the threshold slabs.
#' @param plot_range Numeric vector of length 2: x-axis limits of the latent
#'   scale.
#' @param ci Width of the curvewise credible band around each density curve
#'   (between 0 and 1).
#' @param ndraws Optional number of posterior draws used to compute the
#'   density bands. The default (`NULL`) uses all draws; a value such as
#'   `500` makes the plot considerably faster with little visual change.
#' @param ttl Plot title.
#' @return A [patchwork][patchwork::patchwork-package] object combining the
#'   panels; it can be printed, modified with `&`, or saved with
#'   [ggplot2::ggsave()].
#' @seealso [SDT_distributions_plot()] for `ordinal::clm()` models, and
#'   [bayesian_roc_plot()] for the corresponding ROC curves.
#' @examples
#' \dontrun{
#' # b_fit <- brms::brm(
#' #   brms::bf(value ~ target * time, disc ~ target * time),
#' #   family = brms::cumulative("probit"),
#' #   data = sdt_ratings
#' # )
#' bayesian_SDT_distribution_plot(
#'   b_fit,
#'   var_signal = "target", var_group = "time",
#'   plot_range = c(-4, 6), ndraws = 500
#' )
#' }
#' @export
bayesian_SDT_distribution_plot <- function(b_model,
                                           var_signal,
                                           var_group,
                                           var_facet = NULL,
                                           signal_labels = NULL,
                                           group_labels = NULL,
                                           facet_labels = NULL,
                                           palette = 9,
                                           alpha = 0.6,
                                           plot_range = c(-11, 11),
                                           ci = 0.9,
                                           ndraws = NULL,
                                           ttl = "") {
  check_brmsfit(b_model)
  check_cumulative(b_model)
  check_prob(ci, "ci")
  dens <- link_density(b_model$family$link)

  model_data <- b_model$data
  sig <- check_factor_var(model_data, var_signal, "var_signal")
  grp <- check_factor_var(model_data, var_group, "var_group")
  if (!is.null(var_facet)) {
    check_factor_var(model_data, var_facet, "var_facet")
  }

  signal_labels <- check_labels2(signal_labels, levels(sig), "signal_labels")
  group_labels <- check_labels2(group_labels, levels(grp), "group_labels")

  pars <- tidybayes::get_variables(b_model)
  has_disc <- "b_disc_Intercept" %in% pars

  term_sig <- paste0(var_signal, levels(sig)[2])
  term_grp <- paste0(var_group, levels(grp)[2])
  term_facet <- if (is.null(var_facet)) {
    NULL
  } else {
    paste0(var_facet, levels(model_data[[var_facet]])[2])
  }

  # The signal main effect must exist, otherwise there is nothing to plot.
  match_coef(pars, term_sig, required = TRUE)

  x <- seq(plot_range[1], plot_range[2],
    length.out = ceiling((plot_range[2] - plot_range[1]) / 0.05)
  )
  thr_labels <- threshold_labels(model_data[[response_name(model_data)]])

  build_panel <- function(group_on, facet_on, subtitle, title = NULL) {
    on_terms <- c(
      if (group_on) term_grp,
      if (facet_on) term_facet
    )

    disc_noise <- if (has_disc) {
      c("b_disc_Intercept", match_coefs(pars, term_subsets(on_terms), disc = TRUE))
    } else {
      character(0)
    }
    disc_signal <- if (has_disc) {
      c("b_disc_Intercept", match_coefs(pars, term_subsets(c(term_sig, on_terms)), disc = TRUE))
    } else {
      character(0)
    }
    mean_signal_sets <- lapply(
      c(list(character(0)), term_subsets(on_terms)),
      function(s) c(term_sig, s)
    )
    mean_signal <- match_coefs(pars, mean_signal_sets)
    shift_coefs <- match_coefs(pars, term_subsets(on_terms))

    band_noise <- latent_density_band(
      b_model, character(0), disc_noise, x,
      width = ci, ndraws = ndraws, dens = dens
    )
    band_signal <- latent_density_band(
      b_model, mean_signal, disc_signal, x,
      width = ci, ndraws = ndraws, dens = dens
    )

    th <- threshold_rvars(b_model, shift_coefs)
    criteria <- dplyr::tibble(
      Response = seq_along(th),
      .value = th
    )

    ggplot2::ggplot() +
      # Reference ("noise") distribution
      ggplot2::geom_ribbon(
        ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
        data = band_noise, fill = "grey", alpha = 0.4
      ) +
      ggplot2::geom_line(
        ggplot2::aes(x, d, linetype = signal_labels[1]),
        data = band_noise
      ) +
      # Signal distribution
      ggplot2::geom_ribbon(
        ggplot2::aes(x = x, ymin = .lower, ymax = .upper),
        data = band_signal, fill = "grey", alpha = 0.4
      ) +
      ggplot2::geom_line(
        ggplot2::aes(x, d, linetype = signal_labels[2]),
        data = band_signal
      ) +
      # Thresholds
      tidybayes::stat_slab(
        ggplot2::aes(xdist = .value, fill = ordered(Response)),
        color = "gray", alpha = alpha, key_glyph = "polygon",
        data = criteria
      ) +
      ggplot2::scale_linetype_manual(
        breaks = signal_labels,
        values = c("solid", "dashed")
      ) +
      ggplot2::scale_fill_brewer("Threshold",
        type = "seq", palette = palette,
        labels = thr_labels, na.translate = FALSE
      ) +
      ggplot2::labs(
        color = NULL, linetype = NULL, x = NULL, y = NULL,
        title = title, subtitle = subtitle
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq(plot_range[1], plot_range[2])
      ) +
      # Crop instead of dropping data, so threshold slabs are not clipped
      # with warnings.
      ggplot2::coord_cartesian(xlim = plot_range) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 16, family = "serif", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 12, family = "serif")
      )
  }

  if (is.null(var_facet)) {
    out_plot <- (
      build_panel(FALSE, FALSE, subtitle = group_labels[1]) /
        build_panel(TRUE, FALSE, subtitle = group_labels[2])
    ) +
      patchwork::plot_layout(guides = "collect")
  } else {
    facet_labels <- check_labels2(
      facet_labels, levels(model_data[[var_facet]]), "facet_labels"
    )
    left <- build_panel(FALSE, FALSE, subtitle = group_labels[1], title = facet_labels[1]) /
      build_panel(TRUE, FALSE, subtitle = group_labels[2])
    right <- build_panel(FALSE, TRUE, subtitle = group_labels[1], title = facet_labels[2]) /
      build_panel(TRUE, TRUE, subtitle = group_labels[2])
    out_plot <- (left | right) +
      patchwork::plot_layout(guides = "collect")
  }

  out_plot +
    patchwork::plot_annotation(
      title = ttl,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20, family = "serif", hjust = 0.5)
      )
    )
}
