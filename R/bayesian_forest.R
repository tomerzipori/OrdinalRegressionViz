#' Forest plot for a Bayesian random-effects meta-analysis
#'
#' Draws a forest plot — posterior densities and intervals per study, plus
#' the pooled effect — for a Bayesian random-effects meta-analysis fitted
#' with [brms::brm()]. The model is expected to contain nested grouping
#' terms of the form `(1 | author/study)` (defaults: `Author`, `Study`);
#' densities are colored by the higher-level (author) grouping.
#'
#' The layout is adapted from Harrer, Cuijpers, Furukawa, & Ebert (2021),
#' *Doing Meta-Analysis with R: A Hands-On Guide*
#' (<https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/>).
#'
#' @param b_meta_analysis A `brmsfit` random-effects meta-analysis with
#'   nested author/study grouping terms.
#' @param author_var Name (string) of the higher-level grouping variable.
#' @param study_var Name (string) of the study-level grouping variable
#'   nested in `author_var`.
#' @param palette Integer: number of colors drawn from the Spectral palette
#'   (see [RColorBrewer::brewer.pal()]) before interpolation.
#' @param effect_label Label of the x axis (the summary measure).
#' @param ttl Plot title.
#' @param hjust_ttl Horizontal adjustment of the title.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' \dontrun{
#' # b_meta <- brms::brm(
#' #   yi | se(sei) ~ 1 + (1 | Author / Study),
#' #   data = meta_data
#' # )
#' bayesian_forest(b_meta)
#' }
#' @export
bayesian_forest <- function(b_meta_analysis,
                            author_var = "Author",
                            study_var = "Study",
                            palette = 11,
                            effect_label = "Standardized Mean Difference",
                            ttl = "",
                            hjust_ttl = 1) {
  check_brmsfit(b_meta_analysis, arg = "b_meta_analysis")

  r_author <- paste0("r_", author_var)
  r_nested <- paste0("r_", author_var, ":", study_var)
  pars <- tidybayes::get_variables(b_meta_analysis)
  if (!any(startsWith(pars, paste0(r_author, "["))) ||
    !any(startsWith(pars, paste0(r_nested, "[")))) {
    rlang::abort(
      sprintf(
        paste0(
          "`b_meta_analysis` must contain nested grouping terms `(1 | %s/%s)` ",
          "(random effects \"%s\" and \"%s\")."
        ),
        author_var, study_var, r_author, r_nested
      )
    )
  }

  meta_draws <- tidybayes::spread_draws(
    b_meta_analysis, (!!dplyr::sym(r_author))[Author, ], b_Intercept
  )

  meta_draws_nested <- tidybayes::spread_draws(
    b_meta_analysis, (!!dplyr::sym(r_nested))[Author_id, ]
  ) |>
    tidyr::separate(Author_id, into = c("Author", "Study"), sep = "_")

  meta_combined <- meta_draws_nested |>
    dplyr::left_join(meta_draws, by = c(".draw", "Author")) |>
    dplyr::mutate(
      b_Intercept = b_Intercept + !!dplyr::sym(r_author) + !!dplyr::sym(r_nested)
    )

  pooled_effect_draws <- tidybayes::spread_draws(b_meta_analysis, b_Intercept) |>
    dplyr::mutate(Author = "Pooled Effect")

  forest_data <- dplyr::bind_rows(meta_combined, pooled_effect_draws) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Author = stringr::str_squish(stringr::str_replace_all(Author, "[.]", " ")),
      Study = stringr::str_squish(stringr::str_replace_all(Study, "[.]", " "))
    ) |>
    tidyr::unite("Study", Author:Study, sep = " - ", remove = FALSE) |>
    dplyr::mutate(Study = stats::reorder(Study, b_Intercept)) |>
    dplyr::mutate(Study = factor(dplyr::case_match(
      Study, "Pooled Effect - NA" ~ "Pooled Effect",
      .default = Study
    )))

  forest_data_summary <- dplyr::group_by(forest_data, Study) |>
    tidybayes::mean_qi(b_Intercept) |>
    dplyr::arrange(b_Intercept)

  color_vec <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(palette, "Spectral")
  )(length(unique(forest_data$Author)))

  fixef_summary <- brms::fixef(b_meta_analysis)

  ggplot2::ggplot(
    ggplot2::aes(
      x = b_Intercept,
      y = stats::relevel(Study, "Pooled Effect", after = Inf)
    ),
    data = forest_data
  ) +
    # Pooled effect and its interval
    ggplot2::geom_vline(
      xintercept = fixef_summary[1, 1],
      color = "grey", linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = fixef_summary[1, 3:4],
      color = "grey", linetype = 2
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    # Per-study posterior densities
    ggridges::geom_density_ridges(
      ggplot2::aes(fill = Author),
      rel_min_height = 0.01,
      col = "gray40", scale = 1,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    tidybayes::geom_pointinterval(
      ggplot2::aes(xmin = .lower, xmax = .upper),
      data = forest_data_summary, size = 1
    ) +
    ggplot2::scale_fill_manual(values = color_vec) +
    # Point estimates and intervals as text
    ggplot2::geom_text(
      data = forest_data_summary,
      ggplot2::aes(
        label = sprintf("%.2f [%.2f, %.2f]", b_Intercept, .lower, .upper),
        x = Inf
      ),
      hjust = "inward", size = 3.5
    ) +
    ggplot2::labs(
      x = effect_label,
      y = NULL,
      title = ttl,
      subtitle = "Studies are color-coded"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, family = "serif", hjust = hjust_ttl),
      plot.subtitle = ggplot2::element_text(size = 11, family = "serif", hjust = 0.32)
    )
}
