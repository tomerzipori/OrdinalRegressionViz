#' Forest plot for Bayesian random effects meta-analysis
#'
#' This function makes a forest plot for Bayesian random effects meta-analyses. Code is heavily inspired by Harrer et al. online guide - https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/.
#' Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021). Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London: Chapman & Hall/CRC Press. ISBN 978-0-367-61007-4.
#'
#' Important - this function requires a 'Study' variable for individual studies and an 'Author' clustering variable of studies.
#' Important - In the final plot studies are colored by their author.
#' @param b_meta_analysis Bayesian random effects meta-analysis - a brms object
#' @param palette integer representing a palette in the scale_color_brewer function
#' @param ttl plot's title
#' @param hjust_ttl constant offset of plot title to keep in in the center
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
bayesian_forest <- function(b_meta_analysis,
                            palette = 11,
                            ttl = "",
                            hjust_ttl = 1,
                            filename = NULL,
                            path = getwd(),
                            width = 2450,
                            height = 1446) {

  meta_draws <- tidybayes::spread_draws(b_meta_analysis, r_Author[Author,], b_Intercept)

  meta_draws_nested <- tidybayes::spread_draws(b_meta_analysis, (!!dplyr::sym('r_Author:Study'))[Author_id,]) |>
    tidyr::separate(Author_id, into = c('Author', 'Study'), sep = "_")

  meta_combined <- meta_draws_nested |>
    dplyr::left_join(meta_draws, by = c('.draw', 'Author')) |>
    dplyr::mutate(b_Intercept = b_Intercept + r_Author + !!dplyr::sym('r_Author:Study'))

  pooled_effect_draws <- tidybayes::spread_draws(b_meta_analysis, b_Intercept) |>
    dplyr::mutate(Author = "Pooled Effect")

  forest_data <- dplyr::bind_rows(meta_combined, pooled_effect_draws) |>
    dplyr::ungroup() |>
    dplyr::mutate(Author = stringr::str_squish(stringr::str_replace_all(Author, "[.]", " ")),
           Study = stringr::str_squish(stringr::str_replace_all(Study, "[.]", " "))) |>
    tidyr::unite("Study", Author:Study, sep = " - ", remove = FALSE) |>
    dplyr::mutate(Study = reorder(Study, b_Intercept)) |>
    dplyr::mutate(Study = factor(dplyr::case_match(Study, "Pooled Effect - NA" ~ "Pooled Effect",
                                     .default = Study)))

  forest_data_summary <- dplyr::group_by(forest_data, Study) |>
    tidybayes::mean_qi(b_Intercept) |>
    dplyr::arrange(b_Intercept)

  color_vec <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(palette, "Spectral"))(length(unique(forest_data$Author)))

  out_plot <- ggplot2::ggplot(ggplot2::aes(x = b_Intercept, y = relevel(Study, "Pooled Effect", after = Inf)),
                     data = forest_data) +

    # Add vertical lines for pooled effect and CI
    ggplot2::geom_vline(xintercept = brms::fixef(b_meta_analysis)[1, 1],
               color = "grey", linewidth = 1) +
    ggplot2::geom_vline(xintercept = brms::fixef(b_meta_analysis)[1, 3:4],
               color = "grey", linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, color = "black",
               linewidth = 1) +

    # Add densities
    ggridges::geom_density_ridges(ggplot2::aes(fill = Author),
                                  rel_min_height = 0.01,
                                  col = "gray40", scale = 1,
                                  alpha = 0.8,
                                  show.legend = F) +
    tidybayes::geom_pointinterval(ggplot2::aes(xmin = .lower, xmax = .upper), data = forest_data_summary,
                       size = 1) +

    ggplot2::scale_fill_manual(values = color_vec) +

    # Add text and labels
    ggplot2::geom_text(data = dplyr::mutate_if(forest_data_summary,
                               is.numeric, round, 2),
              ggplot2::aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"),
                  x = Inf), hjust = "inward", size = 3.5) +
    ggplot2::labs(x = "Standardized Mean Difference", # summary measure
         y = ggplot2::element_blank(),
         title = ttl,
         subtitle = "Studies are color-coded") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, family = "serif", hjust = hjust_ttl),
          plot.subtitle = ggplot2::element_text(size = 11, family = "serif", hjust = 0.32))

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

