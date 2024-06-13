#' ROC Curve for two or three 2-level categorical predictor Bayesian ordinal probit regression models
#'
#' This function makes a ROC curve for two or three categorical variables Bayesian ordinal probit regression models.
#' @param b_model Bayesian ordinal probit regression model - a brms object
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param var_facet variable to facet by
#' @param response variable name of the response
#' @param CI credible interval width between 0 and 1
#' @param centrality centrality measure for posterior, either "mean" or "median"
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param ttl plot's title
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
bayesian_roc_plot <- function(b_model,
                              var_signal = "target",
                              var_group = "time",
                              var_facet = NULL,
                              response = "value",
                              CI = 0.95,
                              centrality = "mean",
                              palette = 7,
                              ttl = "",
                              filename = NULL,
                              path = getwd(),
                              width = 2450,
                              height = 1446) {

  grid <- unique(b_model$data[,c(var_signal, var_group, var_facet)])
  probs_rvar <- brms::posterior_epred(b_model, newdata = grid, re_formula = NA) |> posterior::rvar() # re_formula = NA - pop. level predictions
  probs_rvar_df <- as.data.frame(probs_rvar)
  grid <- cbind(grid, probs_rvar_df)

  if (is.null(var_facet)) {
    out_plot <- bayesian_roc_ggplot_2_vars(grid, var_signal = var_signal, var_group = var_group, response = response, CI = CI, centrality = centrality, palette = palette, ttl = ttl)
  }

  if (!is.null(var_facet)) {
    grid1 <- dplyr::filter(grid, !!dplyr::sym(var_facet) == levels(unique(grid[,var_facet]))[1]) |> dplyr::select(-!!dplyr::sym(var_facet))
    plot1 <- bayesian_roc_ggplot_2_vars(grid1, var_signal = var_signal, var_group = var_group, response = response, CI = CI, centrality = centrality, palette = palette, ttl = stringr::str_to_title(levels(unique(grid[,var_facet]))[1]))

    grid2 <- dplyr::filter(grid, !!dplyr::sym(var_facet) == levels(unique(grid[,var_facet]))[2]) |> dplyr::select(-!!dplyr::sym(var_facet))
    plot2 <- bayesian_roc_ggplot_2_vars(grid2, var_signal = var_signal, var_group = var_group, response = response, CI = CI, centrality = centrality, palette = palette, ttl = stringr::str_to_title(levels(unique(grid[,var_facet]))[2]))

    out_plot <- (plot1 + plot2) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(title = ttl, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5, family = "serif")))
  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

