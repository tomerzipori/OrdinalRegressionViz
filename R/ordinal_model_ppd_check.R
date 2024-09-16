#' Visual posterior predictive check for Bayesian ordinal probit models. A barplot
#'
#' This function creates a grouped barplot to check posterior fit to data.
#' @param b_model Bayesian ordinal probit regression model - a brms object
#' @param vars number of variables, either 2 or 3
#' @param response variable name of the response
#' @param ttl plot's title
#' @param filename if not NULL, the name of the png file to save the plot as
#' @param path folder path to save the plot in
#' @param width width of the saved plot - in pixels
#' @param height height of the saved plot - in pixels
#' @return a ggplot plot object
#' @export
ordinal_model_ppd_check <- function(b_model,
                                    vars,
                                    response = "value",
                                    ttl = "",
                                    filename = NULL,
                                    path = getwd(),
                                    width = 2450,
                                    height = 1446) {

  model_data <- insight::get_data(b_model)

  yrep <- brms::posterior_predict(b_model, ndraws = 40, newdata = model_data, re_formula = NULL, allow_new_levels = T)

  y <- as.vector(as.numeric(unlist(model_data[,response])))

  true_vars <- vars + 1

  g <- model_data[,2:true_vars] |>
    purrr::map_dfr(~stringr::str_to_title(.))

  out_plot <- bayesplot::ppc_bars_grouped(y, yrep, interaction(g, sep = ": "), facet_args = list(nrow = 4, ncol = 2), freq = F) +
    ggplot2::scale_x_continuous(breaks = min(y):max(y)) +
    ggplot2::labs(title = ttl) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"))

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

