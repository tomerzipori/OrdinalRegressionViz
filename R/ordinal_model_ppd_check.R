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

  yrep <- posterior_predict(b_model, ndraws = 40, newdata = model_data, re_formula = NULL)

  y <- as.vector(as.numeric(model_data[,response]))

  g <- model_data[,c(2:vars+1)]

  out_plot <- ppc_bars_grouped(y, yrep, interaction(g, sep = ": "), facet_args = list(nrow = 4, ncol = 2)) +
    ggplot2::scale_x_continuous(breaks = min(y):max(y)) +
    labs(title = ttl) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white"))

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

