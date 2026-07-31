#' Grouped posterior predictive check for a Bayesian ordinal model
#'
#' Draws a grouped bar plot comparing the observed response distribution
#' with posterior predictive draws (via [bayesplot::ppc_bars_grouped()]),
#' with one facet per combination of the grouping variables.
#'
#' @inheritParams bayesian_SDT_distribution_plot
#' @param group_vars Character vector with the names of the model variables
#'   whose combinations define the facets (e.g.
#'   `c("target", "time")`).
#' @param ndraws Number of posterior predictive draws to use.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' \dontrun{
#' ordinal_model_ppd_check(b_fit, group_vars = c("target", "time"))
#' }
#' @export
ordinal_model_ppd_check <- function(b_model,
                                    group_vars,
                                    ndraws = 40,
                                    ttl = "") {
  check_brmsfit(b_model)
  model_data <- b_model$data

  if (!is.character(group_vars) || length(group_vars) == 0) {
    rlang::abort("`group_vars` must be a character vector of model variable names.")
  }
  missing_vars <- setdiff(group_vars, names(model_data))
  if (length(missing_vars) > 0) {
    rlang::abort(
      sprintf(
        "`group_vars` not found in the model data: %s. Available variables: %s.",
        paste0("\"", missing_vars, "\"", collapse = ", "),
        paste0("\"", names(model_data), "\"", collapse = ", ")
      )
    )
  }

  yrep <- brms::posterior_predict(b_model, ndraws = ndraws)
  y <- as.numeric(model_data[[response_name(model_data)]])

  # Title-case the group labels while preserving the level order.
  pretty_factor <- function(v) {
    v <- as.factor(v)
    factor(
      stringr::str_to_title(as.character(v)),
      levels = stringr::str_to_title(levels(v))
    )
  }
  g <- interaction(lapply(model_data[group_vars], pretty_factor), sep = ": ")

  n_groups <- nlevels(droplevels(g))
  facet_args <- list(ncol = min(2, n_groups), nrow = ceiling(n_groups / 2))

  p <- bayesplot::ppc_bars_grouped(y, yrep, g, facet_args = facet_args, freq = FALSE)
  # ppc_bars_grouped() already sets an x scale; replacing it emits a message.
  suppressMessages(
    p <- p +
      ggplot2::scale_x_continuous(breaks = min(y):max(y)) +
      ggplot2::labs(title = ttl) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "white", color = "white"),
        panel.background = ggplot2::element_rect(fill = "white", color = "white")
      )
  )
  p
}
