#' Cumulative sum for posterior::rvar objects
#'
#' This function is a simple loop to calculate a cumulative sum of rvar objects
#' @param rvars a posterior rvar object
#' @return numeric vector
#' @keywords internal
rvar_cumsum <- function(rvars) {
  out <- rvars
  for (i in c(2:length(rvars))) {
    out[i] <- out[i-1] + rvars[i]
  }
  return(out)
}


#' Calculating credible interval bands for ROC curves
#'
#' This is a helper function for the main Bayesian ROC curve plotting function, calculating the CI bands around the curve
#' @param xmin lower bound of CI on the x axis
#' @param xmax upper bound of CI on the x axis
#' @param ymin lower bound of CI on the y axis
#' @param ymax upper bound of CI on the y axis
#' @return data frame
#' @keywords internal
bayesian_CI_bands <- function(data,
                              x, xmin, xmax,
                              y, ymin, ymax,
                              group = NULL) {
  bands_low <- data |>
    group_by(pick(any_of(group))) |>
    reframe(
      x = c(pick(all_of(xmin))[[1]], pick(all_of(x))[[1]]),
      y = c(pick(all_of(y))[[1]], pick(all_of(ymin))[[1]]),
    ) |>
    mutate(
      x = coalesce(x, y),
      y = coalesce(y, x)
    ) |>
    arrange(x) |>
    distinct()


  bands_high <- data |>
    group_by(pick(any_of(group))) |>
    reframe(
      x = c(pick(all_of(xmax))[[1]], pick(all_of(x))[[1]]),
      y = c(pick(all_of(y))[[1]], pick(all_of(ymax))[[1]]),
    ) |>
    mutate(
      x = coalesce(x, y),
      y = coalesce(y, x)
    ) |>
    arrange(desc(x)) |>
    distinct()

  bands <- bind_rows(bands_low, bands_high) |>
    ungroup()

  return(bands)
}


#' ROC Curve for two 2-level categorical predictor Bayesian ordinal probit regression models
#'
#' This function makes a ROC curve for two categorical variables Bayesian ordinal probit regression models.
#' @param grid data frame of conditional probabilities posteriors
#' @param var_signal variable analogous to the classic SDT signal, e.g. old/new
#' @param var_group variable name of the co-variate, e.g. experimental group
#' @param response variable name of the response
#' @param CI credible interval width between 0 and 1
#' @param centrality centrality measure for posterior, either "mean" or "median"
#' @param palette integer representing a divergent palette in the scale_color_brewer function
#' @param ttl plot's title
#' @return a ggplot plot object
#' @keywords internal
bayesian_roc_ggplot_2_vars <- function(grid,
                                       var_signal = "target",
                                       var_group = "time",
                                       response = "value",
                                       CI = 0.95,
                                       centrality = "mean",
                                       palette = 7,
                                       ttl = "") {

  grid_long <- grid |>
    pivot_longer(cols = names(grid)[-c(1:2)],
                 names_to = "cut",
                 values_to = "prob") |>
    mutate(cut = as.numeric(cut)) |>
    pivot_wider(values_from = prob, names_from = !!sym(var_signal)) |>
    group_by(!!sym(var_group)) |>
    mutate(
      Sensitivity = lag(rvar_cumsum(!!sym(levels(grid[var_signal][1,])[1])), default = 0),
      Specificity = rev(rvar_cumsum(rev(!!sym(levels(grid[var_signal][1,])[2])))),
      Threshold = paste0(lag(cut), "|", cut)
    )|>
    rows_append(data.frame(Sensitivity = 1, Specificity = 0)) |>
    mutate(
      Threshold = ifelse(mean(Sensitivity) %in% c(0, 1), NA, Threshold)
    )

  grid_long_Sensitivity <- grid_long |>
    group_by(Threshold, !!sym(var_group)) |>
    drop_na() |>
    select(Threshold, !!sym(var_group), Sensitivity) |>
    mutate(Sensitivity_low = quantile(Sensitivity, probs = (1-CI)/2),
           Sensitivity_high = quantile(Sensitivity, probs = (1-CI)/2 + CI),
           Sensitivity = case_when(centrality == "mean" ~ mean(Sensitivity),
                                   centrality == "median" ~ median(Sensitivity)))

  grid_long_Specificity <- grid_long |>
    group_by(Threshold, !!sym(var_group)) |>
    drop_na() |>
    select(Threshold, !!sym(var_group), Specificity) |>
    mutate(Specificity_low = quantile(Specificity, probs = (1-CI)/2),
           Specificity_high = quantile(Specificity, probs = (1-CI)/2 + CI),
           Specificity = case_when(centrality == "mean" ~ mean(Specificity),
                                   centrality == "median" ~ median(Specificity)))

  edges <- data.frame(
    time = rep(levels(unique(grid[var_group])[1,]), each = length(levels(unique(grid[var_group])[1,]))),
    Sensitivity = rep(0:1, times = 2),
    Specificity = rep(1:0, times = 2)
  )

  roc_data_grid <- grid_long_Sensitivity |>
    full_join(grid_long_Specificity, by = join_by(Threshold, !!sym(var_group))) |>
    rows_append(edges) |>
    mutate(FAR = 1 - Specificity,
           FAR_low = 1 - Specificity_low,
           FAR_high = 1 - Specificity_high) |>
    arrange(!!sym(var_group), Sensitivity)

  bands <- roc_data_grid |>
    bayesian_CI_bands("FAR", "FAR_low", "FAR_high",
                      "Sensitivity", "Sensitivity_low", "Sensitivity_high",
                      group = c(var_group))

  ggplot(roc_data_grid, aes(FAR, Sensitivity)) +
    geom_polygon(aes(x, y, fill = !!sym(var_group)), data = bands,
                 alpha = 0.4) +
    scale_fill_brewer(type = "qual", palette = 2) +
    ggnewscale::new_scale_fill() +

    geom_path(aes(linetype = !!sym(var_group)), linewidth = 0.8, show.legend = F) +
    scale_linetype_manual(values = c("dashed", "solid")) +

    geom_linerange(aes(xmin = FAR_low, xmax = FAR_high), color = "grey40") +
    geom_linerange(aes(ymin = Sensitivity_low, ymax = Sensitivity_high), color = "grey40") +
    geom_point(aes(fill = ordered(Threshold)), shape = 21, size = 3,
               data = \(d) drop_na(d, Threshold)) +

    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    expand_limits(x = c(0,1), y = c(0,1)) +
    scale_fill_brewer("Threshold", type = "div", palette = palette,
                      na.translate = FALSE) +
    labs(color = NULL, fill = var_group, y = "HR", title = ttl, subtitle = paste0(as.character(100*CI), "% Credible Interval")) +
    coord_fixed() +
    theme_classic() +
    theme(plot.title = element_text(size = 19, family = "serif", hjust = 0.5),
          plot.subtitle = element_text(size = 15, family = "serif", hjust = 0.5))

}


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
    grid1 <- filter(grid, !!sym(var_facet) == levels(unique(grid[,var_facet]))[1]) |> select(-!!sym(var_facet))
    plot1 <- bayesian_roc_ggplot_2_vars(grid1, var_signal = var_signal, var_group = var_group, response = response, CI = CI, centrality = centrality, palette = palette, ttl = str_to_title(levels(unique(grid[,var_facet]))[1]))

    grid2 <- filter(grid, !!sym(var_facet) == levels(unique(grid[,var_facet]))[2]) |> select(-!!sym(var_facet))
    plot2 <- bayesian_roc_ggplot_2_vars(grid2, var_signal = var_signal, var_group = var_group, response = response, CI = CI, centrality = centrality, palette = palette, ttl = str_to_title(levels(unique(grid[,var_facet]))[2]))

    out_plot <- (plot1 + plot2) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(title = ttl, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5, family = "serif")))
  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, path = path, plot = out_plot, width = width, height = height, units = "px")
  }

  out_plot

}

