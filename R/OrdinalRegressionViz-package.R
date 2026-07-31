#' @keywords internal
#' @importFrom patchwork plot_layout
"_PACKAGE"
# The importFrom above is required (not just cosmetic): combining panels with
# patchwork's `/`, `|` and `+` operators only works once the patchwork
# namespace is loaded, which the import guarantees at package load time.

# Column names used inside dplyr/ggplot2 non-standard evaluation.
utils::globalVariables(c(
  ".draw", ".lower", ".upper", ".value", ".variable",
  "Author", "Author_id", "FAR", "FAR_high", "FAR_low",
  "Response", "Sensitivity", "Sensitivity_high", "Sensitivity_low",
  "Specificity", "Specificity_high", "Specificity_low",
  "Study", "Threshold", "b_Intercept", "asymp.LCL", "asymp.UCL",
  "cumprob", "cut", "cut_label", "d", "exc.prob", "prob", "x", "y"
))
