#' Input validation helpers
#'
#' Small internal checks that fail early with informative error messages.
#'
#' @name orv-checks
#' @keywords internal
NULL

# Assert that `model` is a brmsfit and that brms is available.
check_brmsfit <- function(model, arg = "b_model", call = rlang::caller_env()) {
  if (!inherits(model, "brmsfit")) {
    rlang::abort(
      sprintf(
        "`%s` must be a <brmsfit> object (fitted with brms::brm()), not a <%s>.",
        arg, paste(class(model), collapse = "/")
      ),
      call = call
    )
  }
  rlang::check_installed(
    "brms",
    reason = sprintf("to work with the <brmsfit> object passed to `%s`.", arg)
  )
  invisible(model)
}

# Assert that `model` is an ordinal::clm or ordinal::clmm fit.
check_clm <- function(model, arg = "model", call = rlang::caller_env()) {
  if (!inherits(model, "clm") && !inherits(model, "clmm")) {
    rlang::abort(
      sprintf(
        "`%s` must be a <clm> or <clmm> object (fitted with the ordinal package), not a <%s>.",
        arg, paste(class(model), collapse = "/")
      ),
      call = call
    )
  }
  if (is.null(model$model)) {
    rlang::abort(
      sprintf(
        "`%s` does not contain its model frame. Refit with `model = TRUE` (the default in ordinal::clm()).",
        arg
      ),
      call = call
    )
  }
  invisible(model)
}

# Assert that `var` names a factor column of `data` with exactly `n_levels` levels.
check_factor_var <- function(data, var, arg, n_levels = 2, call = rlang::caller_env()) {
  if (!rlang::is_string(var)) {
    rlang::abort(sprintf("`%s` must be a single variable name (a string).", arg), call = call)
  }
  if (!var %in% names(data)) {
    rlang::abort(
      sprintf(
        "`%s` (\"%s\") is not a variable of the model data. Available variables: %s.",
        arg, var, paste0("\"", names(data), "\"", collapse = ", ")
      ),
      call = call
    )
  }
  v <- data[[var]]
  if (!is.factor(v)) {
    rlang::abort(
      sprintf("`%s` (\"%s\") must be a factor in the model data, not <%s>.", arg, var, class(v)[1]),
      call = call
    )
  }
  if (nlevels(droplevels(v)) != n_levels) {
    rlang::abort(
      sprintf(
        "`%s` (\"%s\") must have exactly %d levels, but has %d (%s).",
        arg, var, n_levels, nlevels(droplevels(v)),
        paste0("\"", levels(droplevels(v)), "\"", collapse = ", ")
      ),
      call = call
    )
  }
  invisible(v)
}

# Assert that `x` is a single probability strictly inside (0, 1).
check_prob <- function(x, arg, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0 || x >= 1) {
    rlang::abort(
      sprintf("`%s` must be a single number strictly between 0 and 1.", arg),
      call = call
    )
  }
  invisible(x)
}

# Assert that a brmsfit is a cumulative-family model.
check_cumulative <- function(model, arg = "b_model", call = rlang::caller_env()) {
  fam <- model$family
  if (!is.null(fam$family) && fam$family != "cumulative") {
    rlang::abort(
      sprintf(
        "`%s` must be a cumulative (ordinal) regression model, but its family is \"%s\".",
        arg, fam$family
      ),
      call = call
    )
  }
  invisible(model)
}

# Optional length-2 character labels (e.g. legend labels); falls back to
# str_to_title() of the factor levels when NULL.
check_labels2 <- function(labels, levels, arg, call = rlang::caller_env()) {
  if (is.null(labels)) {
    return(stringr::str_to_title(levels))
  }
  if (!is.character(labels) || length(labels) != 2) {
    rlang::abort(sprintf("`%s` must be NULL or a character vector of length 2.", arg), call = call)
  }
  labels
}
