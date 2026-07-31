#' Signal detection indices of an ordinal regression model
#'
#' Computes the signal detection theory (SDT) quantities implied by a
#' cumulative regression model, per design cell: the sensitivity index
#' (`d_prime`), the scale ratio of the signal and noise distributions
#' (`sigma_ratio`), the response thresholds (`criterion ...`), and the area
#' under the implied ROC curve (`auc`).
#'
#' Definitions follow the classic SDT conventions, on the latent scale of
#' the model's link function (normal for probit):
#' * `d_prime` is the location difference between the signal (second level
#'   of `var_signal`) and noise (first level) distributions, in units of the
#'   noise distribution's scale.
#' * `sigma_ratio` is the signal scale divided by the noise scale (1 unless
#'   the model includes a `disc` part; `ordinal::clm()` scale models are not
#'   supported).
#' * `criterion <k|k+1>` are the response thresholds of the cell, in raw
#'   latent units (matching [SDT_distributions_plot()]).
#' * `auc` is the probability that a random signal draw exceeds a random
#'   noise draw. It has a closed form under probit and is computed by
#'   numeric integration for other links.
#'
#' For `brmsfit` models the full posterior of every index is used, and the
#' returned `estimate` is the posterior median with a `ci` central credible
#' interval. For `clm`/`clmm` models the indices are point estimates
#' (`lower`/`upper` are `NA`); use the Bayesian interface if you need
#' uncertainty for the derived quantities.
#'
#' @param model A cumulative regression model: `ordinal::clm()`,
#'   `ordinal::clmm()`, or a `brmsfit` from [brms::brm()].
#' @param var_signal Name (string) of the 2-level factor playing the role of
#'   the SDT signal. Its first level is the noise distribution.
#' @param var_group,var_facet Optional names (strings) of additional 2-level
#'   factor covariates; indices are computed for every combination of their
#'   levels.
#' @param ci Width of the credible intervals (`brmsfit` only).
#' @param ndraws Optional number of posterior draws to use (`brmsfit` only);
#'   `NULL` uses all draws.
#' @param ... Passed on to methods.
#' @return A data frame with one row per index and design cell: the levels
#'   of `var_group`/`var_facet` (when supplied), `index`, `estimate`,
#'   `lower`, and `upper`.
#' @seealso [SDT_distributions_plot()] and [bayesian_SDT_distribution_plot()]
#'   for the corresponding latent-distribution displays.
#' @examplesIf requireNamespace("ordinal", quietly = TRUE)
#' fit <- ordinal::clm(value ~ target * time, data = sdt_ratings, link = "probit")
#' sdt_indices(fit, var_signal = "target", var_group = "time")
#' @export
sdt_indices <- function(model, var_signal, var_group = NULL, var_facet = NULL,
                        ci = 0.95, ndraws = NULL, ...) {
  UseMethod("sdt_indices")
}

#' @export
sdt_indices.default <- function(model, var_signal, var_group = NULL,
                                var_facet = NULL, ci = 0.95, ndraws = NULL, ...) {
  rlang::abort(
    sprintf(
      "`model` must be a <clm>, <clmm>, or <brmsfit> object, not a <%s>.",
      paste(class(model), collapse = "/")
    )
  )
}

#' @export
sdt_indices.clm <- function(model, var_signal, var_group = NULL,
                            var_facet = NULL, ci = 0.95, ndraws = NULL, ...) {
  check_clm(model)
  if (length(model$zeta) > 0) {
    rlang::abort(
      "clm models with a scale part are not supported; fit the scale (disc) effects with brms instead."
    )
  }

  model_data <- model$model
  vars <- validate_sdt_vars(model_data, var_signal, var_group, var_facet)
  link <- model$link

  thresholds <- model$alpha
  coefs <- stats::coef(model)
  beta_names <- setdiff(names(coefs), names(thresholds))
  match_term(beta_names, vars$term_sig, required = TRUE)
  crit_labels <- names(thresholds)

  rows <- lapply(seq_len(nrow(vars$cells)), function(i) {
    on_terms <- cell_on_terms(vars, i)
    specs <- cell_specs(vars$term_sig, on_terms)

    mu_signal <- sum_matched_terms(coefs, beta_names, specs$mean_signal)
    shift <- sum_matched_terms(coefs, beta_names, specs$shifts)

    cell_frame(
      vars, i,
      index = c(
        "d_prime", "sigma_ratio", "auc",
        paste("criterion", crit_labels)
      ),
      estimate = c(
        mu_signal,
        1,
        latent_auc(link, 0, 1, mu_signal, 1),
        unname(thresholds) - shift
      ),
      lower = NA_real_,
      upper = NA_real_
    )
  })
  do.call(rbind, rows)
}

#' @export
sdt_indices.clmm <- sdt_indices.clm

#' @export
sdt_indices.brmsfit <- function(model, var_signal, var_group = NULL,
                                var_facet = NULL, ci = 0.95, ndraws = NULL, ...) {
  check_brmsfit(model, arg = "model")
  check_cumulative(model, arg = "model")
  check_prob(ci, "ci")

  model_data <- model$data
  vars <- validate_sdt_vars(model_data, var_signal, var_group, var_facet)
  link <- model$family$link
  link_cdf(link) # validate the link early

  pars <- tidybayes::get_variables(model)
  has_disc <- "b_disc_Intercept" %in% pars
  match_coef(pars, vars$term_sig, required = TRUE)

  draws <- posterior::as_draws_df(model)
  if (!is.null(ndraws) && ndraws < posterior::ndraws(draws)) {
    keep <- sort(sample.int(posterior::ndraws(draws), ndraws))
    draws <- suppressMessages(posterior::subset_draws(draws, draw = keep))
  }
  draws <- as.data.frame(draws)

  intercept_cols <- grep("^b_Intercept\\[", names(draws), value = TRUE)
  crit_labels <- threshold_labels(model_data[[response_name(model_data)]])
  if (length(crit_labels) != length(intercept_cols)) {
    crit_labels <- paste0(seq_along(intercept_cols), "|", seq_along(intercept_cols) + 1)
  }

  col_sums <- function(cols) {
    if (length(cols) == 0) {
      rep(0, nrow(draws))
    } else {
      rowSums(draws[, cols, drop = FALSE])
    }
  }

  summarize_draws_row <- function(x) {
    c(
      estimate = stats::median(x),
      lower = unname(stats::quantile(x, (1 - ci) / 2)),
      upper = unname(stats::quantile(x, (1 + ci) / 2))
    )
  }

  rows <- lapply(seq_len(nrow(vars$cells)), function(i) {
    on_terms <- cell_on_terms(vars, i)
    specs <- cell_specs(vars$term_sig, on_terms)

    mu_signal <- col_sums(match_coefs(pars, specs$mean_signal))
    scale_noise <- if (has_disc) {
      exp(-col_sums(c("b_disc_Intercept", match_coefs(pars, specs$disc_noise, disc = TRUE))))
    } else {
      rep(1, nrow(draws))
    }
    scale_signal <- if (has_disc) {
      exp(-col_sums(c("b_disc_Intercept", match_coefs(pars, specs$disc_signal, disc = TRUE))))
    } else {
      rep(1, nrow(draws))
    }
    shift <- col_sums(match_coefs(pars, specs$shifts))

    index_draws <- c(
      list(
        d_prime = mu_signal / scale_noise,
        sigma_ratio = scale_signal / scale_noise,
        auc = latent_auc(link, 0, scale_noise, mu_signal, scale_signal)
      ),
      stats::setNames(
        lapply(intercept_cols, function(cl) draws[[cl]] - shift),
        paste("criterion", crit_labels)
      )
    )

    summaries <- t(vapply(index_draws, summarize_draws_row, numeric(3)))
    cell_frame(
      vars, i,
      index = names(index_draws),
      estimate = summaries[, "estimate"],
      lower = summaries[, "lower"],
      upper = summaries[, "upper"]
    )
  })
  do.call(rbind, rows)
}


# ---- shared helpers for sdt_indices ---------------------------------------

# Validate the design variables and build the cell grid. Returns a list with
# the level-qualified terms and a `cells` data frame of level indices.
validate_sdt_vars <- function(model_data, var_signal, var_group, var_facet,
                              call = rlang::caller_env()) {
  sig <- check_factor_var(model_data, var_signal, "var_signal", call = call)
  term_grp <- NULL
  term_facet <- NULL
  if (!is.null(var_group)) {
    grp <- check_factor_var(model_data, var_group, "var_group", call = call)
    term_grp <- paste0(var_group, levels(grp)[2])
  }
  if (!is.null(var_facet)) {
    fct <- check_factor_var(model_data, var_facet, "var_facet", call = call)
    term_facet <- paste0(var_facet, levels(fct)[2])
  }

  cells <- expand.grid(
    g = if (is.null(var_group)) NA_integer_ else 1:2,
    f = if (is.null(var_facet)) NA_integer_ else 1:2
  )

  list(
    data = model_data,
    var_signal = var_signal, var_group = var_group, var_facet = var_facet,
    term_sig = paste0(var_signal, levels(sig)[2]),
    term_grp = term_grp, term_facet = term_facet,
    cells = cells
  )
}

# The active (non-reference) terms of cell `i`.
cell_on_terms <- function(vars, i) {
  c(
    if (!is.na(vars$cells$g[i]) && vars$cells$g[i] == 2) vars$term_grp,
    if (!is.na(vars$cells$f[i]) && vars$cells$f[i] == 2) vars$term_facet
  )
}

# A block of result rows for cell `i`, prefixed with the cell's factor
# levels (one column per supplied design variable).
cell_frame <- function(vars, i, index, estimate, lower, upper) {
  out <- data.frame(
    index = index,
    estimate = unname(estimate),
    lower = unname(lower),
    upper = unname(upper)
  )
  if (!is.null(vars$var_facet)) {
    out <- cbind(
      stats::setNames(
        data.frame(levels(vars$data[[vars$var_facet]])[vars$cells$f[i]]),
        vars$var_facet
      ),
      out
    )
  }
  if (!is.null(vars$var_group)) {
    out <- cbind(
      stats::setNames(
        data.frame(levels(vars$data[[vars$var_group]])[vars$cells$g[i]]),
        vars$var_group
      ),
      out
    )
  }
  out
}
