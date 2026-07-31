#' Cumulative sum for posterior::rvar objects
#'
#' Computes the running sum of a [posterior::rvar] vector, preserving the
#' full posterior distribution of each partial sum.
#'
#' @param rvars A [posterior::rvar] vector.
#' @return An rvar vector of the same length holding the cumulative sums.
#' @keywords internal
rvar_cumsum <- function(rvars) {
  out <- rvars
  for (i in seq_along(rvars)[-1]) {
    out[i] <- out[i - 1] + rvars[i]
  }
  out
}


#' Credible/confidence band polygon for ROC curves
#'
#' Builds the polygon that connects the lower and upper interval bounds of a
#' ROC curve so it can be drawn with [ggplot2::geom_polygon()].
#'
#' @param data Data frame holding the curve and its interval bounds.
#' @param x,xmin,xmax Names (strings) of the x coordinate and its bounds.
#' @param y,ymin,ymax Names (strings) of the y coordinate and its bounds.
#' @param group Optional character vector of grouping column names.
#' @return A data frame of polygon vertices with columns `x` and `y` (plus
#'   any grouping columns).
#' @keywords internal
bayesian_CI_bands <- function(data,
                              x, xmin, xmax,
                              y, ymin, ymax,
                              group = NULL) {
  bands_low <- data |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(group))) |>
    dplyr::reframe(
      x = c(dplyr::pick(dplyr::all_of(xmin))[[1]], dplyr::pick(dplyr::all_of(x))[[1]]),
      y = c(dplyr::pick(dplyr::all_of(y))[[1]], dplyr::pick(dplyr::all_of(ymin))[[1]]),
    ) |>
    dplyr::mutate(
      x = dplyr::coalesce(x, y),
      y = dplyr::coalesce(y, x)
    ) |>
    dplyr::arrange(x) |>
    dplyr::distinct()

  bands_high <- data |>
    dplyr::group_by(dplyr::pick(dplyr::any_of(group))) |>
    dplyr::reframe(
      x = c(dplyr::pick(dplyr::all_of(xmax))[[1]], dplyr::pick(dplyr::all_of(x))[[1]]),
      y = c(dplyr::pick(dplyr::all_of(y))[[1]], dplyr::pick(dplyr::all_of(ymax))[[1]]),
    ) |>
    dplyr::mutate(
      x = dplyr::coalesce(x, y),
      y = dplyr::coalesce(y, x)
    ) |>
    dplyr::arrange(dplyr::desc(x)) |>
    dplyr::distinct()

  dplyr::bind_rows(bands_low, bands_high) |>
    dplyr::ungroup()
}


# "1 | 2", "2 | 3", ... labels for the response thresholds, derived from the
# observed response levels.
threshold_labels <- function(response) {
  lev <- if (is.factor(response)) levels(response) else sort(unique(response))
  paste(utils::head(lev, -1), "|", utils::tail(lev, -1))
}

# The response variable of a model frame (its first column, by convention of
# both ordinal::clm and brms::brm model data).
response_name <- function(model_data) {
  names(model_data)[1]
}

# The full set of response categories of a (factor or integer) response.
full_response_levels <- function(x) {
  if (is.factor(x)) levels(x) else sort(unique(x))
}

# Serif title/subtitle styling shared by all plots in the package.
serif_titles <- function(title_size = 20, subtitle_size = 15) {
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = title_size, family = "serif", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = subtitle_size, family = "serif", hjust = 0.5)
  )
}


# ---- Coefficient matching -------------------------------------------------
#
# brms names fixed-effect draws "b_<term>" and distributional (disc) draws
# "b_disc_<term>", where interactions join level-qualified terms with ":".
# These helpers locate coefficients by the *set* of terms they involve, so
# they are robust to the order in which terms appear in the model formula.

# All non-empty subsets of a character vector, as a list.
term_subsets <- function(terms) {
  if (length(terms) == 0) {
    return(list())
  }
  out <- list()
  for (k in seq_along(terms)) {
    out <- c(out, utils::combn(terms, k, simplify = FALSE))
  }
  out
}

# Find, among `coef_names`, the single coefficient built from exactly
# `terms` (a character vector of level-qualified names such as
# "targetfake"), comparing term *sets* so formula order does not matter.
# Returns NA when the model does not contain that term (e.g. an additive
# model without interaction), unless `required` is TRUE.
match_term <- function(coef_names, terms, required = FALSE,
                       call = rlang::caller_env()) {
  parts <- strsplit(coef_names, ":", fixed = TRUE)
  hit <- coef_names[vapply(
    parts,
    function(p) length(p) == length(terms) && setequal(p, terms),
    logical(1)
  )]
  if (length(hit) == 1) {
    return(hit)
  }
  if (length(hit) > 1) {
    rlang::abort(
      sprintf(
        "Multiple coefficients match the term%s %s: %s.",
        if (length(terms) > 1) "s" else "",
        paste0("\"", terms, "\"", collapse = ", "),
        paste0("\"", hit, "\"", collapse = ", ")
      ),
      call = call
    )
  }
  if (required) {
    rlang::abort(
      sprintf(
        "No coefficient for term%s %s was found in the model. Available coefficients: %s.",
        if (length(terms) > 1) "s" else "",
        paste0("\"", terms, "\"", collapse = ", "),
        paste0("\"", coef_names, "\"", collapse = ", ")
      ),
      call = call
    )
  }
  NA_character_
}

# match_term() for brms draws: fixed-effect draws are named "b_<term>" and
# distributional (disc) draws "b_disc_<term>".
match_coef <- function(pars, terms, disc = FALSE, required = FALSE,
                       call = rlang::caller_env()) {
  prefix <- if (disc) "b_disc_" else "b_"
  cand <- pars[startsWith(pars, prefix)]
  if (!disc) {
    cand <- cand[!startsWith(cand, "b_disc_") & !startsWith(cand, "b_Intercept")]
  } else {
    cand <- cand[cand != "b_disc_Intercept"]
  }
  stripped <- substring(cand, nchar(prefix) + 1)
  hit <- match_term(stripped, terms, required = required, call = call)
  if (is.na(hit)) {
    return(NA_character_)
  }
  paste0(prefix, hit)
}

# Sum of the values of the coefficients matching each term set in
# `term_sets`, treating absent terms as zero. Used for clm/clmm fixed
# effects (named numeric vectors).
sum_matched_terms <- function(coefs, coef_names, term_sets) {
  total <- 0
  for (tt in term_sets) {
    nm <- match_term(coef_names, tt)
    if (!is.na(nm)) {
      total <- total + coefs[[nm]]
    }
  }
  total
}

# Resolve a list of term sets to existing coefficient names, dropping terms
# the model does not include.
match_coefs <- function(pars, term_sets, disc = FALSE) {
  if (length(term_sets) == 0) {
    return(character(0))
  }
  out <- vapply(term_sets, function(tt) match_coef(pars, tt, disc = disc), character(1))
  out[!is.na(out)]
}


# ---- Latent distributions implied by cumulative links ---------------------
#
# A cumulative model with link F^{-1} describes a latent variable whose error
# follows the distribution with CDF F: normal for probit, logistic for logit,
# Cauchy for cauchit, and Gumbel (minimum/maximum) for cloglog/loglog. Each
# density below is parameterized by location and scale.

# The latent density function implied by a link, or an error for links with
# no latent-variable representation supported here.
link_density <- function(link, call = rlang::caller_env()) {
  switch(link,
    probit = ,
    probit_approx = function(x, location = 0, scale = 1) {
      stats::dnorm(x, mean = location, sd = scale)
    },
    logit = function(x, location = 0, scale = 1) {
      stats::dlogis(x, location = location, scale = scale)
    },
    cauchit = function(x, location = 0, scale = 1) {
      stats::dcauchy(x, location = location, scale = scale)
    },
    cloglog = function(x, location = 0, scale = 1) {
      z <- (x - location) / scale
      exp(z - exp(z)) / scale
    },
    loglog = function(x, location = 0, scale = 1) {
      z <- -(x - location) / scale
      exp(z - exp(z)) / scale
    },
    rlang::abort(
      sprintf(
        paste0(
          "Link \"%s\" has no supported latent-distribution representation. ",
          "Supported links: \"probit\", \"logit\", \"cauchit\", \"cloglog\", \"loglog\"."
        ),
        link
      ),
      call = call
    )
  )
}

# Human-readable name of the latent distribution implied by a link.
link_distribution_name <- function(link) {
  switch(link,
    probit = ,
    probit_approx = "normal",
    logit = "logistic",
    cauchit = "Cauchy",
    cloglog = "Gumbel (minimum)",
    loglog = "Gumbel (maximum)",
    link
  )
}

# The latent CDF implied by a link, parameterized like link_density().
link_cdf <- function(link, call = rlang::caller_env()) {
  switch(link,
    probit = ,
    probit_approx = function(x, location = 0, scale = 1) {
      stats::pnorm(x, mean = location, sd = scale)
    },
    logit = function(x, location = 0, scale = 1) {
      stats::plogis(x, location = location, scale = scale)
    },
    cauchit = function(x, location = 0, scale = 1) {
      stats::pcauchy(x, location = location, scale = scale)
    },
    cloglog = function(x, location = 0, scale = 1) {
      1 - exp(-exp((x - location) / scale))
    },
    loglog = function(x, location = 0, scale = 1) {
      exp(-exp(-(x - location) / scale))
    },
    rlang::abort(
      sprintf(
        paste0(
          "Link \"%s\" has no supported latent-distribution representation. ",
          "Supported links: \"probit\", \"logit\", \"cauchit\", \"cloglog\", \"loglog\"."
        ),
        link
      ),
      call = call
    )
  )
}

# Area under the ROC curve implied by two latent distributions:
# P(signal draw > noise draw). Closed form for probit; numeric integration
# of f_signal(x) * F_noise(x) otherwise. All arguments are vectorized.
latent_auc <- function(link, mu_noise, scale_noise, mu_signal, scale_signal) {
  if (link %in% c("probit", "probit_approx")) {
    return(stats::pnorm(
      (mu_signal - mu_noise) / sqrt(scale_noise^2 + scale_signal^2)
    ))
  }
  dens <- link_density(link)
  cdf <- link_cdf(link)
  n <- max(length(mu_noise), length(mu_signal))
  mu_noise <- rep_len(mu_noise, n)
  scale_noise <- rep_len(scale_noise, n)
  mu_signal <- rep_len(mu_signal, n)
  scale_signal <- rep_len(scale_signal, n)
  vapply(seq_len(n), function(i) {
    stats::integrate(
      function(x) dens(x, mu_signal[i], scale_signal[i]) * cdf(x, mu_noise[i], scale_noise[i]),
      lower = -Inf, upper = Inf
    )$value
  }, numeric(1))
}


# ---- Design-cell term sets -------------------------------------------------
#
# For one design cell (a combination of the levels of var_group/var_facet,
# encoded by which non-reference terms are "on"), the sets of model terms
# entering each SDT quantity:
# * mean_signal: terms shifting the signal distribution's location (the
#   signal main effect and all its interactions with active terms),
# * disc_noise / disc_signal: terms entering the disc (log inverse scale)
#   of each distribution,
# * shifts: terms shifting the response thresholds (criteria).
cell_specs <- function(term_sig, on_terms) {
  list(
    mean_signal = lapply(
      c(list(character(0)), term_subsets(on_terms)),
      function(s) c(term_sig, s)
    ),
    disc_noise = term_subsets(on_terms),
    disc_signal = term_subsets(c(term_sig, on_terms)),
    shifts = term_subsets(on_terms)
  )
}


# ---- Empirical ROC points --------------------------------------------------

# Observed cumulative response proportions per group, in the same
# orientation as the model-based ROC: at threshold k, the hit rate is
# P(response <= k | first level of var_signal) and the false-alarm rate is
# P(response <= k | second level of var_signal).
#
# `response_levels` is the full set of response categories; pass it whenever
# `data` may be a subset that does not realize every category (e.g. one
# facet's rows), otherwise the thresholds would silently misalign.
empirical_roc_points <- function(data, response, var_signal, var_group,
                                 response_levels = NULL) {
  resp_raw <- data[[response]]
  if (is.null(response_levels)) {
    response_levels <- if (is.factor(resp_raw)) {
      levels(resp_raw)
    } else {
      sort(unique(resp_raw))
    }
  }
  resp <- as.integer(factor(resp_raw, levels = response_levels))
  n_cat <- length(response_levels)
  sig_levels <- levels(data[[var_signal]])
  group_levels <- levels(data[[var_group]])

  out <- lapply(group_levels, function(g) {
    in_group <- data[[var_group]] == g
    cumprop <- function(sig_level) {
      r <- resp[in_group & data[[var_signal]] == sig_level]
      (cumsum(tabulate(r, nbins = n_cat)) / length(r))[-n_cat]
    }
    data.frame(
      group = g,
      Sensitivity = cumprop(sig_levels[1]),
      FAR = cumprop(sig_levels[2])
    )
  })
  out <- do.call(rbind, out)
  out$group <- factor(out$group, levels = group_levels)
  names(out)[1] <- var_group
  out
}


# ---- Posterior latent-distribution machinery ------------------------------

# Posterior density band of a latent distribution whose location is the sum
# of the draws of `mean_coefs` and whose scale is exp(-sum of `disc_coefs`
# draws) (i.e. the product of 1/exp() of each disc coefficient). `dens` is a
# density function of the form function(x, location, scale), see
# link_density(). Returns a data frame with columns x, d, .lower, .upper via
# tidybayes::curve_interval().
latent_density_band <- function(b_model, mean_coefs, disc_coefs, x,
                                width = 0.9, ndraws = NULL,
                                dens = link_density("probit")) {
  vars <- c(mean_coefs, disc_coefs)

  if (length(vars) == 0) {
    # Fully constant cell (e.g. reference distribution of an equal-variance
    # model): a standard latent distribution with no posterior uncertainty.
    d <- dens(x)
    return(dplyr::tibble(x = x, d = d, .lower = d, .upper = d))
  }

  draws <- posterior::as_draws_df(b_model, variable = vars)
  if (!is.null(ndraws) && ndraws < posterior::ndraws(draws)) {
    keep <- sort(sample.int(posterior::ndraws(draws), ndraws))
    # subset_draws() messages about merging chains; that is expected here.
    draws <- suppressMessages(posterior::subset_draws(draws, draw = keep))
  }
  draws_df <- as.data.frame(draws)

  mu <- if (length(mean_coefs)) {
    rowSums(draws_df[, mean_coefs, drop = FALSE])
  } else {
    rep(0, nrow(draws_df))
  }
  sigma <- if (length(disc_coefs)) {
    exp(-rowSums(draws_df[, disc_coefs, drop = FALSE]))
  } else {
    rep(1, nrow(draws_df))
  }

  dmat <- mapply(function(m, s) dens(x, location = m, scale = s), mu, sigma)
  long <- dplyr::tibble(
    .draw = rep(seq_along(mu), each = length(x)),
    x = rep(x, times = length(mu)),
    d = as.vector(dmat)
  )
  tidybayes::curve_interval(long, .along = x, .width = width)
}

# Posterior thresholds (criteria) of a cumulative model as an rvar vector,
# shifted by the coefficients in `shift_coefs` (criterion shifts appear in
# the mean part of the model).
threshold_rvars <- function(b_model, shift_coefs) {
  vars <- c("b_Intercept", shift_coefs)
  drv <- posterior::as_draws_rvars(b_model, variable = vars)
  th <- drv$b_Intercept
  for (cf in shift_coefs) {
    th <- th - drv[[cf]]
  }
  th
}
