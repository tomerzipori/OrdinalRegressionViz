test_that("rvar_cumsum matches cumsum draw by draw", {
  draws <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2) # 2 draws x 3 variables
  rv <- posterior::rvar(draws)
  out <- rvar_cumsum(rv)
  expected <- t(apply(draws, 1, cumsum))
  expect_equal(posterior::draws_of(out), posterior::draws_of(posterior::rvar(expected)))
})

test_that("rvar_cumsum leaves length-1 input unchanged", {
  rv <- posterior::rvar(matrix(c(1.5, 2.5), nrow = 2))
  expect_equal(rvar_cumsum(rv), rv)
})

test_that("threshold_labels derives labels from factors and numerics", {
  expect_equal(
    threshold_labels(factor(c("1", "2", "3"), levels = c("1", "2", "3"))),
    c("1 | 2", "2 | 3")
  )
  expect_equal(threshold_labels(c(3, 1, 2, 1)), c("1 | 2", "2 | 3"))
})

test_that("term_subsets returns all non-empty subsets", {
  expect_equal(term_subsets(character(0)), list())
  expect_equal(term_subsets("a"), list("a"))
  out <- term_subsets(c("a", "b"))
  expect_length(out, 3)
  expect_true(list(c("a", "b")) %in% out || any(vapply(out, identical, logical(1), y = c("a", "b"))))
})

test_that("match_term matches by term set, order-independently", {
  coefs <- c("targetfake", "timepost", "timepost:targetfake")
  expect_equal(match_term(coefs, "targetfake"), "targetfake")
  expect_equal(match_term(coefs, c("targetfake", "timepost")), "timepost:targetfake")
  expect_true(is.na(match_term(coefs, "conditionb")))
  expect_error(match_term(coefs, "conditionb", required = TRUE), "No coefficient")
})

test_that("match_coef distinguishes mean and disc coefficients", {
  pars <- c(
    "b_Intercept[1]", "b_Intercept[2]",
    "b_targetfake", "b_timepost", "b_timepost:targetfake",
    "b_disc_Intercept", "b_disc_targetfake", "b_disc_timepost:targetfake"
  )
  expect_equal(match_coef(pars, "targetfake"), "b_targetfake")
  expect_equal(match_coef(pars, "targetfake", disc = TRUE), "b_disc_targetfake")
  expect_equal(
    match_coef(pars, c("targetfake", "timepost"), disc = TRUE),
    "b_disc_timepost:targetfake"
  )
  expect_true(is.na(match_coef(pars, "timepost", disc = TRUE)))
})

test_that("sum_matched_terms treats absent terms as zero", {
  coefs <- c(targetfake = 1.5, timepost = -0.5)
  nms <- names(coefs)
  expect_equal(sum_matched_terms(coefs, nms, list("targetfake")), 1.5)
  expect_equal(
    sum_matched_terms(coefs, nms, list("targetfake", c("targetfake", "timepost"))),
    1.5
  )
  expect_equal(sum_matched_terms(coefs, nms, list()), 0)
})

test_that("link_density returns correct standard densities", {
  x <- seq(-3, 3, by = 0.5)
  expect_equal(link_density("probit")(x), stats::dnorm(x))
  expect_equal(link_density("probit_approx")(x), stats::dnorm(x))
  expect_equal(link_density("logit")(x), stats::dlogis(x))
  expect_equal(link_density("cauchit")(x), stats::dcauchy(x))
  # Gumbel densities integrate to ~1
  for (link in c("cloglog", "loglog")) {
    dens <- link_density(link)
    expect_equal(
      stats::integrate(dens, -Inf, Inf)$value, 1,
      tolerance = 1e-6
    )
  }
  expect_error(link_density("loggamma"), "no supported latent-distribution")
})

test_that("link_density handles location and scale", {
  x <- c(-1, 0, 2)
  expect_equal(
    link_density("logit")(x, location = 1, scale = 2),
    stats::dlogis(x, location = 1, scale = 2)
  )
  d <- link_density("cloglog")
  expect_equal(d(x, location = 1, scale = 2), d((x - 1) / 2) / 2)
})

test_that("bayesian_CI_bands builds a closed polygon data frame", {
  df <- data.frame(
    g = "a",
    x = c(0.2, 0.5), x_lo = c(0.1, 0.4), x_hi = c(0.3, 0.6),
    y = c(0.4, 0.8), y_lo = c(0.3, 0.7), y_hi = c(0.5, 0.9)
  )
  bands <- bayesian_CI_bands(
    df, "x", "x_lo", "x_hi", "y", "y_lo", "y_hi",
    group = "g"
  )
  expect_true(all(c("g", "x", "y") %in% names(bands)))
  expect_gt(nrow(bands), nrow(df))
  expect_false(anyNA(bands$x))
  expect_false(anyNA(bands$y))
})
