# These tests use small pre-fitted brms fixtures (see
# data-raw/make_test_fixtures.R) and skip when the fixtures or brms are not
# available (e.g. in the built package on CRAN).

test_that("bayesian_roc_plot returns a ggplot for a 2-variable model", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  p <- bayesian_roc_plot(b_fit, var_signal = "target", var_group = "time")
  expect_s3_class(p, "ggplot")
  expect_renders(p)
})

test_that("bayesian_roc_plot facets for a 3-variable model", {
  b_fit <- read_fixture("brms_sdt_3var.rds")
  p <- bayesian_roc_plot(
    b_fit,
    var_signal = "target", var_group = "time", var_facet = "condition"
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("bayesian_roc_plot supports centrality and CI options", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  p <- bayesian_roc_plot(
    b_fit,
    var_signal = "target", var_group = "time",
    CI = 0.8, centrality = "median"
  )
  expect_s3_class(p, "ggplot")
  expect_error(
    bayesian_roc_plot(b_fit, "target", "time", centrality = "mode"),
    "should be one of"
  )
})

test_that("bayesian_SDT_distribution_plot works for unequal-variance models", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  p <- bayesian_SDT_distribution_plot(
    b_fit,
    var_signal = "target", var_group = "time",
    plot_range = c(-5, 7), ndraws = 100
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("bayesian_SDT_distribution_plot works for equal-variance models", {
  b_fit <- read_fixture("brms_sdt_eqvar.rds")
  p <- bayesian_SDT_distribution_plot(
    b_fit,
    var_signal = "target", var_group = "time",
    plot_range = c(-5, 7), ndraws = 100
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("bayesian_SDT_distribution_plot facets for a 3-variable model", {
  b_fit <- read_fixture("brms_sdt_3var.rds")
  p <- bayesian_SDT_distribution_plot(
    b_fit,
    var_signal = "target", var_group = "time", var_facet = "condition",
    plot_range = c(-5, 7), ndraws = 100
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("bayesian functions validate their inputs", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  expect_error(
    bayesian_SDT_distribution_plot(lm(mpg ~ hp, mtcars), "a", "b"),
    "brmsfit"
  )
  expect_error(
    bayesian_SDT_distribution_plot(b_fit, "nope", "time"),
    "not a variable"
  )
  expect_error(
    bayesian_SDT_distribution_plot(b_fit, "target", "time", ci = 1.5),
    "between 0 and 1"
  )
  expect_error(bayesian_roc_plot(b_fit, "target", "value"), "exactly 2 levels")
})

test_that("non-cumulative models are rejected", {
  b_meta <- read_fixture("brms_meta.rds")
  expect_error(bayesian_roc_plot(b_meta, "a", "b"), "cumulative")
  expect_error(bayesian_SDT_distribution_plot(b_meta, "a", "b"), "cumulative")
})

test_that("ordinal_model_ppd_check returns a ggplot", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  p <- ordinal_model_ppd_check(b_fit, group_vars = c("target", "time"), ndraws = 10)
  expect_s3_class(p, "ggplot")
  expect_renders(p)
})

test_that("ordinal_model_ppd_check validates group_vars", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  expect_error(ordinal_model_ppd_check(b_fit, group_vars = 2), "character vector")
  expect_error(
    ordinal_model_ppd_check(b_fit, group_vars = c("target", "nope")),
    "not found in the model data"
  )
})

test_that("bayesian_forest returns a ggplot for a nested meta-analysis", {
  b_meta <- read_fixture("brms_meta.rds")
  p <- bayesian_forest(b_meta)
  expect_s3_class(p, "ggplot")
  expect_renders(p)
})

test_that("bayesian_forest rejects models without nested grouping terms", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  expect_error(bayesian_forest(b_fit), "nested grouping terms")
})
