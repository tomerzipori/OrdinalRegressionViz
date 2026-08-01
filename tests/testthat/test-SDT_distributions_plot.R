test_that("SDT_distributions_plot returns a patchwork for a 2-variable model", {
  p <- SDT_distributions_plot(clm_2var(), var_signal = "target", var_group = "time")
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("SDT_distributions_plot facets for a 3-variable model", {
  p <- SDT_distributions_plot(
    clm_3var(),
    var_signal = "target", var_group = "time", var_facet = "condition"
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("SDT_distributions_plot handles additive models (no interaction)", {
  p <- SDT_distributions_plot(
    clm_2var_additive(),
    var_signal = "target", var_group = "time"
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("SDT_distributions_plot supports non-probit links", {
  p <- SDT_distributions_plot(
    clm_2var(link = "logit"),
    var_signal = "target", var_group = "time"
  )
  expect_s3_class(p, "patchwork")
})

test_that("SDT_distributions_plot draws the fitted thresholds", {
  fit <- clm_2var()
  p <- SDT_distributions_plot(fit, var_signal = "target", var_group = "time")
  # The first panel's threshold lines are the model's alpha estimates.
  xints <- ggplot2::layer_data(p[[1]], 3)$xintercept
  expect_equal(sort(xints), sort(unname(fit$alpha)), tolerance = 1e-8)
})

test_that("SDT_distributions_plot validates its inputs", {
  fit <- clm_2var()
  expect_error(SDT_distributions_plot(lm(mpg ~ hp, mtcars), "a", "b"), "clm")
  expect_error(SDT_distributions_plot(fit, "nope", "time"), "not a variable")
  expect_error(
    SDT_distributions_plot(fit, "target", "time", signal_labels = "x"),
    "length 2"
  )
})
