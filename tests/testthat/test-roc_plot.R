test_that("roc_plot returns a ggplot for a 2-variable model", {
  p <- roc_plot(clm_2var(), var_signal = "target", var_group = "time")
  expect_s3_class(p, "ggplot")
  expect_renders(p)
})

test_that("roc_plot facets into a patchwork for a 3-variable model", {
  p <- roc_plot(
    clm_3var(),
    var_signal = "target", var_group = "time", var_facet = "condition"
  )
  expect_s3_class(p, "patchwork")
  expect_renders(p)
})

test_that("roc_plot works for non-probit links", {
  p <- roc_plot(clm_2var(link = "logit"), var_signal = "target", var_group = "time")
  expect_s3_class(p, "ggplot")
})

test_that("roc_plot accepts custom labels and CI", {
  p <- roc_plot(
    clm_2var(),
    var_signal = "target", var_group = "time",
    CI = 0.8, group_labels = c("Before", "After"), ttl = "ROC"
  )
  expect_s3_class(p, "ggplot")
  expect_renders(p)
})

test_that("roc_plot validates its inputs", {
  fit <- clm_2var()
  expect_error(roc_plot(lm(mpg ~ hp, mtcars), "a", "b"), "clm")
  expect_error(roc_plot(fit, "nope", "time"), "not a variable")
  expect_error(roc_plot(fit, "target", "value"), "exactly 2 levels")
  expect_error(roc_plot(fit, "target", "time", CI = 2), "between 0 and 1")
  expect_error(
    roc_plot(fit, "target", "time", group_labels = "one"),
    "length 2"
  )
})
