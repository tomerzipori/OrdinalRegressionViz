test_that("check_brmsfit rejects non-brmsfit objects", {
  expect_error(check_brmsfit(lm(mpg ~ hp, mtcars)), "brmsfit")
  expect_error(check_brmsfit(NULL), "brmsfit")
})

test_that("check_clm rejects non-clm objects", {
  expect_error(check_clm(lm(mpg ~ hp, mtcars)), "clm")
  expect_error(check_clm(list()), "clm")
})

test_that("check_clm requires a stored model frame", {
  skip_if_not_installed("ordinal")
  fit <- ordinal::clm(value ~ target, data = control_ratings(), model = FALSE)
  expect_error(check_clm(fit), "model frame")
})

test_that("check_factor_var validates name, type, and levels", {
  df <- data.frame(
    f2 = factor(rep(c("a", "b"), 3)),
    f3 = factor(rep(c("a", "b", "c"), 2)),
    num = 1:6
  )
  expect_silent(check_factor_var(df, "f2", "var_signal"))
  expect_error(check_factor_var(df, c("f2", "f3"), "var_signal"), "single variable name")
  expect_error(check_factor_var(df, "nope", "var_signal"), "not a variable")
  expect_error(check_factor_var(df, "num", "var_signal"), "must be a factor")
  expect_error(check_factor_var(df, "f3", "var_signal"), "exactly 2 levels")
})

test_that("check_prob validates probabilities", {
  expect_silent(check_prob(0.9, "ci"))
  expect_error(check_prob(0, "ci"), "between 0 and 1")
  expect_error(check_prob(1, "ci"), "between 0 and 1")
  expect_error(check_prob(c(0.1, 0.2), "ci"), "between 0 and 1")
  expect_error(check_prob("a", "ci"), "between 0 and 1")
})

test_that("check_labels2 defaults to title-cased levels", {
  expect_equal(check_labels2(NULL, c("pre", "post"), "x"), c("Pre", "Post"))
  expect_equal(check_labels2(c("A", "B"), c("pre", "post"), "x"), c("A", "B"))
  expect_error(check_labels2("only-one", c("pre", "post"), "x"), "length 2")
  expect_error(check_labels2(1:2, c("pre", "post"), "x"), "length 2")
})
