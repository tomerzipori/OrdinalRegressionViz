test_that("empirical_roc_points computes cumulative proportions per group", {
  df <- data.frame(
    value = factor(c(1, 1, 2, 3, 1, 2, 3, 3), levels = 1:3),
    target = factor(rep(c("fake", "true"), each = 4), levels = c("fake", "true")),
    time = factor("pre", levels = c("pre", "post"))
  )
  df <- df[rep(1:8, 2), ]
  df$time <- factor(rep(c("pre", "post"), each = 8), levels = c("pre", "post"))

  out <- empirical_roc_points(df, "value", "target", "time")
  expect_named(out, c("time", "Sensitivity", "FAR"))
  # 2 groups x (3 - 1) thresholds
  expect_equal(nrow(out), 4)
  pre <- out[out$time == "pre", ]
  # fake (signal): 2/4 at rating 1, 3/4 at <= 2
  expect_equal(pre$Sensitivity, c(0.5, 0.75))
  # true (noise): 1/4 at rating 1, 2/4 at <= 2
  expect_equal(pre$FAR, c(0.25, 0.5))
})

test_that("show_empirical adds a point layer to roc_plot", {
  p_plain <- roc_plot(clm_2var(), var_signal = "target", var_group = "time")
  p_emp <- roc_plot(
    clm_2var(),
    var_signal = "target", var_group = "time",
    show_empirical = TRUE
  )
  expect_equal(length(p_emp$layers), length(p_plain$layers) + 1)
  expect_renders(p_emp)
})

test_that("show_empirical works on faceted and Bayesian ROC plots", {
  p <- roc_plot(
    clm_3var(),
    var_signal = "target", var_group = "time", var_facet = "condition",
    show_empirical = TRUE
  )
  expect_renders(p)

  b_fit <- read_fixture("brms_sdt_2var.rds")
  p_b <- bayesian_roc_plot(
    b_fit,
    var_signal = "target", var_group = "time",
    show_empirical = TRUE
  )
  expect_renders(p_b)
})
