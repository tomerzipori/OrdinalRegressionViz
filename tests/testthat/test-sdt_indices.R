test_that("sdt_indices.clm recovers the simulation parameters", {
  out <- sdt_indices(clm_2var(), var_signal = "target", var_group = "time")
  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "index", "estimate", "lower", "upper"))
  # 2 cells x (3 scalar indices + 5 criteria)
  expect_equal(nrow(out), 2 * 8)

  pre <- out[out$time == "pre", ]
  d_pre <- pre$estimate[pre$index == "d_prime"]
  # Generating d' = 1.3 with SD ~1.2 for signal and a subject intercept;
  # the fitted marginal d' lands near 1.2.
  expect_gt(d_pre, 0.8)
  expect_lt(d_pre, 1.6)
  expect_equal(pre$estimate[pre$index == "sigma_ratio"], 1)
  auc_pre <- pre$estimate[pre$index == "auc"]
  expect_gt(auc_pre, 0.7)
  expect_lt(auc_pre, 0.95)
  # probit closed form is consistent with its inputs
  expect_equal(auc_pre, unname(stats::pnorm(d_pre / sqrt(2))), tolerance = 1e-8)
  # criteria are ordered
  crit <- pre$estimate[startsWith(pre$index, "criterion")]
  expect_length(crit, 5)
  expect_true(all(diff(crit) > 0))
  # clm gives point estimates only
  expect_true(all(is.na(out$lower)))
})

test_that("sdt_indices.clm integrates the AUC for non-probit links", {
  out <- sdt_indices(clm_2var(link = "logit"), var_signal = "target", var_group = "time")
  auc <- out$estimate[out$index == "auc"]
  expect_true(all(auc > 0.5 & auc < 1))
})

test_that("sdt_indices works without a grouping variable", {
  skip_if_not_installed("ordinal")
  fit <- ordinal::clm(value ~ target, data = control_ratings(), link = "probit")
  out <- sdt_indices(fit, var_signal = "target")
  expect_named(out, c("index", "estimate", "lower", "upper"))
  expect_equal(nrow(out), 8)
})

test_that("sdt_indices handles a facet variable", {
  out <- sdt_indices(
    clm_3var(),
    var_signal = "target", var_group = "time", var_facet = "condition"
  )
  expect_named(out, c("time", "condition", "index", "estimate", "lower", "upper"))
  expect_equal(nrow(out), 4 * 8)
  # Post-inoculation discrimination is higher than pre (d' 1.9 vs 1.3 in
  # the simulation).
  d <- function(t, cond) {
    out$estimate[out$index == "d_prime" & out$time == t & out$condition == cond]
  }
  expect_gt(d("post", "inoculation"), d("pre", "inoculation"))
})

test_that("sdt_indices.brmsfit returns posterior summaries", {
  b_fit <- read_fixture("brms_sdt_2var.rds")
  out <- sdt_indices(b_fit, var_signal = "target", var_group = "time", ci = 0.9)
  expect_named(out, c("time", "index", "estimate", "lower", "upper"))
  expect_equal(nrow(out), 2 * 8)
  expect_true(all(out$lower < out$estimate & out$estimate < out$upper))

  # The unequal-variance fixture should detect sigma_ratio away from 1
  # (generating value ~1.2) and agree with the clm point estimates.
  pre <- out[out$time == "pre", ]
  expect_gt(pre$estimate[pre$index == "sigma_ratio"], 0.9)
  expect_lt(pre$estimate[pre$index == "sigma_ratio"], 1.6)

  clm_out <- sdt_indices(clm_2var(), var_signal = "target", var_group = "time")
  d_brms <- pre$estimate[pre$index == "d_prime"]
  d_clm <- clm_out$estimate[clm_out$index == "d_prime" & clm_out$time == "pre"]
  expect_equal(d_brms, d_clm, tolerance = 0.25)
})

test_that("sdt_indices.brmsfit supports equal-variance models and ndraws", {
  b_fit <- read_fixture("brms_sdt_eqvar.rds")
  out <- sdt_indices(b_fit, var_signal = "target", var_group = "time", ndraws = 100)
  expect_true(all(out$estimate[out$index == "sigma_ratio"] == 1))
})

test_that("sdt_indices validates its inputs", {
  expect_error(sdt_indices(lm(mpg ~ hp, mtcars), "a"), "clm.*brmsfit")
  fit <- clm_2var()
  expect_error(sdt_indices(fit, "nope"), "not a variable")
  b_meta <- read_fixture("brms_meta.rds")
  expect_error(sdt_indices(b_meta, "Author"), "cumulative")
})
