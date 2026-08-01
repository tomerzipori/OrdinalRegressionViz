# Shared model fits and fixtures for the test suite.
#
# clm models are cheap and fitted once per test run (lazily, cached).
# brms models are expensive: they are pre-fitted by data-raw/make_test_fixtures.R
# and stored in tests/testthat/fixtures/; tests skip when a fixture is absent
# (e.g. on CI or in the built package).

.test_fits <- new.env(parent = emptyenv())

control_ratings <- function() {
  droplevels(subset(sdt_ratings, condition == "control"))
}

clm_2var <- function(link = "probit") {
  testthat::skip_if_not_installed("ordinal")
  key <- paste0("clm_2var_", link)
  if (is.null(.test_fits[[key]])) {
    .test_fits[[key]] <- ordinal::clm(
      value ~ target * time,
      data = control_ratings(), link = link
    )
  }
  .test_fits[[key]]
}

clm_2var_additive <- function() {
  testthat::skip_if_not_installed("ordinal")
  if (is.null(.test_fits$clm_2var_additive)) {
    .test_fits$clm_2var_additive <- ordinal::clm(
      value ~ target + time,
      data = control_ratings(), link = "probit"
    )
  }
  .test_fits$clm_2var_additive
}

clm_3var <- function() {
  testthat::skip_if_not_installed("ordinal")
  if (is.null(.test_fits$clm_3var)) {
    .test_fits$clm_3var <- ordinal::clm(
      value ~ target * time * condition,
      data = sdt_ratings, link = "probit"
    )
  }
  .test_fits$clm_3var
}

read_fixture <- function(name) {
  path <- testthat::test_path("fixtures", name)
  if (!file.exists(path)) {
    testthat::skip(sprintf("brms fixture %s not available", name))
  }
  testthat::skip_if_not_installed("brms")
  readRDS(path)
}

# Render a ggplot/patchwork object on a null device so that layer and scale
# errors surface.
expect_renders <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::expect_no_error(print(p))
}
