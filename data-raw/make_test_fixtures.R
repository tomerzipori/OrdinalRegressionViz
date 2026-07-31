# Fit the small brms models used by the test suite and (optionally) the
# vignette figures. Run manually; the resulting .rds files live in
# tests/testthat/fixtures/ and tests skip when they are absent.
#
# The models are deliberately small (few chains/iterations) - they exist to
# exercise the plotting code, not to make inferences.

library(brms)

load("data/sdt_ratings.rda")

fixture_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

opts <- list(
  chains = 2, iter = 600, warmup = 400, cores = 2,
  seed = 20260731, refresh = 0,
  backend = "rstan"
)

# 1. Unequal-variance SDT model, 2 predictors (control condition only)
b_sdt_2var <- do.call(brm, c(list(
  formula = bf(value ~ target * time, disc ~ target * time),
  family = cumulative("probit"),
  data = droplevels(subset(sdt_ratings, condition == "control"))
), opts))
saveRDS(b_sdt_2var, file.path(fixture_dir, "brms_sdt_2var.rds"), compress = "xz")

# 1b. Equal-variance variant (no disc part)
b_sdt_eqvar <- do.call(brm, c(list(
  formula = value ~ target * time,
  family = cumulative("probit"),
  data = droplevels(subset(sdt_ratings, condition == "control"))
), opts))
saveRDS(b_sdt_eqvar, file.path(fixture_dir, "brms_sdt_eqvar.rds"), compress = "xz")

# 2. Same model with a third (faceting) predictor
b_sdt_3var <- do.call(brm, c(list(
  formula = bf(
    value ~ target * time * condition,
    disc ~ target * time * condition
  ),
  family = cumulative("probit"),
  data = sdt_ratings
), opts))
saveRDS(b_sdt_3var, file.path(fixture_dir, "brms_sdt_3var.rds"), compress = "xz")

# 3. Random-effects meta-analysis with nested Author/Study grouping
set.seed(20260731)
n_authors <- 4
studies_per_author <- 3
meta_data <- expand.grid(
  Author = paste0("Author", seq_len(n_authors)),
  StudyNum = seq_len(studies_per_author)
)
meta_data$Study <- paste0("Study", meta_data$StudyNum)
author_effect <- rnorm(n_authors, sd = 0.2)
meta_data$sei <- runif(nrow(meta_data), 0.08, 0.25)
meta_data$yi <- 0.4 +
  author_effect[as.integer(factor(meta_data$Author))] +
  rnorm(nrow(meta_data), sd = 0.15) +
  rnorm(nrow(meta_data), sd = meta_data$sei)

b_meta <- do.call(brm, c(list(
  formula = yi | se(sei) ~ 1 + (1 | Author / Study),
  data = meta_data,
  control = list(adapt_delta = 0.95)
), opts))
saveRDS(b_meta, file.path(fixture_dir, "brms_meta.rds"), compress = "xz")

cat("Fixture sizes (bytes):\n")
print(file.size(list.files(fixture_dir, full.names = TRUE)))
