# Simulate the `sdt_ratings` example dataset.
#
# Trial-level truth ratings of true and fake news headlines, in a 2x2x2
# signal-detection design inspired by inoculation studies: participants in
# two conditions (control / inoculation) rate headlines (true / fake) before
# and after the intervention (pre / post) on a 6-point perceived-truth scale
# (1 = "definitely fake" ... 6 = "definitely true").
#
# Generative model: an unequal-variance probit SDT model on a latent
# "perceived truth" scale. Fake headlines are centered at 0 (SD 1); true
# headlines are shifted by d' (SD 1.2). Inoculation raises post-intervention
# discrimination and makes the response criteria stricter.

set.seed(20260731)

n_subj_per_cond <- 25
n_items <- 20 # per target type per time point

conditions <- c("control", "inoculation")
times <- c("pre", "post")
targets <- c("fake", "true")

thresholds <- c(-1.2, -0.4, 0.35, 1.1, 1.9)
sd_true <- 1.2 # unequal-variance SDT
sd_subj <- 0.35

d_prime <- function(condition, time) {
  if (time == "pre") {
    1.3
  } else if (condition == "control") {
    1.2
  } else {
    1.9
  }
}

criterion_shift <- function(condition, time) {
  # Inoculated participants become stricter (rate everything as less true).
  if (condition == "inoculation" && time == "post") 0.3 else 0
}

cells <- expand.grid(
  subj = seq_len(n_subj_per_cond),
  condition = conditions,
  time = times,
  target = targets,
  item = seq_len(n_items),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

subj_effect <- stats::rnorm(n_subj_per_cond * length(conditions), sd = sd_subj)
names(subj_effect) <- paste(
  rep(conditions, each = n_subj_per_cond),
  rep(seq_len(n_subj_per_cond), times = length(conditions))
)

latent <- mapply(
  function(subj, condition, time, target) {
    mu <- if (target == "true") d_prime(condition, time) else 0
    sdv <- if (target == "true") sd_true else 1
    mu + subj_effect[[paste(condition, subj)]] + stats::rnorm(1, sd = sdv)
  },
  cells$subj, cells$condition, cells$time, cells$target
)

shift <- mapply(criterion_shift, cells$condition, cells$time)
value <- vapply(
  seq_along(latent),
  function(i) findInterval(latent[i], thresholds + shift[i]) + 1,
  numeric(1)
)

sdt_ratings <- data.frame(
  id = factor(sprintf("S%02d", cells$subj + (cells$condition == "inoculation") * n_subj_per_cond)),
  condition = factor(cells$condition, levels = conditions),
  time = factor(cells$time, levels = times),
  target = factor(cells$target, levels = targets),
  value = factor(value, levels = 1:6, ordered = TRUE)
)
sdt_ratings <- sdt_ratings[order(sdt_ratings$id, sdt_ratings$time, sdt_ratings$target), ]
rownames(sdt_ratings) <- NULL

save(sdt_ratings, file = "data/sdt_ratings.rda", compress = "bzip2")
