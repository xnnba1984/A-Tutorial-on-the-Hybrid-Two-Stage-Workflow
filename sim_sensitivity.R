## Focused sensitivity analysis for sample size, candidate covariate dimension,
## and Stage 2 learner complexity in the SBR major revision.

suppressPackageStartupMessages({
  library(dplyr)
  library(grf)
  library(sandwich)
})

ROOT <- normalizePath(getwd())
OUT_DIR <- Sys.getenv(
  "SBR_SENS_OUT_DIR",
  file.path(ROOT, "result", "information_sensitivity")
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("SBR_SENS_B", "300"))
NUM_TREES <- as.integer(Sys.getenv("SBR_SENS_TREES", "1000"))
NUM_CORES <- as.integer(Sys.getenv(
  "SBR_SENS_CORES",
  as.character(max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L)))
))

N_VALUES <- c(500L, 1000L, 2000L)
P_VALUES <- c(3L, 20L)
SCENARIOS <- c(
  "Constant benefit without HTE",
  "Strong quantitative HTE",
  "Strong qualitative HTE"
)

expit <- function(z) 1 / (1 + exp(-z))

scenario_tau <- function(x1, scenario) {
  switch(
    scenario,
    "Constant benefit without HTE" = rep(0.04, length(x1)),
    "Strong quantitative HTE" = 0.06 + 0.050 * tanh(x1),
    "Strong qualitative HTE" = 0.180 * tanh(x1),
    stop("Unknown scenario: ", scenario)
  )
}

simulate_trial <- function(n, p, scenario, seed) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x[, 3] <- rbinom(n, 1, 0.5)
  colnames(x) <- paste0("X", seq_len(p))

  p0 <- 0.25 + 0.35 * expit(
    -0.4 + 0.45 * x[, 1] - 0.25 * x[, 2] + 0.35 * x[, 3]
  )
  tau <- scenario_tau(x[, 1], scenario)
  p1 <- p0 + tau
  stopifnot(all(p1 >= 0), all(p1 <= 1))

  a <- rbinom(n, 1, 0.5)
  y0 <- rbinom(n, 1, p0)
  y1 <- rbinom(n, 1, p1)
  y <- ifelse(a == 1, y1, y0)

  list(
    X = as.data.frame(x),
    A = a,
    Y = y,
    p0 = p0,
    p1 = p1,
    tau = tau
  )
}

stratified_split <- function(a, seed) {
  set.seed(seed)
  split_id <- rep(NA_character_, length(a))
  for (arm in sort(unique(a))) {
    idx <- sample(which(a == arm))
    n_arm <- length(idx)
    n_train <- floor(0.50 * n_arm)
    n_tune <- floor(0.25 * n_arm)
    split_id[idx[seq_len(n_train)]] <- "train"
    split_id[idx[(n_train + 1):(n_train + n_tune)]] <- "tune"
    split_id[idx[(n_train + n_tune + 1):n_arm]] <- "test"
  }
  list(
    train = which(split_id == "train"),
    tune = which(split_id == "tune"),
    test = which(split_id == "test")
  )
}

stage1_test <- function(dat) {
  x_names <- names(dat$X)
  df <- data.frame(Y = dat$Y, A = dat$A, dat$X)
  form <- as.formula(paste("Y ~ A * (", paste(x_names, collapse = " + "), ")"))
  fit <- lm(form, data = df)
  interaction_terms <- grep("^A:", names(coef(fit)), value = TRUE)
  beta <- coef(fit)[interaction_terms]
  vc <- sandwich::vcovHC(fit, type = "HC3")[
    interaction_terms,
    interaction_terms,
    drop = FALSE
  ]
  stat <- tryCatch(
    as.numeric(crossprod(beta, solve(vc, beta))),
    error = function(e) NA_real_
  )
  p_value <- if (is.finite(stat)) {
    pchisq(stat, df = length(beta), lower.tail = FALSE)
  } else {
    NA_real_
  }
  c(p_value = p_value, reject = as.integer(p_value < 0.05))
}

fit_outcome_nuisance <- function(dat, train_idx) {
  x_names <- names(dat$X)
  dtr <- data.frame(Y = dat$Y[train_idx], A = dat$A[train_idx], dat$X[train_idx, ])
  form <- as.formula(paste("Y ~", paste(x_names, collapse = " + ")))
  list(
    m0 = glm(form, data = dtr[dtr$A == 0, ], family = quasibinomial()),
    m1 = glm(form, data = dtr[dtr$A == 1, ], family = quasibinomial())
  )
}

predict_outcome_nuisance <- function(fits, x) {
  list(
    mu0 = pmin(pmax(plogis(predict(fits$m0, newdata = x)), 1e-4), 1 - 1e-4),
    mu1 = pmin(pmax(plogis(predict(fits$m1, newdata = x)), 1e-4), 1 - 1e-4)
  )
}

dr_pseudo <- function(a, y, mu1, mu0) {
  mu1 - mu0 + 2 * a * (y - mu1) - 2 * (1 - a) * (y - mu0)
}

dr_value <- function(policy, a, y, mu1, mu0) {
  m_policy <- ifelse(policy == 1, mu1, mu0)
  m_observed <- ifelse(a == 1, mu1, mu0)
  mean(m_policy + ifelse(a == policy, 2 * (y - m_observed), 0))
}

centered_auqc <- function(score, dr_score) {
  ord <- order(score, decreasing = TRUE)
  n <- length(ord)
  frac <- seq_len(n) / n
  centered <- cumsum(dr_score[ord]) - frac * sum(dr_score)
  frac0 <- c(0, frac)
  centered0 <- c(0, centered)
  area <- sum(
    diff(frac0) *
      (centered0[-1] + centered0[-length(centered0)]) / 2
  )
  area / n
}

threshold_grid <- function(score) {
  unique(c(
    min(score) - 1e-6,
    as.numeric(quantile(score, probs = seq(0.01, 0.99, length.out = 99), names = FALSE)),
    max(score) + 1e-6
  ))
}

fit_scores <- function(dat, split, seed) {
  x_train <- as.matrix(dat$X[split$train, ])
  x_tune <- as.matrix(dat$X[split$tune, ])
  x_test <- as.matrix(dat$X[split$test, ])
  a_train <- dat$A[split$train]
  y_train <- dat$Y[split$train]

  cf <- causal_forest(
    x_train,
    y_train,
    a_train,
    num.trees = NUM_TREES,
    num.threads = 1L,
    seed = seed
  )
  cf_tune <- predict(cf, x_tune)$predictions
  cf_test <- predict(cf, x_test)$predictions

  x_names <- names(dat$X)
  dtr <- data.frame(Y = y_train, A = a_train, dat$X[split$train, ])
  linear_form <- as.formula(paste(
    "Y ~ A * (",
    paste(x_names, collapse = " + "),
    ")"
  ))
  linear_fit <- lm(linear_form, data = dtr)

  linear_predict <- function(x) {
    d0 <- data.frame(A = 0, x)
    d1 <- data.frame(A = 1, x)
    as.numeric(predict(linear_fit, newdata = d1) - predict(linear_fit, newdata = d0))
  }

  list(
    causal_forest = list(tune = cf_tune, test = cf_test),
    linear_interaction = list(
      tune = linear_predict(dat$X[split$tune, ]),
      test = linear_predict(dat$X[split$test, ])
    )
  )
}

evaluate_learner <- function(dat, split, scores, nuisance) {
  tune_nuis <- predict_outcome_nuisance(nuisance, dat$X[split$tune, ])
  test_nuis <- predict_outcome_nuisance(nuisance, dat$X[split$test, ])

  a_tune <- dat$A[split$tune]
  y_tune <- dat$Y[split$tune]
  a_test <- dat$A[split$test]
  y_test <- dat$Y[split$test]
  dr_test <- dr_pseudo(a_test, y_test, test_nuis$mu1, test_nuis$mu0)

  t_grid <- threshold_grid(scores$tune)
  tune_values <- vapply(t_grid, function(t) {
    policy <- as.integer(scores$tune >= t)
    dr_value(policy, a_tune, y_tune, tune_nuis$mu1, tune_nuis$mu0)
  }, numeric(1))
  threshold <- t_grid[which.max(tune_values)]

  tune_all <- dr_value(
    rep(1L, length(a_tune)),
    a_tune,
    y_tune,
    tune_nuis$mu1,
    tune_nuis$mu0
  )
  tune_none <- dr_value(
    rep(0L, length(a_tune)),
    a_tune,
    y_tune,
    tune_nuis$mu1,
    tune_nuis$mu0
  )
  fixed_treat <- as.integer(tune_all >= tune_none)

  test_policy <- as.integer(scores$test >= threshold)
  estimated_policy_value <- dr_value(
    test_policy,
    a_test,
    y_test,
    test_nuis$mu1,
    test_nuis$mu0
  )
  estimated_fixed_value <- dr_value(
    rep(fixed_treat, length(a_test)),
    a_test,
    y_test,
    test_nuis$mu1,
    test_nuis$mu0
  )

  p0_test <- dat$p0[split$test]
  p1_test <- dat$p1[split$test]
  true_policy_value <- mean(ifelse(test_policy == 1, p1_test, p0_test))
  true_best_fixed_value <- max(mean(p0_test), mean(p1_test))
  true_selected_fixed_value <- if (fixed_treat == 1) mean(p1_test) else mean(p0_test)
  true_oracle_policy_value <- mean(pmax(p0_test, p1_test))

  c(
    centered_auqc = centered_auqc(scores$test, dr_test),
    estimated_incremental_value = estimated_policy_value - estimated_fixed_value,
    true_incremental_value_best_fixed = true_policy_value - true_best_fixed_value,
    true_incremental_value_selected_fixed = true_policy_value - true_selected_fixed_value,
    oracle_personalization_gain = true_oracle_policy_value - true_best_fixed_value,
    learned_policy_regret = true_oracle_policy_value - true_policy_value,
    treated_fraction = mean(test_policy),
    selected_fixed_treat_all = fixed_treat,
    threshold = threshold
  )
}

setting_grid <- expand.grid(
  scenario = SCENARIOS,
  n = N_VALUES,
  p = P_VALUES,
  stringsAsFactors = FALSE
)
setting_grid$setting_id <- seq_len(nrow(setting_grid))

tasks <- do.call(rbind, lapply(seq_len(nrow(setting_grid)), function(i) {
  cbind(setting_grid[rep(i, B), ], replicate = seq_len(B))
}))

run_task <- function(task_index) {
  task <- tasks[task_index, ]
  seed <- 202608230L + 100000L * task$setting_id + task$replicate
  dat <- simulate_trial(task$n, task$p, task$scenario, seed)
  split <- stratified_split(dat$A, seed + 1L)
  st1 <- stage1_test(dat)
  nuisance <- fit_outcome_nuisance(dat, split$train)
  scores <- fit_scores(dat, split, seed + 2L)

  out <- lapply(names(scores), function(learner) {
    metrics <- evaluate_learner(dat, split, scores[[learner]], nuisance)
    data.frame(
      setting_id = task$setting_id,
      scenario = task$scenario,
      n = task$n,
      p = task$p,
      replicate = task$replicate,
      learner = learner,
      stage1_p_value = unname(st1["p_value"]),
      stage1_reject = unname(st1["reject"]),
      n_train = length(split$train),
      n_tune = length(split$tune),
      n_test = length(split$test),
      t(metrics),
      check.names = FALSE
    )
  })
  bind_rows(out)
}

cat(sprintf(
  "[Sensitivity] settings=%d, B=%d, tasks=%d, trees=%d, cores=%d\n",
  nrow(setting_grid), B, nrow(tasks), NUM_TREES, NUM_CORES
))

if (.Platform$OS.type == "unix" && NUM_CORES > 1L) {
  result_list <- parallel::mclapply(
    seq_len(nrow(tasks)),
    run_task,
    mc.cores = NUM_CORES,
    mc.preschedule = FALSE
  )
} else {
  result_list <- lapply(seq_len(nrow(tasks)), run_task)
}

raw <- bind_rows(result_list)

mc_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

summary <- raw |>
  group_by(scenario, n, p, learner) |>
  summarise(
    B = n(),
    valid_B = sum(
      is.finite(stage1_reject) &
        is.finite(centered_auqc) &
        is.finite(estimated_incremental_value) &
        is.finite(true_incremental_value_best_fixed)
    ),
    stage1_rejection_rate = mean(stage1_reject, na.rm = TRUE),
    se_stage1_rejection_rate = mc_se(stage1_reject),
    mean_centered_auqc = mean(centered_auqc, na.rm = TRUE),
    se_centered_auqc = mc_se(centered_auqc),
    mean_estimated_incremental_value = mean(estimated_incremental_value, na.rm = TRUE),
    se_estimated_incremental_value = mc_se(estimated_incremental_value),
    mean_true_incremental_value_best_fixed = mean(true_incremental_value_best_fixed, na.rm = TRUE),
    se_true_incremental_value_best_fixed = mc_se(true_incremental_value_best_fixed),
    mean_true_incremental_value_selected_fixed = mean(true_incremental_value_selected_fixed, na.rm = TRUE),
    se_true_incremental_value_selected_fixed = mc_se(true_incremental_value_selected_fixed),
    mean_oracle_personalization_gain = mean(oracle_personalization_gain, na.rm = TRUE),
    se_oracle_personalization_gain = mc_se(oracle_personalization_gain),
    mean_learned_policy_regret = mean(learned_policy_regret, na.rm = TRUE),
    se_learned_policy_regret = mc_se(learned_policy_regret),
    positive_estimated_gain_rate = mean(estimated_incremental_value > 0, na.rm = TRUE),
    mean_treated_fraction = mean(treated_fraction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(scenario, n, p, learner)

write.csv(raw, file.path(OUT_DIR, "sim_information_sensitivity_raw.csv"), row.names = FALSE)
write.csv(summary, file.path(OUT_DIR, "sim_information_sensitivity_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sim_information_sensitivity_session_info.txt"))

cat("\n[Sensitivity summary]\n")
print(summary, n = nrow(summary), width = Inf)
cat("\n[Saved]\n")
cat(file.path(OUT_DIR, "sim_information_sensitivity_raw.csv"), "\n")
cat(file.path(OUT_DIR, "sim_information_sensitivity_summary.csv"), "\n")
cat(file.path(OUT_DIR, "sim_information_sensitivity_session_info.txt"), "\n")
