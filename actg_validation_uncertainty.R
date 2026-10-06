## Repeated split validation for the ACTG 175 IPCW case study.
## This script addresses optimism in the original cross-fitted ACTG analysis by
## selecting the CATE threshold on a tuning split and estimating final policy
## performance on an independent test split.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(grf)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(car)
})

set.seed(20260517)

ROOT <- getwd()
if (!file.exists(file.path(ROOT, "data/ACTG175.csv"))) {
  stop("Run this script from the repository root containing data/ACTG175.csv.")
}

OUT_DIR <- Sys.getenv(
  "ACTG_VALIDATION_OUT_DIR",
  file.path(ROOT, "result", "actg_validation")
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

HORIZON_DAYS <- 96 * 7
ALPHA_HARM <- 0.10
BENEFIT_MARGIN <- 0.0
B <- as.integer(Sys.getenv("ACTG_VALIDATION_B", "200"))
B_STAGE1_BOOT <- as.integer(Sys.getenv("ACTG_STAGE1_BOOT_B", "1000"))
NUM_TREES <- as.integer(Sys.getenv("ACTG_VALIDATION_TREES", "1500"))
NUM_THREADS <- as.integer(Sys.getenv(
  "ACTG_VALIDATION_THREADS",
  as.character(max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L)))
))

rename_if_present <- function(dat, from, to) {
  for (i in seq_along(from)) {
    if (from[i] %in% names(dat) && !(to[i] %in% names(dat))) {
      names(dat)[names(dat) == from[i]] <- to[i]
    }
  }
  dat
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

extract_chisq_test <- function(lh) {
  out <- as.data.frame(lh)
  data.frame(
    chisq = out$Chisq[2],
    df = out$Df[2],
    p_value = out$`Pr(>Chisq)`[2]
  )
}

prepare_actg <- function() {
  df_raw <- fread("data/ACTG175.csv", na.strings = c("", "NA", ".", "NaN")) |>
    as.data.frame()

  df_raw <- rename_if_present(
    df_raw,
    from = c("time", "cid", "trt", "cd4.baseline", "cd8.baseline"),
    to = c("days", "cens", "arms", "cd40", "cd80")
  )

  if (!("treat" %in% names(df_raw)) && ("arms" %in% names(df_raw))) {
    df_raw <- df_raw |>
      mutate(treat = ifelse(arms == 0, 0L, 1L))
  }

  baseline_candidates <- c(
    "age", "wtkg", "hemo", "homo", "drugs", "karnof",
    "oprior", "z30", "zprior", "preanti", "race", "gender",
    "symptom", "cd40", "cd80", "str2", "strat"
  )

  df <- df_raw |>
    mutate(
      A = as.integer(treat),
      event_by_h = as.integer(!is.na(days) & !is.na(cens) &
                                days <= HORIZON_DAYS & cens == 1L),
      early_censored = as.integer(!is.na(days) & !is.na(cens) &
                                    days < HORIZON_DAYS & cens == 0L),
      R_h = as.integer(event_by_h == 1L | days >= HORIZON_DAYS),
      Y = ifelse(event_by_h == 1L, 0,
                 ifelse(days >= HORIZON_DAYS, 1, NA_real_)),
      censor_time = pmin(days, HORIZON_DAYS),
      censor_event = early_censored
    ) |>
    filter(arms != 3)

  X_names <- setdiff(intersect(baseline_candidates, names(df)), "zprior")
  keep <- complete.cases(df[, c(
    "A", "days", "cens", "censor_time", "censor_event", "R_h", X_names
  )])

  dat_all <- df[keep, c(
    "A", "Y", "days", "cens", "arms", "event_by_h", "early_censored",
    "R_h", "censor_time", "censor_event", X_names
  )]

  factor_names <- intersect(
    c("hemo", "homo", "drugs", "race", "gender", "oprior",
      "z30", "symptom", "str2", "strat"),
    X_names
  )
  for (nm in factor_names) dat_all[[nm]] <- factor(dat_all[[nm]])

  list(df_raw = df_raw, dat_all = dat_all, X_names = X_names)
}

prep <- prepare_actg()
dat_all <- prep$dat_all
X_names <- prep$X_names
rhs <- paste(X_names, collapse = " + ")
censor_names <- intersect(
  c("A", "age", "wtkg", "karnof", "preanti", "cd40", "cd80", "strat"),
  names(dat_all)
)
censor_rhs <- paste(censor_names, collapse = " + ")
X_num <- model.matrix(as.formula(paste("~", rhs, "- 1")), data = dat_all)

fit_censor_model <- function(train) {
  f_cens <- as.formula(paste("Surv(censor_time, censor_event) ~", censor_rhs))
  fit <- survival::coxph(f_cens, data = train, ties = "efron", x = TRUE)
  bh <- survival::basehaz(fit, centered = FALSE)
  list(fit = fit, basehaz = bh)
}

predict_censor_survival <- function(censor_model, eval) {
  bh <- censor_model$basehaz
  eval_time <- pmin(eval$days, HORIZON_DAYS)
  H0 <- approx(
    bh$time, bh$hazard, xout = eval_time,
    method = "constant", f = 0, rule = 2, yleft = 0
  )$y
  lp <- predict(
    censor_model$fit,
    newdata = eval,
    type = "lp",
    reference = "zero"
  )
  pmax(exp(-H0 * exp(lp)), 0.02)
}

fit_ipcw <- function(train, eval) {
  predict_censor_survival(fit_censor_model(train), eval)
}

fit_stage1_model <- function(data) {
  dat_stage1 <- data |>
    mutate(
      G_hat = fit_ipcw(data, data),
      ipcw = ifelse(R_h == 1L, 1 / G_hat, 0)
    ) |>
    filter(R_h == 1L)

  f_int <- as.formula(paste("Y ~ A * (", rhs, ")"))
  lm_int <- lm(f_int, data = dat_stage1, weights = ipcw)
  interaction_terms <- grep("^A:", names(coef(lm_int)), value = TRUE)
  interaction_terms <- interaction_terms[!is.na(coef(lm_int)[interaction_terms])]
  list(
    fit = lm_int,
    interaction_terms = interaction_terms,
    beta = coef(lm_int)[interaction_terms]
  )
}

stage1_intervals <- function() {
  point <- fit_stage1_model(dat_all)
  interaction_terms <- point$interaction_terms
  beta <- point$beta

  set.seed(202608230)
  bootstrap_indices <- replicate(
    B_STAGE1_BOOT,
    sample.int(nrow(dat_all), nrow(dat_all), replace = TRUE),
    simplify = FALSE
  )
  bootstrap_one <- function(index) {
    tryCatch({
      fitted <- fit_stage1_model(dat_all[index, , drop = FALSE])
      out <- rep(NA_real_, length(interaction_terms))
      names(out) <- interaction_terms
      common <- intersect(interaction_terms, names(fitted$beta))
      out[common] <- fitted$beta[common]
      out
    }, error = function(e) {
      out <- rep(NA_real_, length(interaction_terms))
      names(out) <- interaction_terms
      out
    })
  }
  if (.Platform$OS.type == "unix" && NUM_THREADS > 1L) {
    bootstrap_list <- parallel::mclapply(
      bootstrap_indices,
      bootstrap_one,
      mc.cores = NUM_THREADS,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    bootstrap_list <- lapply(bootstrap_indices, bootstrap_one)
  }
  bootstrap_beta <- do.call(rbind, bootstrap_list)
  colnames(bootstrap_beta) <- interaction_terms
  valid <- complete.cases(bootstrap_beta)
  bootstrap_valid <- sum(valid)
  if (bootstrap_valid < max(100L, 0.8 * B_STAGE1_BOOT)) {
    stop("Too few valid Stage 1 bootstrap replicates: ", bootstrap_valid)
  }
  bootstrap_cov <- stats::cov(bootstrap_beta[valid, , drop = FALSE])
  global_chisq <- as.numeric(crossprod(beta, solve(bootstrap_cov, beta)))
  global_df <- length(beta)
  global_summary <- data.frame(
    chisq = global_chisq,
    df = global_df,
    p_value = pchisq(global_chisq, global_df, lower.tail = FALSE),
    bootstrap_B = B_STAGE1_BOOT,
    bootstrap_valid = bootstrap_valid,
    uncertainty_method = "patient-level bootstrap refitting censoring and outcome models"
  )

  prespec <- intersect(c("cd40", "karnof"), X_names)
  out <- lapply(prespec, function(x) {
    term <- paste0("A:", x)
    if (!(term %in% names(beta))) {
      return(data.frame(
        term = x, scale_label = NA_character_, coefficient = NA_real_,
        se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        scaled_coefficient = NA_real_, scaled_ci_low = NA_real_,
        scaled_ci_high = NA_real_, p_raw = NA_real_
      ))
    }
    mult <- ifelse(x == "cd40", 100, ifelse(x == "karnof", 10, 1))
    label <- ifelse(x == "cd40", "per 100 cells/mm^3", ifelse(x == "karnof", "per 10 points", "per unit"))
    est <- beta[term]
    boot_values <- bootstrap_beta[valid, term]
    se <- stats::sd(boot_values)
    ci <- stats::quantile(boot_values, c(0.025, 0.975), names = FALSE)
    p_raw <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
    data.frame(
      term = x,
      scale_label = label,
      coefficient = est,
      se = se,
      ci_low = ci[1],
      ci_high = ci[2],
      scaled_coefficient = est * mult,
      scaled_ci_low = ci[1] * mult,
      scaled_ci_high = ci[2] * mult,
      p_raw = p_raw,
      uncertainty_method = "patient-level bootstrap"
    )
  }) |>
    bind_rows() |>
    mutate(p_holm = p.adjust(p_raw, method = "holm"))

  list(
    global = global_summary,
    option_c = out,
    bootstrap_coefficients = bootstrap_beta,
    bootstrap_covariance = bootstrap_cov
  )
}

fit_split_models <- function(train_idx, seed) {
  dtr <- dat_all[train_idx, , drop = FALSE]

  censor_model <- fit_censor_model(dtr)
  G_tr <- predict_censor_survival(censor_model, dtr)
  ipcw_tr <- ifelse(dtr$R_h == 1L, 1 / G_tr, 0)

  e_fit <- glm(
    as.formula(paste("A ~", rhs)),
    data = dtr,
    family = binomial()
  )

  obs_tr <- dtr$R_h == 1L
  f_y <- as.formula(paste("Y ~", rhs))
  m0 <- glm(
    f_y,
    data = dtr[obs_tr & dtr$A == 0L, ],
    family = quasibinomial(),
    weights = ipcw_tr[obs_tr & dtr$A == 0L]
  )
  m1 <- glm(
    f_y,
    data = dtr[obs_tr & dtr$A == 1L, ],
    family = quasibinomial(),
    weights = ipcw_tr[obs_tr & dtr$A == 1L]
  )

  cf_obs_idx <- train_idx[obs_tr]
  cf <- causal_forest(
    X_num[cf_obs_idx, , drop = FALSE],
    dat_all$Y[cf_obs_idx],
    dat_all$A[cf_obs_idx],
    sample.weights = ipcw_tr[obs_tr],
    num.trees = NUM_TREES,
    num.threads = NUM_THREADS,
    seed = seed
  )

  list(
    censor_model = censor_model,
    e_fit = e_fit,
    m0 = m0,
    m1 = m1,
    causal_forest = cf
  )
}

predict_split <- function(model_bundle, eval_idx) {
  dev <- dat_all[eval_idx, , drop = FALSE]

  G_ev <- predict_censor_survival(model_bundle$censor_model, dev)
  e_hat <- plogis(predict(model_bundle$e_fit, newdata = dev))
  m0_hat <- plogis(predict(model_bundle$m0, newdata = dev))
  m1_hat <- plogis(predict(model_bundle$m1, newdata = dev))
  tau_hat <- predict(
    model_bundle$causal_forest,
    X_num[eval_idx, , drop = FALSE]
  )$predictions

  R_ev <- dev$R_h
  Y_ev <- ifelse(R_ev == 1L, dev$Y, 0)
  A_ev <- dev$A
  e_ev <- pmin(pmax(e_hat, 0.02), 0.98)

  tau_dr <- (m1_hat - m0_hat) +
    R_ev / G_ev * (
      (A_ev / e_ev) * (Y_ev - m1_hat) -
        ((1 - A_ev) / (1 - e_ev)) * (Y_ev - m0_hat)
    )

  data.frame(
    A = A_ev,
    Y = Y_ev,
    R_h = R_ev,
    G = G_ev,
    e_hat = e_ev,
    m0_hat = m0_hat,
    m1_hat = m1_hat,
    tau_hat = tau_hat,
    tau_dr = tau_dr
  )
}

policy_value <- function(dat, threshold) {
  if (is.infinite(threshold) && threshold < 0) {
    a_hat <- rep(1L, nrow(dat))
  } else if (is.infinite(threshold) && threshold > 0) {
    a_hat <- rep(0L, nrow(dat))
  } else {
    a_hat <- as.integer(dat$tau_hat >= threshold)
  }
  m_hat <- ifelse(a_hat == 1L, dat$m1_hat, dat$m0_hat)
  m_A <- ifelse(dat$A == 1L, dat$m1_hat, dat$m0_hat)
  p_A <- ifelse(dat$A == 1L, dat$e_hat, 1 - dat$e_hat)
  aug <- ifelse(
    dat$R_h == 1L & dat$A == a_hat,
    (1 / dat$G) / p_A * (dat$Y - m_A),
    0
  )
  mean(m_hat + aug)
}

centered_auqc <- function(dat) {
  ord <- order(dat$tau_hat, decreasing = TRUE, na.last = NA)
  n <- length(ord)
  fraction <- seq_len(n) / n
  uplift <- cumsum(dat$tau_dr[ord]) / n
  random_uplift <- fraction * mean(dat$tau_dr[ord])
  centered <- uplift - random_uplift
  fraction0 <- c(0, fraction)
  centered0 <- c(0, centered)
  sum(
    diff(fraction0) *
      (centered0[-1] + centered0[-length(centered0)]) / 2
  )
}

operating_metrics <- function(dat, threshold) {
  if (is.infinite(threshold) && threshold < 0) {
    a_hat <- rep(1L, nrow(dat))
  } else if (is.infinite(threshold) && threshold > 0) {
    a_hat <- rep(0L, nrow(dat))
  } else {
    a_hat <- as.integer(dat$tau_hat >= threshold)
  }
  treated <- a_hat == 1L
  harm_sur <- as.integer(dat$tau_dr < BENEFIT_MARGIN)
  bene_sur <- as.integer(dat$tau_dr >= BENEFIT_MARGIN)
  c(
    treat_fraction = mean(treated),
    surrogate_harm = ifelse(sum(treated) == 0L, NA_real_, mean(harm_sur[treated])),
    surrogate_benefit_capture = mean(bene_sur * a_hat),
    surrogate_ppv_benefit = ifelse(sum(treated) == 0L, NA_real_, mean(bene_sur[treated]))
  )
}

make_split <- function(seed) {
  set.seed(seed)
  split_id <- rep(NA_character_, nrow(dat_all))
  for (a in sort(unique(dat_all$A))) {
    idx <- which(dat_all$A == a)
    idx <- sample(idx)
    n_a <- length(idx)
    n_tr <- floor(0.50 * n_a)
    n_tu <- floor(0.25 * n_a)
    split_id[idx[seq_len(n_tr)]] <- "train"
    split_id[idx[(n_tr + 1):(n_tr + n_tu)]] <- "tune"
    split_id[idx[(n_tr + n_tu + 1):n_a]] <- "test"
  }
  list(
    train = which(split_id == "train"),
    tune = which(split_id == "tune"),
    test = which(split_id == "test")
  )
}

validate_one <- function(b) {
  seed <- 175000 + b
  split <- make_split(seed)
  model_bundle <- fit_split_models(split$train, seed + 10000)
  tune <- predict_split(model_bundle, split$tune)
  test <- predict_split(model_bundle, split$test)

  t_grid <- unique(as.numeric(quantile(
    tune$tau_hat,
    probs = seq(0.05, 0.95, by = 0.05),
    na.rm = TRUE,
    names = FALSE
  )))
  tune_values <- sapply(t_grid, function(t) policy_value(tune, t))
  threshold <- t_grid[which.max(tune_values)]
  tune_threshold_value <- max(tune_values)

  tune_all <- policy_value(tune, -Inf)
  tune_none <- policy_value(tune, Inf)
  fixed_rule <- ifelse(tune_all >= tune_none, "treat all", "treat none")

  test_threshold_value <- policy_value(test, threshold)
  test_all <- policy_value(test, -Inf)
  test_none <- policy_value(test, Inf)
  test_fixed <- if (fixed_rule == "treat all") test_all else test_none
  test_value_gain <- test_threshold_value - test_fixed
  selected_ops <- operating_metrics(test, threshold)

  tune_np <- t(sapply(t_grid, function(t) operating_metrics(tune, t)))
  ok <- which(
    !is.na(tune_np[, "surrogate_harm"]) &
      tune_np[, "surrogate_harm"] <= ALPHA_HARM
  )
  np_feasible_tune <- length(ok) > 0
  if (np_feasible_tune) {
    np_idx <- ok[which.max(tune_values[ok])]
  } else {
    np_idx <- which.min(tune_np[, "surrogate_harm"])
  }
  np_threshold <- t_grid[np_idx]
  np_test_ops <- operating_metrics(test, np_threshold)

  data.frame(
    split = b,
    n_train = length(split$train),
    n_tune = length(split$tune),
    n_test = length(split$test),
    threshold = threshold,
    fixed_rule = fixed_rule,
    tune_threshold_value = tune_threshold_value,
    tune_treat_all_value = tune_all,
    tune_treat_none_value = tune_none,
    test_threshold_value = test_threshold_value,
    test_treat_all_value = test_all,
    test_treat_none_value = test_none,
    test_fixed_value = test_fixed,
    test_value_gain_vs_selected_fixed = test_value_gain,
    test_centered_auqc = centered_auqc(test),
    test_treat_fraction = selected_ops["treat_fraction"],
    test_surrogate_harm = selected_ops["surrogate_harm"],
    test_surrogate_benefit_capture = selected_ops["surrogate_benefit_capture"],
    test_surrogate_ppv_benefit = selected_ops["surrogate_ppv_benefit"],
    np_feasible_tune = np_feasible_tune,
    np_threshold = np_threshold,
    np_test_treat_fraction = np_test_ops["treat_fraction"],
    np_test_surrogate_harm = np_test_ops["surrogate_harm"],
    np_test_surrogate_benefit_capture = np_test_ops["surrogate_benefit_capture"],
    np_test_surrogate_ppv_benefit = np_test_ops["surrogate_ppv_benefit"]
  )
}

summarise_metric <- function(x) {
  x <- x[is.finite(x)]
  data.frame(
    mean = mean(x),
    sd = stats::sd(x),
    p025 = as.numeric(stats::quantile(x, 0.025, names = FALSE)),
    median = stats::median(x),
    p975 = as.numeric(stats::quantile(x, 0.975, names = FALSE))
  )
}

cat(sprintf(
  "[ACTG validation] B=%d repeated splits, trees=%d, threads=%d\n",
  B, NUM_TREES, NUM_THREADS
))
cat(sprintf(
  "[ACTG validation] n=%d, observed outcomes=%d, early censored=%d\n",
  nrow(dat_all), sum(dat_all$R_h), sum(dat_all$early_censored)
))

st1 <- stage1_intervals()
stage1_path <- file.path(OUT_DIR, "actg175_ipcw_stage1_effect_intervals.csv")
write.csv(st1$option_c, stage1_path, row.names = FALSE)
stage1_global_path <- file.path(OUT_DIR, "actg175_ipcw_stage1_global.csv")
write.csv(st1$global, stage1_global_path, row.names = FALSE)
stage1_boot_path <- file.path(OUT_DIR, "actg175_ipcw_stage1_bootstrap_coefficients.csv")
write.csv(st1$bootstrap_coefficients, stage1_boot_path, row.names = FALSE)
stage1_cov_path <- file.path(OUT_DIR, "actg175_ipcw_stage1_bootstrap_covariance.csv")
write.csv(st1$bootstrap_covariance, stage1_cov_path)

split_path <- file.path(OUT_DIR, "actg175_ipcw_validation_splits.csv")
if (file.exists(split_path)) file.remove(split_path)

res_list <- vector("list", B)
for (b in seq_len(B)) {
  if (b == 1L || b %% 10L == 0L) {
    cat(sprintf("[ACTG validation] split %d/%d\n", b, B))
  }
  res_list[[b]] <- validate_one(b)
  if (b %% 10L == 0L || b == B) {
    write.csv(bind_rows(res_list[seq_len(b)]), split_path, row.names = FALSE)
  }
}

res <- bind_rows(res_list)
write.csv(res, split_path, row.names = FALSE)

summary_metrics <- bind_rows(lapply(
  c(
    "threshold",
    "test_threshold_value",
    "test_treat_all_value",
    "test_treat_none_value",
    "test_fixed_value",
    "test_value_gain_vs_selected_fixed",
    "test_centered_auqc",
    "test_treat_fraction",
    "test_surrogate_harm",
    "test_surrogate_benefit_capture",
    "test_surrogate_ppv_benefit",
    "np_test_treat_fraction",
    "np_test_surrogate_harm",
    "np_test_surrogate_benefit_capture",
    "np_test_surrogate_ppv_benefit"
  ),
  function(metric) cbind(metric = metric, summarise_metric(res[[metric]]))
))

summary_extra <- data.frame(
  metric = c(
    "B",
    "num_trees",
    "num_threads",
    "analysis_n",
    "analysis_observed_outcomes",
    "analysis_early_censored",
    "stage1_global_chisq",
    "stage1_global_df",
    "stage1_global_p",
    "fixed_rule_treat_all_rate",
    "value_gain_positive_rate",
    "np_feasible_tune_rate"
  ),
  mean = c(
    B,
    NUM_TREES,
    NUM_THREADS,
    nrow(dat_all),
    sum(dat_all$R_h),
    sum(dat_all$early_censored),
    st1$global$chisq,
    st1$global$df,
    st1$global$p_value,
    mean(res$fixed_rule == "treat all"),
    mean(res$test_value_gain_vs_selected_fixed > 0),
    mean(res$np_feasible_tune)
  ),
  sd = NA_real_,
  p025 = NA_real_,
  median = NA_real_,
  p975 = NA_real_
)

summary_df <- bind_rows(summary_metrics, summary_extra)
summary_path <- file.path(OUT_DIR, "actg175_ipcw_validation_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

session_path <- file.path(OUT_DIR, "actg175_ipcw_validation_session_info.txt")
capture.output(
  sessionInfo(),
  file = session_path
)

cat("\n[ACTG validation summary]\n")
print(summary_df)
cat("\n[Stage 1 prespecified interaction intervals]\n")
print(st1$option_c)
cat("\n[Saved]\n")
cat("- ", split_path, "\n", sep = "")
cat("- ", summary_path, "\n", sep = "")
cat("- ", stage1_path, "\n", sep = "")
cat("- ", stage1_global_path, "\n", sep = "")
cat("- ", stage1_boot_path, "\n", sep = "")
cat("- ", stage1_cov_path, "\n", sep = "")
cat("- ", session_path, "\n", sep = "")
