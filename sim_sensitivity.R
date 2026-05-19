setwd("~/Library/CloudStorage/Box-Box/Xi/HT_vs_ML")
set.seed(2025)

suppressPackageStartupMessages({
  library(grf)
  library(ggplot2)
  library(pROC)
  library(PRROC)
  library(dplyr)
  library(scales)
})

logit  <- function(p) log(p / (1 - p))
expit  <- function(z) 1 / (1 + exp(-z))

## ---- 0-helpers-uplift ----------------------------------------------------
uplift_curve_from_DR <- function(score, W_DR) {
  ord <- order(score, decreasing = TRUE)
  cum_gain <- cumsum(W_DR[ord])
  n <- length(score)
  frac <- seq_len(n) / n
  random_ref <- frac * sum(W_DR)
  centered_gain <- cum_gain - random_ref
  data.frame(
    frac = frac,
    cum_gain = cum_gain,
    random_ref = random_ref,
    centered_gain = centered_gain
  )
}

auqc_from_curve <- function(frac, gain) {
  frac0 <- c(0, frac)
  gain0 <- c(0, gain)
  sum(diff(frac0) * (gain0[-1] + gain0[-length(gain0)]) / 2) / length(frac)
}

auqc_from_DR <- function(score, W_DR) {
  curve <- uplift_curve_from_DR(score, W_DR)
  c(
    raw = auqc_from_curve(curve$frac, curve$cum_gain),
    centered = auqc_from_curve(curve$frac, curve$centered_gain)
  )
}

## ---- 0-helpers-stepp -----------------------------------------------------
stepp_df <- function(X1, A, Y, w = 150L, step = 20L) {
  ord <- order(X1); X1o <- X1[ord]; Ao <- A[ord]; Yo <- Y[ord]
  idxs <- lapply(seq(1, length(X1) - w + 1, by = step), function(i) i:(i + w - 1))
  
  te_hat <- sapply(idxs, function(id) {
    yt <- Yo[id][Ao[id] == 1]
    yc <- Yo[id][Ao[id] == 0]
    if (length(yt) < 10 || length(yc) < 10) return(NA_real_)
    mean(yt) - mean(yc)
  })
  x_mid <- sapply(idxs, function(id) mean(X1o[id]))
  
  data.frame(x_mid = x_mid, te_hat = te_hat)
}

## ---- 0-helpers-dr-value --------------------------------------------------
dr_value_vec <- function(t_grid, score, A, Y, mu1, mu0, ehat) {
  sapply(t_grid, function(t) {
    pi_t <- as.integer(score >= t)       # policy: treat if score >= t
    m_hat <- ifelse(pi_t == 1, mu1, mu0) # outcome model at recommended arm
    m_A   <- ifelse(A == 1, mu1, mu0)    # outcome model at observed arm
    w_ipw <- ifelse(A == 1, 1 / pmax(ehat, 1e-4),
                    1 / pmax(1 - ehat, 1e-4))
    aug   <- ifelse(A == pi_t, w_ipw * (Y - m_A), 0)
    mean(m_hat + aug)
  })
}

## DR value for an arbitrary (possibly non-threshold) policy
dr_value_policy <- function(a_hat, A, Y, mu1, mu0, ehat) {
  stopifnot(length(a_hat) == length(A))
  a_hat <- as.integer(a_hat)
  
  # m_hat(X, a_hat)
  m_hat <- ifelse(a_hat == 1L, mu1, mu0)
  
  # DR augmentation term
  w <- ifelse(A == 1L,
              1 / pmax(ehat, 1e-4),
              1 / pmax(1 - ehat, 1e-4))
  
  aug <- ifelse(
    A == a_hat,
    w * (Y - ifelse(A == 1L, mu1, mu0)),
    0
  )
  
  mean(m_hat + aug)
}


## ---- 0-helpers-np-single -------------------------------------------------
np_choose <- function(score, W_DR, delta = 0.03,
                      alpha_harm = 0.10, conf = 0.95,
                      t_grid = NULL, min_treat = 100,
                      ci_method = c("clopper-pearson", "wilson")) {
  
  ci_method <- match.arg(ci_method)
  n <- length(score)
  if (is.null(t_grid)) t_grid <- seq(min(score), max(score), length.out = 200)
  
  harm_upper <- det_rate <- rep(NA_real_, length(t_grid))
  n_treat    <- integer(length(t_grid))
  
  best <- list(t = NA_real_, det = -Inf, harm_upper = Inf, feasible = FALSE)
  
  # Surrogate harm / benefit at individual level
  harm_sur <- as.integer(W_DR <= 0)
  bene_sur <- as.integer(W_DR >  delta)
  
  for (j in seq_along(t_grid)) {
    t <- t_grid[j]
    treat <- as.integer(score >= t)
    n_t <- sum(treat == 1)
    n_treat[j] <- n_t
    if (n_t < min_treat) next
    
    m_harm <- sum(treat == 1 & harm_sur == 1)
    m_det  <- sum(treat == 1 & bene_sur == 1)
    
    # one-sided upper bound on harm among treated
    if (ci_method == "clopper-pearson") {
      upper <- binom.test(m_harm, n_t, conf.level = conf,
                          alternative = "less")$conf.int[2]
    } else {
      pt <- prop.test(m_harm, n_t, conf.level = conf,
                      alternative = "less", correct = FALSE)
      upper <- pt$conf.int[2]
    }
    
    harm_upper[j] <- upper
    det_rate[j]   <- m_det / n
    
    if (upper <= alpha_harm && det_rate[j] > best$det) {
      best <- list(t = t, det = det_rate[j],
                   harm_upper = upper, feasible = TRUE)
    }
  }
  
  # Fallback if no threshold satisfies the bound
  if (!best$feasible) {
    j_min <- which.min(harm_upper)
    best$t          <- t_grid[j_min]
    best$det        <- det_rate[j_min]
    best$harm_upper <- harm_upper[j_min]
  }
  
  list(
    t_star     = best$t,
    det        = best$det,
    harm_upper = best$harm_upper,
    feasible   = best$feasible,
    curve      = data.frame(
      t = t_grid,
      n_treat = n_treat,
      harm_upper = harm_upper,
      det = det_rate
    )
  )
}

## ---- 0-theme -------------------------------------------------------------
theme_pub <- function(base_size = 10.5, base_family = NULL) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      axis.text        = element_text(color = "black"),
      axis.ticks       = element_line(color = "black"),
      plot.title       = element_text(face = "bold", hjust = 0, size = 12),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = element_blank()
    )
}

## ---- 1-dgp-generic -------------------------------------------------------
baseline_risk <- function(X) {
  0.25 + 0.35 * expit(-0.4 + 0.45 * X$X1 - 0.25 * X$X2 + 0.35 * X$X3)
}

scenario_tau <- function(X, scenario) {
  switch(
    scenario,
    "No treatment effect" = rep(0, nrow(X)),
    "No HTE" = rep(0.04, nrow(X)),
    "Weak quantitative HTE" = 0.04 + 0.015 * tanh(X$X1 / 1.2),
    "Strong quantitative HTE" = 0.06 + 0.050 * tanh(X$X1),
    "Weak qualitative HTE" = 0.080 * tanh(X$X1 / 1.1),
    "Strong qualitative HTE" = 0.180 * tanh(X$X1),
    stop("Unknown simulation scenario: ", scenario)
  )
}

sim_one_generic <- function(n = 2000, delta = 0.03,
                            scenario = "No HTE") {
  # Covariates
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)
  X  <- data.frame(X1 = X1, X2 = X2, X3 = X3)
  
  # Treatment (randomized 1:1)
  A  <- rbinom(n, 1, 0.5)
  
  # Risks and treatment effect are defined on the absolute benefit scale.
  p0 <- baseline_risk(X)
  tau_true <- scenario_tau(X, scenario)
  p1 <- p0 + tau_true
  if (any(p1 < 0 | p1 > 1)) {
    stop("Invalid treated risk outside [0, 1] in scenario: ", scenario)
  }
  
  # Potential outcomes and observed outcome
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y  <- ifelse(A == 1, Y1, Y0)
  
  # True CATE and benefiter indicator on the risk-difference scale
  Z_true   <- as.integer(tau_true > delta)
  
  list(
    X = X, X1 = X1, scenario = scenario,
    A = A, Y = Y,
    p0 = p0, p1 = p1,
    tau_true = tau_true,
    Z_true   = Z_true,
    delta    = delta
  )
}

## ---- 1-dgp-no-hte --------------------------------------------------------
sim_no_treatment_effect <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "No treatment effect")
}

sim_no_hte <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "No HTE")
}

## ---- 1-dgp-weak-hte ------------------------------------------------------
sim_weak_quant_hte <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "Weak quantitative HTE")
}

## ---- 1-dgp-strong-hte ----------------------------------------------------
sim_strong_quant_hte <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "Strong quantitative HTE")
}

sim_weak_qual_hte <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "Weak qualitative HTE")
}

sim_strong_qual_hte <- function(n = 2000, delta = 0.03) {
  sim_one_generic(n = n, delta = delta, scenario = "Strong qualitative HTE")
}

sim_weak_hte <- sim_weak_quant_hte
sim_strong_hte <- sim_strong_qual_hte

truth_metrics <- function(dat) {
  value_none <- mean(dat$p0)
  value_all <- mean(dat$p1)
  fixed_best <- max(value_none, value_all)
  treat_oracle <- dat$tau_true > 0
  value_oracle <- mean(ifelse(treat_oracle, dat$p1, dat$p0))
  
  data.frame(
    truth_ate = mean(dat$tau_true),
    truth_tau_sd = sd(dat$tau_true),
    truth_benefit_rate = mean(dat$tau_true > dat$delta),
    truth_harm_rate = mean(dat$tau_true < 0),
    truth_value_treat_none = value_none,
    truth_value_treat_all = value_all,
    truth_value_oracle = value_oracle,
    truth_oracle_value_gain = value_oracle - fixed_best,
    truth_oracle_treat_fraction = mean(treat_oracle)
  )
}

## ---- 2-stage1 ------------------------------------------------------------
run_stage1 <- function(dat) {
  X <- dat$X
  A <- dat$A
  Y <- dat$Y
  
  df_glm <- cbind(X, A = A, Y = Y)
  
  # Risk-difference scale models. Stage 1 tests treatment-covariate
  # interactions on the same absolute benefit scale used by the DGP.
  fit_full <- lm(Y ~ A * (X1 + X2 + X3), data = df_glm)
  
  # Option A: omnibus robust Wald test for all A:X interactions
  coefs <- coef(fit_full)
  idx_ax <- grep("^A:", names(coefs))
  beta_ax <- coefs[idx_ax]
  vc <- sandwich::vcovHC(fit_full, type = "HC3")[idx_ax, idx_ax, drop = FALSE]
  stat_global <- as.numeric(t(beta_ax) %*% solve(vc, beta_ax))
  p_global <- pchisq(stat_global, df = length(beta_ax), lower.tail = FALSE)
  
  # Option C: prespecified interactions with Holm adjustment
  se_ax <- sqrt(diag(vc))
  p_int <- 2 * pnorm(abs(beta_ax / se_ax), lower.tail = FALSE)
  names(p_int) <- sub("^A:", "", names(beta_ax))
  p_holm <- p.adjust(p_int, method = "holm")
  
  # Confirmatory Stage 1 gate: proceed only if the omnibus Option A test rejects.
  # Holm-adjusted Option C tests are reported to localize the signal after the gate opens.
  proceed <- (p_global < 0.05)
  
  # STEPP-like windowed risk difference along X1
  df_stepp <- stepp_df(X1 = dat$X1, A = A, Y = Y)
  
  list(
    p_global = p_global,
    p_int    = p_int,
    p_holm   = p_holm,
    proceed  = proceed,
    df_stepp = df_stepp
  )
}

## ---- 3-stage2 ------------------------------------------------------------
make_stage2_splits <- function(n, train_frac = 0.50, tune_frac = 0.25) {
  idx <- sample(seq_len(n))
  n_train <- floor(train_frac * n)
  n_tune <- floor(tune_frac * n)
  
  list(
    train = sort(idx[seq_len(n_train)]),
    tune = sort(idx[(n_train + 1):(n_train + n_tune)]),
    test = sort(idx[(n_train + n_tune + 1):n])
  )
}

threshold_grid <- function(score, n_grid = 200) {
  rng <- range(score, na.rm = TRUE)
  span <- diff(rng)
  pad <- max(span * 1e-4, 1e-6)
  seq(rng[1] - pad, rng[2] + pad, length.out = n_grid)
}

dr_pseudo <- function(A, Y, mu1, mu0, ehat) {
  ehat <- pmin(pmax(ehat, 1e-3), 1 - 1e-3)
  mu1 - mu0 +
    A * (Y - mu1) / ehat -
    (1 - A) * (Y - mu0) / (1 - ehat)
}

fit_nuisance_models <- function(X_train, A_train, Y_train) {
  df_train <- data.frame(Y = Y_train, A = A_train, X_train)
  list(
    e_fit = glm(A ~ X1 + X2 + X3, family = binomial(), data = df_train),
    m1_fit = glm(Y ~ X1 + X2 + X3, family = binomial(),
                 data = df_train[df_train$A == 1, ]),
    m0_fit = glm(Y ~ X1 + X2 + X3, family = binomial(),
                 data = df_train[df_train$A == 0, ])
  )
}

predict_nuisance <- function(fits, X_new) {
  ehat <- predict(fits$e_fit, newdata = X_new, type = "response")
  list(
    ehat = pmin(pmax(ehat, 1e-3), 1 - 1e-3),
    mu1 = predict(fits$m1_fit, newdata = X_new, type = "response"),
    mu0 = predict(fits$m0_fit, newdata = X_new, type = "response")
  )
}

np_eval_threshold <- function(score, W_DR, threshold, delta = 0.03,
                              conf = 0.95,
                              ci_method = c("clopper-pearson", "wilson")) {
  ci_method <- match.arg(ci_method)
  treat <- as.integer(score >= threshold)
  n_t <- sum(treat == 1)
  harm_sur <- as.integer(W_DR <= 0)
  bene_sur <- as.integer(W_DR > delta)
  m_harm <- sum(treat == 1 & harm_sur == 1)
  m_det <- sum(treat == 1 & bene_sur == 1)
  
  if (n_t == 0) {
    return(list(harm_upper = NA_real_, det = 0, n_treat = 0))
  }
  
  if (ci_method == "clopper-pearson") {
    upper <- binom.test(m_harm, n_t, conf.level = conf,
                        alternative = "less")$conf.int[2]
  } else {
    upper <- prop.test(m_harm, n_t, conf.level = conf,
                       alternative = "less", correct = FALSE)$conf.int[2]
  }
  
  list(harm_upper = upper, det = m_det / length(score), n_treat = n_t)
}

evaluate_stage2_score <- function(score_tune, score_test,
                                  A_tune, Y_tune, mu1_tune, mu0_tune, ehat_tune,
                                  A_test, Y_test, mu1_test, mu0_test, ehat_test,
                                  Z_test, p0_test, p1_test, delta,
                                  alpha_harm = 0.10,
                                  ci_method = "wilson") {
  W_tune <- dr_pseudo(A_tune, Y_tune, mu1_tune, mu0_tune, ehat_tune)
  W_test <- dr_pseudo(A_test, Y_test, mu1_test, mu0_test, ehat_test)
  
  auqc <- auqc_from_DR(score_test, W_test)
  
  if (length(unique(Z_test)) == 2L) {
    roc_true <- roc(Z_test, score_test, quiet = TRUE)
    auroc_true <- as.numeric(auc(roc_true))
    
    pr_true <- pr.curve(
      scores.class0 = score_test[Z_test == 1],
      scores.class1 = score_test[Z_test == 0],
      curve = FALSE
    )
    auprc_true <- pr_true$auc.integral
  } else {
    auroc_true <- NA_real_
    auprc_true <- NA_real_
  }
  
  t_grid <- threshold_grid(score_tune)
  val_grid_tune <- dr_value_vec(t_grid, score_tune, A_tune, Y_tune,
                                mu1_tune, mu0_tune, ehat_tune)
  t_star_value <- t_grid[which.max(val_grid_tune)]
  
  a_star <- as.integer(score_test >= t_star_value)
  val_star <- dr_value_policy(a_star, A_test, Y_test, mu1_test, mu0_test, ehat_test)
  
  val_all1_tune <- dr_value_policy(rep(1L, length(Y_tune)), A_tune, Y_tune,
                                   mu1_tune, mu0_tune, ehat_tune)
  val_all0_tune <- dr_value_policy(rep(0L, length(Y_tune)), A_tune, Y_tune,
                                   mu1_tune, mu0_tune, ehat_tune)
  selected_const_treat <- as.integer(val_all1_tune >= val_all0_tune)
  
  a_all1 <- rep(1L, length(Y_test))
  a_all0 <- rep(0L, length(Y_test))
  val_all1 <- dr_value_policy(a_all1, A_test, Y_test, mu1_test, mu0_test, ehat_test)
  val_all0 <- dr_value_policy(a_all0, A_test, Y_test, mu1_test, mu0_test, ehat_test)
  val_const_selected <- if (selected_const_treat == 1L) val_all1 else val_all0
  value_gain_all <- val_star - val_const_selected
  
  truth_value_star <- mean(ifelse(a_star == 1, p1_test, p0_test))
  truth_fixed_best <- max(mean(p1_test), mean(p0_test))
  truth_selected_const <- if (selected_const_treat == 1L) mean(p1_test) else mean(p0_test)
  truth_value_gain_learned <- truth_value_star - truth_fixed_best
  truth_value_gain_selected_const <- truth_value_star - truth_selected_const
  
  np_tune <- np_choose(
    score = score_tune,
    W_DR = W_tune,
    delta = delta,
    alpha_harm = alpha_harm,
    conf = 0.95,
    t_grid = t_grid,
    min_treat = max(25, 0.05 * length(score_tune)),
    ci_method = ci_method
  )
  np_test <- np_eval_threshold(
    score = score_test,
    W_DR = W_test,
    threshold = np_tune$t_star,
    delta = delta,
    conf = 0.95,
    ci_method = ci_method
  )
  np_tune$harm_upper_tune <- np_tune$harm_upper
  np_tune$det_tune <- np_tune$det
  np_tune$harm_upper <- np_test$harm_upper
  np_tune$det <- np_test$det
  np_tune$n_treat_test <- np_test$n_treat
  
  list(
    tau_hat = score_test,
    W_DR = W_test,
    AUQC = unname(auqc["centered"]),
    AUQC_raw = unname(auqc["raw"]),
    t_grid = t_grid,
    val_grid = val_grid_tune,
    t_star_value = t_star_value,
    value_max = val_star,
    value_all = val_all1,
    value_none = val_all0,
    np_out = np_tune,
    auroc_true = auroc_true,
    auprc_true = auprc_true,
    value_gain_all = value_gain_all,
    truth_value_gain_learned = truth_value_gain_learned,
    truth_value_gain_selected_const = truth_value_gain_selected_const,
    selected_const_treat = selected_const_treat,
    learned_treat_fraction = mean(a_star)
  )
}

run_stage2 <- function(dat, train_frac = 0.50, tune_frac = 0.25,
                       alpha_harm = 0.10,
                       ci_method = "wilson") {
  
  X <- dat$X
  A <- dat$A
  Y <- dat$Y
  n <- nrow(X)
  
  idx <- make_stage2_splits(n, train_frac = train_frac, tune_frac = tune_frac)
  
  X_train <- X[idx$train, , drop = FALSE]
  X_tune <- X[idx$tune, , drop = FALSE]
  X_test <- X[idx$test, , drop = FALSE]
  A_train <- A[idx$train]
  Y_train <- Y[idx$train]
  A_tune <- A[idx$tune]
  Y_tune <- Y[idx$tune]
  A_test <- A[idx$test]
  Y_test <- Y[idx$test]
  
  ## 3.1 CATE learner: train, tune the threshold, then evaluate on final test patients
  cf <- causal_forest(X_train, Y_train, A_train)
  tau_hat_tune <- predict(cf, X_tune)$predictions
  tau_hat_test <- predict(cf, X_test)$predictions
  
  ## 3.2 Nuisance predictions for tuning and final testing
  nuisance <- fit_nuisance_models(X_train, A_train, Y_train)
  tune_nuis <- predict_nuisance(nuisance, X_tune)
  test_nuis <- predict_nuisance(nuisance, X_test)
  
  out <- evaluate_stage2_score(
    score_tune = tau_hat_tune,
    score_test = tau_hat_test,
    A_tune = A_tune,
    Y_tune = Y_tune,
    mu1_tune = tune_nuis$mu1,
    mu0_tune = tune_nuis$mu0,
    ehat_tune = tune_nuis$ehat,
    A_test = A_test,
    Y_test = Y_test,
    mu1_test = test_nuis$mu1,
    mu0_test = test_nuis$mu0,
    ehat_test = test_nuis$ehat,
    Z_test = dat$Z_true[idx$test],
    p0_test = dat$p0[idx$test],
    p1_test = dat$p1[idx$test],
    delta = dat$delta,
    alpha_harm = alpha_harm,
    ci_method = ci_method
  )
  out$eval_index <- idx$test
  out$tune_index <- idx$tune
  out$train_index <- idx$train
  out$ehat <- test_nuis$ehat
  out$mu1 <- test_nuis$mu1
  out$mu0 <- test_nuis$mu0
  out
}

## ---- 3-stage2-xlearner ----------------------------------------------------
# Alternative Stage 2 using a DR X-learner instead of a causal forest CATE
# The interface and outputs are kept parallel to run_stage2() so that
# downstream code can compare them easily.
run_stage2_xlearner <- function(dat, train_frac = 0.50, tune_frac = 0.25,
                                alpha_harm = 0.10,
                                ci_method = "wilson") {
  
  X <- dat$X
  A <- dat$A
  Y <- dat$Y
  n <- nrow(X)
  
  idx <- make_stage2_splits(n, train_frac = train_frac, tune_frac = tune_frac)
  
  X_train <- X[idx$train, , drop = FALSE]
  X_tune <- X[idx$tune, , drop = FALSE]
  X_test <- X[idx$test, , drop = FALSE]
  A_train <- A[idx$train]
  Y_train <- Y[idx$train]
  A_tune <- A[idx$tune]
  Y_tune <- Y[idx$tune]
  A_test <- A[idx$test]
  Y_test <- Y[idx$test]
  
  ## 3.1 Nuisance predictions for training, tuning, and final testing
  nuisance <- fit_nuisance_models(X_train, A_train, Y_train)
  train_nuis <- predict_nuisance(nuisance, X_train)
  tune_nuis <- predict_nuisance(nuisance, X_tune)
  test_nuis <- predict_nuisance(nuisance, X_test)
  
  ## 3.2 CATE learner: DR X-learner implemented with regression forests
  # Step 1: construct X-learner pseudo-outcomes for treated and controls.
  
  # For treated units (A=1): how much better than predicted control?
  D1 <- Y_train[A_train == 1] - train_nuis$mu0[A_train == 1]
  
  # For control units (A=0): how much worse than predicted treatment?
  D0 <- train_nuis$mu1[A_train == 0] - Y_train[A_train == 0]
  
  # Step 2: fit regression forests for E[D1 | X] and E[D0 | X].
  # We use grf::regression_forest to stay within your existing stack.
  
  rf1 <- regression_forest(X_train[A_train == 1, , drop = FALSE], D1)
  rf0 <- regression_forest(X_train[A_train == 0, , drop = FALSE], D0)
  
  tau1_hat_tune <- predict(rf1, X_tune)$predictions
  tau0_hat_tune <- predict(rf0, X_tune)$predictions
  tau1_hat_test <- predict(rf1, X_test)$predictions
  tau0_hat_test <- predict(rf0, X_test)$predictions
  
  # Step 3: combine the two conditional estimates using ehat(X) as weights.
  # This is the standard X-learner combination: w(X)*tau0 + (1-w(X))*tau1.
  
  tau_hat_tune <- tune_nuis$ehat * tau0_hat_tune +
    (1 - tune_nuis$ehat) * tau1_hat_tune
  tau_hat_test <- test_nuis$ehat * tau0_hat_test +
    (1 - test_nuis$ehat) * tau1_hat_test
  
  out <- evaluate_stage2_score(
    score_tune = tau_hat_tune,
    score_test = tau_hat_test,
    A_tune = A_tune,
    Y_tune = Y_tune,
    mu1_tune = tune_nuis$mu1,
    mu0_tune = tune_nuis$mu0,
    ehat_tune = tune_nuis$ehat,
    A_test = A_test,
    Y_test = Y_test,
    mu1_test = test_nuis$mu1,
    mu0_test = test_nuis$mu0,
    ehat_test = test_nuis$ehat,
    Z_test = dat$Z_true[idx$test],
    p0_test = dat$p0[idx$test],
    p1_test = dat$p1[idx$test],
    delta = dat$delta,
    alpha_harm = alpha_harm,
    ci_method = ci_method
  )
  out$eval_index <- idx$test
  out$tune_index <- idx$tune
  out$train_index <- idx$train
  out$ehat <- test_nuis$ehat
  out$mu1 <- test_nuis$mu1
  out$mu0 <- test_nuis$mu0
  out
}

## ---- 4-single-trial-visualization ----------------------------------------
plot_single_trial <- function(sim_fun,
                              n = 2000,
                              delta = 0.03,
                              alpha_harm = 0.10,
                              main_label = "") {
  
  # 4.1 Generate one dataset
  dat <- sim_fun(n = n, delta = delta)
  
  # 4.2 Run Stage 1 and Stage 2
  st1 <- run_stage1(dat)
  st2 <- run_stage2(dat, alpha_harm = alpha_harm)
  
  ## 4.3 Panel A: STEPP-style X1 exploration
  pA <- ggplot(st1$df_stepp, aes(x = x_mid, y = te_hat)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey30") +
    geom_line(linewidth = 0.55, color = "black") +
    labs(
      x = expression(X[1] ~ "(ordered windows)"),
      y = "Windowed risk difference\n(treated – control)",
      title = "A"
    ) +
    coord_cartesian(xlim = c(-2, 2)) +
    theme_pub()
  
  ## 4.4 Panel B: Uplift curve
  df_uplift <- uplift_curve_from_DR(st2$tau_hat, st2$W_DR)
  
  pB <- ggplot(df_uplift, aes(frac, centered_gain)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey40") +
    geom_line(linewidth = 0.55, color = "black") +
    annotate(
      "label",
      x = 0.98, y = max(df_uplift$centered_gain, na.rm = TRUE) * 0.98,
      label = sprintf("cAUQC = %.3f", st2$AUQC),
      hjust = 1, vjust = 1, size = 3
    ) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    labs(
      x = "Top fraction by CATE",
      y = "Centered cumulative DR uplift",
      title = "B"
    ) +
    theme_pub()
  
  ## 4.5 Panel C: Policy value vs threshold
  df_value <- data.frame(
    threshold = st2$t_grid,
    value     = st2$val_grid
  )
  
  pC <- ggplot(df_value, aes(threshold, value)) +
    geom_line(linewidth = 0.55, colour = "black") +
    geom_vline(xintercept = st2$t_star_value,
               linetype = 2, colour = "grey30") +
    annotate(
      "label",
      x = max(df_value$threshold, na.rm = TRUE),
      y = max(df_value$value,    na.rm = TRUE),
      hjust = 1, vjust = 1,
      label = sprintf("t* = %.3f", st2$t_star_value),
      size = 3
    ) +
    labs(
      x = "Threshold on CATE score",
      y = "Estimated policy value (DR)",
      title = "C"
    ) +
    theme_pub()
  
  ## 4.6 Panel D: NP-ROC curve (harm vs benefit-capture)
  df_np <- na.omit(st2$np_out$curve)
  t_np  <- st2$np_out$t_star
  pt_np <- df_np[which.min(abs(df_np$t - t_np)), , drop = FALSE]
  
  pD <- ggplot(df_np, aes(harm_upper, det)) +
    annotate("rect", xmin = -Inf, xmax = alpha_harm,
             ymin = -Inf, ymax = Inf, fill = "grey95", color = NA) +
    geom_vline(xintercept = alpha_harm,
               linetype = 2, color = "grey30") +
    geom_path(linewidth = 0.55, color = "black") +
    geom_point(data = pt_np, colour = "red", size = 2) +
    annotate(
      "text",
      x = pt_np$harm_upper, y = pt_np$det,
      label = "Best attainable",
      vjust = -1, hjust = 1, size = 3, colour = "red"
    ) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(
      x = "Harm rate",
      y = "Benefit-capture rate",
      title = "D"
    ) +
    theme_pub()
  
  list(dat = dat, st1 = st1, st2 = st2,
       pA = pA, pB = pB, pC = pC, pD = pD)
}

## ---- 5-mc-wrapper --------------------------------------------------------
run_mc_scenario <- function(sim_fun,
                            n = 2000,
                            delta = 0.03,
                            B = 500,
                            alpha_harm = 0.10,
                            scenario_label = "scenario") {
  
  res_list <- vector("list", B)
  
  for (b in 1:B) {
    dat <- sim_fun(n = n, delta = delta)
    
    st1 <- run_stage1(dat)
    st2 <- run_stage2(dat, alpha_harm = alpha_harm)
    
    res_list[[b]] <- data.frame(
      scenario       = scenario_label,
      p_global       = st1$p_global,
      min_p_holm     = min(st1$p_holm),
      proceed        = as.integer(st1$proceed),
      AUQC           = st2$AUQC,
      AUQC_raw       = st2$AUQC_raw,
      value_gain_all = st2$value_gain_all,
      NP_feasible    = as.integer(st2$np_out$feasible),
      NP_harm        = st2$np_out$harm_upper,
      NP_bene        = st2$np_out$det,
      AUROC_true     = st2$auroc_true,
      AUPRC_true     = st2$auprc_true
    )
  }
  
  do.call(rbind, res_list)
}

## ---- 5-mc-wrapper-sensitivity ---------------------------------------------
# Monte Carlo wrapper that runs Stage 2 with two learners:
#   - causal forest (original run_stage2)
#   - DR X-learner (run_stage2_xlearner)
#
# Stage 1 is shared; we record the same gate information for both.
run_mc_scenario_sens <- function(sim_fun,
                                 n = 2000,
                                 delta = 0.03,
                                 B = 500,
                                 alpha_harm = 0.10,
                                 scenario_label = "scenario") {
  
  res_list <- vector("list", B)
  
  for (b in 1:B) {
    if (b == 1 || b %% 25 == 0 || b == B) {
      message(sprintf("[%s] replicate %d/%d", scenario_label, b, B))
    }
    # 1. Generate one trial replicate
    dat <- sim_fun(n = n, delta = delta)
    truth <- truth_metrics(dat)
    
    # 2. Stage 1: population-level inference (shared across learners)
    st1 <- run_stage1(dat)
    
    # 3. Stage 2: learner 1 = causal forest (original)
    st2_cf <- run_stage2(dat, alpha_harm = alpha_harm)
    
    # 4. Stage 2: learner 2 = DR X-learner
    st2_x  <- run_stage2_xlearner(dat, alpha_harm = alpha_harm)
    
    # 5. Combine results into two rows (one per learner) for this replicate
    res_list[[b]] <- rbind(
      cbind(data.frame(
        scenario       = scenario_label,
        learner        = "causal_forest",
        n              = n,
        replicate      = b,
        p_global       = st1$p_global,
      min_p_holm     = min(st1$p_holm),
      proceed        = as.integer(st1$proceed),
      AUQC           = st2_cf$AUQC,
      AUQC_raw       = st2_cf$AUQC_raw,
      value_gain_all = st2_cf$value_gain_all,
        NP_feasible    = as.integer(st2_cf$np_out$feasible),
        NP_harm        = st2_cf$np_out$harm_upper,
        NP_bene        = st2_cf$np_out$det,
        AUROC_true     = st2_cf$auroc_true,
        AUPRC_true     = st2_cf$auprc_true,
        truth_value_gain_learned = st2_cf$truth_value_gain_learned,
        truth_value_gain_selected_const = st2_cf$truth_value_gain_selected_const,
        selected_const_treat = st2_cf$selected_const_treat,
        learned_treat_fraction = st2_cf$learned_treat_fraction
      ), truth),
      cbind(data.frame(
        scenario       = scenario_label,
        learner        = "x_learner",
        n              = n,
        replicate      = b,
        p_global       = st1$p_global,
      min_p_holm     = min(st1$p_holm),
      proceed        = as.integer(st1$proceed),
      AUQC           = st2_x$AUQC,
      AUQC_raw       = st2_x$AUQC_raw,
      value_gain_all = st2_x$value_gain_all,
        NP_feasible    = as.integer(st2_x$np_out$feasible),
        NP_harm        = st2_x$np_out$harm_upper,
        NP_bene        = st2_x$np_out$det,
        AUROC_true     = st2_x$auroc_true,
        AUPRC_true     = st2_x$auprc_true,
        truth_value_gain_learned = st2_x$truth_value_gain_learned,
        truth_value_gain_selected_const = st2_x$truth_value_gain_selected_const,
        selected_const_treat = st2_x$selected_const_treat,
        learned_treat_fraction = st2_x$learned_treat_fraction
      ), truth)
    )
  }
  
  # Stack all replicates and learners
  do.call(rbind, res_list)
}

## ---- 6-run-sensitivity ----------------------------------------------------
scenario_funs <- list(
  "No treatment effect" = sim_no_treatment_effect,
  "No HTE" = sim_no_hte,
  "Weak quantitative HTE" = sim_weak_quant_hte,
  "Strong quantitative HTE" = sim_strong_quant_hte,
  "Weak qualitative HTE" = sim_weak_qual_hte,
  "Strong qualitative HTE" = sim_strong_qual_hte
)

B_mc <- as.integer(Sys.getenv("SIM_B_SENS", "500"))
n_mc <- as.integer(Sys.getenv("SIM_N_SENS", "2000"))

mc_all_sens <- bind_rows(lapply(names(scenario_funs), function(label) {
  run_mc_scenario_sens(
    sim_fun = scenario_funs[[label]],
    scenario_label = label,
    n = n_mc,
    B = B_mc
  )
}))

mc_se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

# Quick summary: compare learners within each scenario
mc_sens_summary <- mc_all_sens %>%
  group_by(scenario, learner) %>%
  summarise(
    n = first(n),
    B = n(),
    proceed_rate   = mean(proceed),
    se_proceed_rate = mc_se(proceed),
    any_holm_rate  = mean(min_p_holm < 0.05),
    mean_AUQC      = mean(AUQC),
    se_AUQC = mc_se(AUQC),
    mean_AUQC_raw = mean(AUQC_raw),
    se_AUQC_raw = mc_se(AUQC_raw),
    mean_valuegain = mean(value_gain_all),
    se_valuegain = mc_se(value_gain_all),
    mean_truth_valuegain_learned = mean(truth_value_gain_learned),
    se_truth_valuegain_learned = mc_se(truth_value_gain_learned),
    mean_truth_valuegain_selected_const = mean(truth_value_gain_selected_const),
    se_truth_valuegain_selected_const = mc_se(truth_value_gain_selected_const),
    mean_selected_const_treat = mean(selected_const_treat),
    mean_learned_treat_fraction = mean(learned_treat_fraction),
    NP_feasible    = mean(NP_feasible),
    mean_NP_harm   = mean(NP_harm, na.rm = TRUE),
    se_NP_harm = mc_se(NP_harm),
    mean_NP_bene   = mean(NP_bene, na.rm = TRUE),
    se_NP_bene = mc_se(NP_bene),
    mean_AUROC_true = mean(AUROC_true, na.rm = TRUE),
    se_AUROC_true = mc_se(AUROC_true),
    mean_AUPRC_true = mean(AUPRC_true, na.rm = TRUE),
    se_AUPRC_true = mc_se(AUPRC_true),
    truth_ate = mean(truth_ate),
    truth_tau_sd = mean(truth_tau_sd),
    truth_benefit_rate = mean(truth_benefit_rate),
    truth_harm_rate = mean(truth_harm_rate),
    se_truth_oracle_value_gain = mc_se(truth_oracle_value_gain),
    truth_oracle_value_gain = mean(truth_oracle_value_gain),
    .groups = "drop"
  )

print(mc_sens_summary)

write.csv(mc_all_sens, "result/sim_sensitivity_SiM.csv", row.names = FALSE)
write.csv(mc_sens_summary, "result/sim_sensitivity_SiM_summary.csv", row.names = FALSE)

# ## Example (this is just demonstration; you can run it yourself)
# mc_no   <- run_mc_scenario(sim_no_hte,   scenario_label = "No HTE",   B = 200)
# mc_weak <- run_mc_scenario(sim_weak_hte, scenario_label = "Weak HTE", B = 200)
# mc_str  <- run_mc_scenario(sim_strong_hte, scenario_label = "Strong HTE", B = 200)
# 
# mc_all <- bind_rows(mc_no, mc_weak, mc_str)
# 
# # Gate proceed rate by scenario
# mc_all %>%
#   group_by(scenario) %>%
#   summarise(
#     proceed_rate   = mean(proceed),
#     mean_AUQC      = mean(AUQC),
#     mean_valuegain = mean(value_gain_all),
#     NP_feasible    = mean(NP_feasible)
#   )
# 
# write.csv(mc_all, 'result/mc_all.csv', row.names = F)

if (interactive()) {
  plot_single_trial(sim_strong_qual_hte, main_label = "Strong qualitative HTE: ")
  plot_single_trial(sim_no_hte, main_label = "No HTE: ")
  plot_single_trial(sim_weak_quant_hte, main_label = "Weak quantitative HTE: ")
}
