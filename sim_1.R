set.seed(2025)

logit  <- function(p) log(p / (1 - p))
expit  <- function(z) 1 / (1 + exp(-z))

auqc_from_DR <- function(score, W_DR) {
  ord <- order(score, decreasing = TRUE)
  cum_gain <- cumsum(W_DR[ord])
  n <- length(score)
  # Trapezoidal integral over fraction treated (0–1), normalized by n
  sum((cum_gain[-1] + cum_gain[-n]) / 2) / n
}

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
theme_pub <- function(base_size = 15, base_family = NULL) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", size = 0.3),
      axis.text        = element_text(color = "black"),
      axis.ticks       = element_line(color = "black"),
      plot.title       = element_text(face = "bold", hjust = 0),
      panel.border     = element_rect(color = "black", fill = NA, size = 0.6),
      strip.background = element_blank()
    )
}

## ---- 1-dgp-generic -------------------------------------------------------
sim_one_generic <- function(n = 2000, delta = 0.03,
                            gamma0, gamma1) {
  # Covariates
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)
  X  <- data.frame(X1 = X1, X2 = X2, X3 = X3)
  
  # Treatment (randomized 1:1)
  A  <- rbinom(n, 1, 0.5)
  
  # Control arm risk
  eta0 <- -0.6 + 0.6 * X1 - 0.2 * X2 + 0.3 * X3
  p0   <- expit(eta0)
  
  # Treatment log-odds increment (possibly heterogeneous)
  eta1 <- eta0 + gamma0 + gamma1 * X1
  p1   <- expit(eta1)
  
  # Potential outcomes and observed outcome
  Y0 <- rbinom(n, 1, p0)
  Y1 <- rbinom(n, 1, p1)
  Y  <- ifelse(A == 1, Y1, Y0)
  
  # True CATE and benefiter indicator
  tau_true <- p1 - p0
  Z_true   <- as.integer(tau_true > delta)
  
  list(
    X = X, X1 = X1,
    A = A, Y = Y,
    p0 = p0, p1 = p1,
    tau_true = tau_true,
    Z_true   = Z_true,
    delta    = delta
  )
}

## ---- 1-dgp-no-hte --------------------------------------------------------
sim_no_hte <- function(n = 2000, delta = 0.03) {
  # Constant treatment effect: gamma1 = 0
  # Choose gamma0 > 0 so treatment is beneficial but non-heterogeneous
  sim_one_generic(n = n, delta = delta,
                  gamma0 = 0.4,   # moderate positive effect for everyone
                  gamma1 = 0.0)   # no HTE
}

## ---- 1-dgp-weak-hte ------------------------------------------------------
sim_weak_hte <- function(n = 2000, delta = 0.03) {
  # Same baseline as strong scenario but milder heterogeneity:
  # gamma1 smaller -> tau(X) less spread, many subjects near the margin δ
  sim_one_generic(n = n, delta = delta,
                  gamma0 = -0.05,
                  gamma1 = 0.3)   # weaker slope than strong HTE
}

## ---- 1-dgp-strong-hte ----------------------------------------------------
sim_strong_hte <- function(n = 2000, delta = 0.03) {
  # Your current simulation: strong effect modification by X1
  sim_one_generic(n = n, delta = delta,
                  gamma0 = -0.05,
                  gamma1 = 1.0)
}

## ---- 2-stage1 ------------------------------------------------------------
run_stage1 <- function(dat) {
  X <- dat$X
  A <- dat$A
  Y <- dat$Y
  
  df_glm <- cbind(X, A = A, Y = Y)
  
  # Base logistic model (no interactions)
  fit_base <- glm(Y ~ A + X1 + X2 + X3,
                  family = binomial(), data = df_glm)
  
  # Full model with A×X interactions
  fit_full <- glm(Y ~ A * (X1 + X2 + X3),
                  family = binomial(), data = df_glm)
  
  # Option A: omnibus LRT for "any interaction"
  lrt_tab   <- anova(fit_base, fit_full, test = "LRT")
  p_global  <- lrt_tab[2, "Pr(>Chi)"]
  
  # Option C: prespecified interactions (here: all A:Xi terms with Holm)
  coefs  <- summary(fit_full)$coefficients
  idx_ax <- grep("^A:", rownames(coefs))
  p_int  <- coefs[idx_ax, "Pr(>|z|)"]
  names(p_int) <- sub("^A:", "", rownames(coefs)[idx_ax])
  p_holm <- p.adjust(p_int, method = "holm")
  
  # Gate: proceed to Stage 2 if global or any Holm-adjusted test is significant
  proceed <- (p_global < 0.05) || any(p_holm < 0.05, na.rm = TRUE)
  
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
run_stage2 <- function(dat, K = 5,
                       alpha_harm = 0.10,
                       ci_method = "wilson") {
  
  X <- dat$X
  A <- dat$A
  Y <- dat$Y
  n <- nrow(X)
  
  ## 3.1 CATE learner: causal forest
  cf <- causal_forest(X, Y, A)
  tau_hat <- predict(cf)$predictions
  
  ## 3.2 Cross-fitted nuisance models and DR pseudo-outcomes
  fold <- sample(rep(1:K, length.out = n))
  mu1  <- mu0 <- ehat <- rep(NA_real_, n)
  
  for (k in 1:K) {
    idx_tr <- which(fold != k)
    idx_te <- which(fold == k)
    
    df_tr  <- data.frame(Y = Y[idx_tr],
                         A = A[idx_tr],
                         X[idx_tr, , drop = FALSE])
    
    # Propensity model
    mod_e <- glm(A ~ X1 + X2 + X3, family = binomial(), data = df_tr)
    ehat[idx_te] <- predict(mod_e,
                            newdata = X[idx_te, , drop = FALSE],
                            type = "response")
    
    # Outcome models
    df_tr_treated  <- df_tr[df_tr$A == 1, ]
    df_tr_control  <- df_tr[df_tr$A == 0, ]
    
    m1_fit <- glm(Y ~ X1 + X2 + X3, family = binomial(), data = df_tr_treated)
    m0_fit <- glm(Y ~ X1 + X2 + X3, family = binomial(), data = df_tr_control)
    
    mu1[idx_te] <- predict(m1_fit,
                           newdata = X[idx_te, , drop = FALSE],
                           type = "response")
    mu0[idx_te] <- predict(m0_fit,
                           newdata = X[idx_te, , drop = FALSE],
                           type = "response")
  }
  
  # Clip propensity scores
  ehat <- pmin(pmax(ehat, 1e-3), 1 - 1e-3)
  
  # DR pseudo-outcome
  W_DR <- mu1 - mu0 +
    A * (Y - mu1) / ehat -
    (1 - A) * (Y - mu0) / (1 - ehat)
  
  ## 3.3 Ranking performance: uplift curve and AUQC
  AUQC <- auqc_from_DR(tau_hat, W_DR)
  
  # (Simulation-only check) ROC/PR vs true Z*
  if (length(unique(dat$Z_true)) == 2L) {
    # ROC
    roc_true   <- roc(dat$Z_true, tau_hat)
    auroc_true <- as.numeric(auc(roc_true))
    
    # PR curve: scores.class0 = positives (Z_true==1) is PRROC's convention
    pr_true <- pr.curve(
      scores.class0 = tau_hat[dat$Z_true == 1],
      scores.class1 = tau_hat[dat$Z_true == 0],
      curve = FALSE
    )
    auprc_true <- pr_true$auc.integral
  } else {
    # Degenerate case: no true benefiters or everyone is a benefiter
    auroc_true <- NA_real_
    auprc_true <- NA_real_
  }
  
  ## 3.4 Policy evaluation: value curve and t*
  t_grid   <- seq(min(tau_hat), max(tau_hat), length.out = 200)
  val_grid <- dr_value_vec(t_grid, tau_hat, A, Y, mu1, mu0, ehat)
  
  best_idx     <- which.max(val_grid)
  t_star_value <- t_grid[best_idx]
  value_max    <- val_grid[best_idx]
  
  ## DR value of the personalized rule at t*  (this should equal value_max up to noise)
  a_star  <- as.integer(tau_hat > t_star_value)
  val_star <- dr_value_policy(a_star, A, Y, mu1, mu0, ehat)
  
  ## DR value of constant treat-all and treat-none policies
  a_all1  <- rep(1L, length(Y))
  a_all0  <- rep(0L, length(Y))
  
  val_all1 <- dr_value_policy(a_all1, A, Y, mu1, mu0, ehat)
  val_all0 <- dr_value_policy(a_all0, A, Y, mu1, mu0, ehat)
  
  ## Best constant policy value
  val_const_best <- max(val_all1, val_all0)
  
  ## Incremental value of personalization vs best constant rule
  value_gain_all <- val_star - val_const_best
  
  # Treat-all and treat-none values for reference
  val_all  <- mean(Y)                                  # everyone treated or not? (randomized trial)
  val_none <- mean(Y[A == 0])                          # crude; optional benchmark
  
  ## 3.5 NP safety dial
  np_out <- np_choose(
    score      = tau_hat,
    W_DR       = W_DR,
    delta      = dat$delta,
    alpha_harm = alpha_harm,
    conf       = 0.95,
    t_grid     = t_grid,
    min_treat  = max(100, 0.05 * n),
    ci_method  = ci_method
  )
  
  list(
    tau_hat      = tau_hat,
    W_DR         = W_DR,
    ehat         = ehat,
    mu1          = mu1,
    mu0          = mu0,
    AUQC         = AUQC,
    t_grid       = t_grid,
    val_grid     = val_grid,
    t_star_value = t_star_value,
    value_max    = value_max,
    value_all    = val_all,
    value_none   = val_none,
    np_out       = np_out,
    auroc_true   = auroc_true,
    auprc_true   = auprc_true,
    value_gain_all = value_gain_all
  )
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
    geom_line(size = 0.6, color = "black") +
    labs(
      x = expression(X[1] ~ "(ordered windows)"),
      y = "Windowed risk difference\n(treated – control)",
      title = paste0(main_label, "A. STEPP")
    ) +
    coord_cartesian(xlim = c(-2, 2)) +
    theme_pub()
  
  ## 4.4 Panel B: Uplift curve
  ord <- order(st2$tau_hat, decreasing = TRUE)
  cum_gain <- cumsum(st2$W_DR[ord])
  frac     <- seq_along(cum_gain) / length(cum_gain)
  df_uplift <- data.frame(frac = frac, cum_gain = cum_gain)
  
  pB <- ggplot(df_uplift, aes(frac, cum_gain)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey40") +
    geom_line(size = 0.6, color = "black") +
    annotate(
      "label",
      x = 0.98, y = max(cum_gain, na.rm = TRUE) * 0.98,
      label = sprintf("AUQC = %.2f", st2$AUQC),
      hjust = 1, vjust = 1, size = 4
    ) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    labs(
      x = "Top fraction by CATE",
      y = "Cumulative DR uplift",
      title = paste0(main_label, "B. Uplift")
    ) +
    theme_pub()
  
  ## 4.5 Panel C: Policy value vs threshold
  df_value <- data.frame(
    threshold = st2$t_grid,
    value     = st2$val_grid
  )
  
  pC <- ggplot(df_value, aes(threshold, value)) +
    geom_line(size = 0.6, colour = "black") +
    geom_vline(xintercept = st2$t_star_value,
               linetype = 2, colour = "grey30") +
    annotate(
      "label",
      x = max(df_value$threshold, na.rm = TRUE),
      y = max(df_value$value,    na.rm = TRUE),
      hjust = 1, vjust = 1,
      label = sprintf("t* = %.3f", st2$t_star_value),
      size = 4
    ) +
    labs(
      x = "Threshold on CATE score",
      y = "Estimated policy value (DR)",
      title = paste0(main_label, "C. Policy value")
    ) +
    theme_pub()
  
  ## 4.6 Panel D: NP-ROC curve (harm vs benefit-capture)
  df_np <- na.omit(st2$np_out$curve)
  t_np  <- st2$np_out$t_star
  pt_np <- df_np[which.min(abs(df_np$t - t_np)), , drop = FALSE]
  
  pD <- ggplot(df_np, aes(harm_upper, det)) +
    geom_rect(aes(xmin = -Inf, xmax = alpha_harm,
                  ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = "grey95", color = NA) +
    geom_vline(xintercept = alpha_harm,
               linetype = 2, color = "grey30") +
    geom_path(size = 0.6, color = "black") +
    geom_point(data = pt_np, colour = "red", size = 2) +
    annotate(
      "text",
      x = pt_np$harm_upper, y = pt_np$det,
      label = "Best attainable",
      vjust = -1, hjust = 1, size = 4, colour = "red"
    ) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(
      x = "Harm rate",
      y = "Benefit-capture rate",
      title = paste0(main_label, "D. NP safety dial")
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

## Example (this is just demonstration; you can run it yourself)
mc_no   <- run_mc_scenario(sim_no_hte,   scenario_label = "No HTE",   B = 200)
mc_weak <- run_mc_scenario(sim_weak_hte, scenario_label = "Weak HTE", B = 200)
mc_str  <- run_mc_scenario(sim_strong_hte, scenario_label = "Strong HTE", B = 200)

mc_all <- bind_rows(mc_no, mc_weak, mc_str)

# Gate proceed rate by scenario
mc_all %>%
  group_by(scenario) %>%
  summarise(
    proceed_rate   = mean(proceed),
    mean_AUQC      = mean(AUQC),
    mean_valuegain = mean(value_gain_all),
    NP_feasible    = mean(NP_feasible)
  )

write.csv(mc_all, 'result/mc_all.csv', row.names = F)

plot_single_trial(sim_strong_hte, main_label = "Strong HTE: ")
plot_single_trial(sim_no_hte, main_label = "No HTE: ")
plot_single_trial(sim_weak_hte, main_label = "Weak HTE: ")

