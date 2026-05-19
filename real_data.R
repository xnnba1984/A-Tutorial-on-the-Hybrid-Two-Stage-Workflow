## ACTG 175 real-data analysis for the two-stage HTE-to-policy workflow.
## This version uses IPCW for the 96-week event-free endpoint so that
## early censoring is handled explicitly rather than counted as event-free.

set.seed(175)

required_packages <- c(
  "data.table", "dplyr", "ggplot2", "grf", "survival",
  "sandwich", "lmtest", "car", "patchwork", "scales", "gridExtra"
)
to_install <- setdiff(required_packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grf)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(car)
  library(patchwork)
  library(scales)
  library(gridExtra)
})

if (!file.exists("data/ACTG175.csv")) {
  stop("Run real_data.R from the project root containing data/ACTG175.csv.")
}

dir.create("result", showWarnings = FALSE)
dir.create("result/figures", recursive = TRUE, showWarnings = FALSE)

## ---- 0-helpers ------------------------------------------------------------

theme_pub <- function(base_size = 10.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(7, 10, 7, 14)
    )
}

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

pct_lab <- function(x) paste0(round(100 * x), "%")

extract_chisq_test <- function(lh) {
  out <- as.data.frame(lh)
  data.frame(
    chisq = out$Chisq[2],
    df = out$Df[2],
    p_value = out$`Pr(>Chisq)`[2]
  )
}

## ---- 1-data-and-ipcw ------------------------------------------------------

HORIZON_DAYS <- 96 * 7
ALPHA_HARM <- 0.10
BENEFIT_MARGIN <- 0.0

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

stopifnot(all(c("days", "cens", "arms", "treat") %in% names(df_raw)))

df <- df_raw |>
  mutate(
    A = as.integer(treat),
    event_by_h = as.integer(!is.na(days) & !is.na(cens) &
                              days <= HORIZON_DAYS & cens == 1L),
    early_censored = as.integer(!is.na(days) & !is.na(cens) &
                                  days <= HORIZON_DAYS & cens == 0L),
    R_h = as.integer(event_by_h == 1L | days > HORIZON_DAYS),
    Y = ifelse(event_by_h == 1L, 0,
               ifelse(days > HORIZON_DAYS, 1, NA_real_)),
    censor_time = pmin(days, HORIZON_DAYS),
    censor_event = early_censored
  ) |>
  filter(arms != 3)

baseline_candidates <- c(
  "age", "wtkg", "hemo", "homo", "drugs", "karnof",
  "oprior", "z30", "zprior", "preanti", "race", "gender",
  "symptom", "cd40", "cd80", "str2", "strat"
)

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

## Parsimonious censoring model. This avoids unstable separation from rare
## censoring events while conditioning on treatment and major baseline severity.
censor_names <- intersect(
  c("A", "age", "wtkg", "karnof", "preanti", "cd40", "cd80", "strat"),
  names(dat_all)
)
censor_rhs <- paste(censor_names, collapse = " + ")

fit_ipcw <- function(train, eval) {
  f_cens <- as.formula(paste(
    "Surv(censor_time, censor_event) ~", censor_rhs
  ))
  fit <- survival::coxph(f_cens, data = train, ties = "efron", x = TRUE)
  bh <- survival::basehaz(fit, centered = FALSE)
  eval_time <- pmin(eval$days, HORIZON_DAYS)
  H0 <- approx(
    bh$time, bh$hazard, xout = eval_time,
    method = "constant", f = 0, rule = 2, yleft = 0
  )$y
  lp <- predict(fit, newdata = eval, type = "lp", reference = "zero")
  pmax(exp(-H0 * exp(lp)), 0.02)
}

dat_all$G_hat <- fit_ipcw(dat_all, dat_all)
dat_all$ipcw <- ifelse(dat_all$R_h == 1L, 1 / dat_all$G_hat, 0)

cat(sprintf(
  "\n[Data] n=%d, observed horizon outcomes=%d, early censored before 96 weeks=%d\n",
  nrow(dat_all), sum(dat_all$R_h), sum(dat_all$early_censored)
))
cat(sprintf(
  "[Data] events by 96 weeks=%d, event-free through 96 weeks=%d\n",
  sum(dat_all$event_by_h), sum(dat_all$R_h == 1L & dat_all$Y == 1)
))
cat(sprintf(
  "[IPCW] G range %.3f to %.3f; observed weights %.3f to %.3f\n",
  min(dat_all$G_hat), max(dat_all$G_hat),
  min(dat_all$ipcw[dat_all$R_h == 1L]),
  max(dat_all$ipcw[dat_all$R_h == 1L])
))

## ---- 2-stage-1 ------------------------------------------------------------

rhs <- paste(X_names, collapse = " + ")
stage1_dat <- dat_all |>
  filter(R_h == 1L)

f_main <- as.formula(paste("Y ~ A +", rhs))
f_int <- as.formula(paste("Y ~ A * (", rhs, ")"))

lm_main <- lm(f_main, data = stage1_dat, weights = ipcw)
lm_int <- lm(f_int, data = stage1_dat, weights = ipcw)
vc_hc3 <- sandwich::vcovHC(lm_int, type = "HC3")

interaction_terms <- grep("^A:", names(coef(lm_int)), value = TRUE)
interaction_terms <- interaction_terms[!is.na(coef(lm_int)[interaction_terms])]

gate_test <- car::linearHypothesis(
  lm_int,
  paste0(interaction_terms, " = 0"),
  vcov. = vc_hc3,
  test = "Chisq",
  singular.ok = TRUE
)
gate_summary <- extract_chisq_test(gate_test)

coef_test <- lmtest::coeftest(lm_int, vcov. = vc_hc3)
prespec <- intersect(c("cd40", "karnof"), X_names)
option_c <- lapply(prespec, function(x) {
  term <- paste0("A:", x)
  if (!(term %in% rownames(coef_test))) {
    return(data.frame(term = x, coefficient = NA_real_, se = NA_real_,
                      p_raw = NA_real_))
  }
  data.frame(
    term = x,
    coefficient = coef_test[term, "Estimate"],
    se = coef_test[term, "Std. Error"],
    p_raw = coef_test[term, "Pr(>|t|)"]
  )
}) |>
  bind_rows() |>
  mutate(p_holm = p.adjust(p_raw, method = "holm"))

cat(sprintf(
  "\n[Stage 1] IPCW robust Wald gate: chi^2 = %.2f (df=%d), p = %.4g\n",
  gate_summary$chisq, gate_summary$df, gate_summary$p_value
))
cat("[Stage 1] Prespecified Option C terms with Holm adjustment:\n")
print(option_c)

## ---- 2-figures-stage-1 ----------------------------------------------------

stepp_dat <- stage1_dat |>
  arrange(cd40)
n_stepp <- nrow(stepp_dat)
win <- max(50, floor(n_stepp * 0.07))
step <- floor(win / 3)
idx_starts <- seq(1, n_stepp - win + 1, by = step)

stepp_cd4 <- lapply(idx_starts, function(s) {
  e <- s + win - 1
  sub <- stepp_dat[s:e, ]
  if (length(unique(sub$A)) < 2) return(NULL)
  p1 <- weighted_mean(sub$Y[sub$A == 1], sub$ipcw[sub$A == 1])
  p0 <- weighted_mean(sub$Y[sub$A == 0], sub$ipcw[sub$A == 0])
  data.frame(
    cd4_mid = median(sub$cd40, na.rm = TRUE),
    rd = p1 - p0,
    n_win = nrow(sub)
  )
}) |>
  bind_rows()

stepp_karnof <- stage1_dat |>
  filter(!is.na(karnof)) |>
  group_by(karnof) |>
  summarise(
    n1 = sum(A == 1),
    n0 = sum(A == 0),
    p1 = weighted_mean(Y[A == 1], ipcw[A == 1]),
    p0 = weighted_mean(Y[A == 0], ipcw[A == 0]),
    rd = p1 - p0,
    n_group = n(),
    .groups = "drop"
  )

p_cd4 <- ggplot(stepp_cd4, aes(cd4_mid, rd)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 0.7) +
  coord_cartesian(ylim = c(-0.30, 0.30)) +
  labs(
    x = "Baseline CD4 (cells/mm^3)",
    y = "IPCW risk difference\n(treated - control)"
  ) +
  theme_pub()

p_karnof <- ggplot(stepp_karnof, aes(karnof, rd)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  labs(
    x = "Baseline Karnofsky score",
    y = "IPCW risk difference\n(treated - control)"
  ) +
  theme_pub()

## ---- 3-stage-2-cross-fitting ---------------------------------------------

K <- 5
n <- nrow(dat_all)
fold_id <- sample(rep(1:K, length.out = n))
X_num <- model.matrix(as.formula(paste("~", rhs, "- 1")), data = dat_all)

tau_hat <- e_hat <- m0_hat <- m1_hat <- tau_dr <- G_cf <- rep(NA_real_, n)

for (fold in 1:K) {
  idx_tr <- which(fold_id != fold)
  idx_te <- which(fold_id == fold)
  dtr <- dat_all[idx_tr, , drop = FALSE]
  dte <- dat_all[idx_te, , drop = FALSE]
  
  G_tr <- fit_ipcw(dtr, dtr)
  G_te <- fit_ipcw(dtr, dte)
  ipcw_tr <- ifelse(dtr$R_h == 1L, 1 / G_tr, 0)
  
  e_fit <- glm(
    as.formula(paste("A ~", rhs)),
    data = dtr,
    family = binomial()
  )
  e_hat[idx_te] <- plogis(predict(e_fit, newdata = dte))
  
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
  m0_hat[idx_te] <- plogis(predict(m0, newdata = dte))
  m1_hat[idx_te] <- plogis(predict(m1, newdata = dte))
  
  cf_obs_idx <- idx_tr[obs_tr]
  cf <- causal_forest(
    X_num[cf_obs_idx, , drop = FALSE],
    dat_all$Y[cf_obs_idx],
    dat_all$A[cf_obs_idx],
    sample.weights = ipcw_tr[obs_tr],
    num.trees = 2000,
    seed = 175 + fold
  )
  tau_hat[idx_te] <- predict(cf, X_num[idx_te, , drop = FALSE])$predictions
  G_cf[idx_te] <- G_te
  
  R_te <- dte$R_h
  Y_te <- ifelse(R_te == 1L, dte$Y, 0)
  A_te <- dte$A
  e_te <- pmin(pmax(e_hat[idx_te], 0.02), 0.98)
  m0_te <- m0_hat[idx_te]
  m1_te <- m1_hat[idx_te]
  
  tau_dr[idx_te] <- (m1_te - m0_te) +
    R_te / G_te * (
      (A_te / e_te) * (Y_te - m1_te) -
        ((1 - A_te) / (1 - e_te)) * (Y_te - m0_te)
    )
}

eval_df <- dat_all |>
  mutate(
    fold_id = fold_id,
    tau_hat = tau_hat,
    e_hat = e_hat,
    m0_hat = m0_hat,
    m1_hat = m1_hat,
    tau_dr = tau_dr,
    G_cf = G_cf
  )

## ---- 3-2-uplift -----------------------------------------------------------

ord <- order(eval_df$tau_hat, decreasing = TRUE, na.last = NA)
q_grid <- seq(0.05, 1, by = 0.05)
nq <- pmax(1, floor(q_grid * n))
U_q <- cumsum(eval_df$tau_dr[ord])[nq] / n
U_random <- q_grid * mean(eval_df$tau_dr, na.rm = TRUE)
U_centered <- U_q - U_random
uplift_df <- data.frame(q = q_grid, U = U_q, U_random = U_random, U_centered = U_centered)
AUQC_raw <- sum(U_q) * 0.05
AUQC_centered <- sum(U_centered) * 0.05

p_uplift <- ggplot(uplift_df, aes(q, U_centered)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate(
    "label",
    x = 0.98,
    y = max(uplift_df$U_centered, na.rm = TRUE),
    hjust = 1,
    vjust = 1,
    label = sprintf("cAUQC = %.3f", AUQC_centered),
    size = 3
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = pct_lab
  ) +
  labs(x = "Top fraction by CATE", y = "Centered cumulative DR uplift") +
  theme_pub()

cat(sprintf(
  "\n[Stage 2] centered AUQC = %.3f; raw AUQC = %.3f; U(1) = %.3f\n",
  AUQC_centered, AUQC_raw, tail(U_q, 1)
))

## ---- 3-3-policy-value -----------------------------------------------------

dr_value <- function(t) {
  a_hat <- as.integer(eval_df$tau_hat >= t)
  m_hat <- ifelse(a_hat == 1L, eval_df$m1_hat, eval_df$m0_hat)
  m_A <- ifelse(eval_df$A == 1L, eval_df$m1_hat, eval_df$m0_hat)
  e <- pmin(pmax(eval_df$e_hat, 0.02), 0.98)
  p_A <- ifelse(eval_df$A == 1L, e, 1 - e)
  aug <- ifelse(
    eval_df$R_h == 1L & eval_df$A == a_hat,
    (1 / eval_df$G_cf) / p_A * (eval_df$Y - m_A),
    0
  )
  mean(m_hat + aug)
}

t_grid <- quantile(na.omit(eval_df$tau_hat), probs = seq(0.05, 0.95, by = 0.05))
val_grid <- sapply(t_grid, dr_value)
best_idx <- which.max(val_grid)
t_best <- as.numeric(t_grid[best_idx])
v_best <- as.numeric(val_grid[best_idx])
v_all <- dr_value(-Inf)
v_none <- dr_value(Inf)
v_fixed <- max(v_all, v_none)
fixed_rule <- ifelse(v_all >= v_none, "treat all", "treat none")
value_gain <- v_best - v_fixed

pv_df <- data.frame(
  threshold = as.numeric(t_grid),
  value = as.numeric(val_grid),
  treat_fraction = sapply(t_grid, function(t) mean(eval_df$tau_hat >= t))
)

p_policy <- ggplot(pv_df, aes(threshold, value)) +
  geom_line(linewidth = 0.7, colour = "black") +
  geom_vline(xintercept = t_best, linetype = 2, colour = "grey30") +
  geom_hline(yintercept = v_all, linetype = 3, colour = "grey50") +
  annotate(
    "label",
    x = max(pv_df$threshold, na.rm = TRUE),
    y = max(pv_df$value, na.rm = TRUE),
    hjust = 1,
    vjust = 1,
    label = sprintf("t* = %.3f", t_best),
    size = 3
  ) +
  labs(x = "Threshold on CATE", y = "Policy value") +
  theme_pub()

cat(sprintf(
  "[Stage 2] Best threshold t = %.4f; value = %.4f; %s value = %.4f; gain = %.4f\n",
  t_best, v_best, fixed_rule, v_fixed, value_gain
))

## ---- 3-4-harm-frontier ----------------------------------------------------

eval_df <- eval_df |>
  mutate(
    harm_sur = as.integer(tau_dr < BENEFIT_MARGIN),
    bene_sur = as.integer(tau_dr >= BENEFIT_MARGIN)
  )

np_metrics <- function(t) {
  a_hat <- as.integer(eval_df$tau_hat >= t)
  treated <- a_hat == 1L
  c(
    harm = ifelse(sum(treated) == 0L, NA_real_, mean(eval_df$harm_sur[treated])),
    bene = mean(eval_df$bene_sur * a_hat),
    treat = mean(a_hat),
    ppv = ifelse(sum(treated) == 0L, NA_real_, mean(eval_df$bene_sur[treated]))
  )
}

np_mat <- t(sapply(t_grid, np_metrics))
np_df <- data.frame(
  threshold = as.numeric(t_grid),
  harm = np_mat[, "harm"],
  bene = np_mat[, "bene"],
  treat = np_mat[, "treat"],
  ppv = np_mat[, "ppv"]
)

ok <- which(!is.na(np_df$harm) & np_df$harm <= ALPHA_HARM)
np_feasible <- length(ok) > 0
if (np_feasible) {
  idx_np <- ok[which.max(np_df$bene[ok])]
} else {
  idx_np <- which.min(np_df$harm)
}

best_np <- np_df[idx_np, ]

p_np <- ggplot(np_df, aes(harm, bene)) +
  geom_vline(xintercept = ALPHA_HARM, linetype = 2) +
  geom_path(linewidth = 0.7) +
  geom_point(size = 1) +
  geom_point(data = best_np, aes(harm, bene), colour = "red", size = 2) +
  annotate(
    "text",
    x = best_np$harm - 0.003,
    y = best_np$bene,
    label = ifelse(np_feasible, "Selected", "Best attainable"),
    colour = "red",
    hjust = 1,
    vjust = -0.2
  ) +
  scale_x_continuous(name = "Harm rate", labels = pct_lab) +
  scale_y_continuous(name = "Benefit capture", labels = pct_lab) +
  theme_pub()

if (np_feasible) {
  cat(sprintf(
    "[Stage 2] NP rule feasible: threshold %.4f; harm %.3f; benefit capture %.3f\n",
    best_np$threshold, best_np$harm, best_np$bene
  ))
} else {
  cat(sprintf(
    "[Stage 2] NP rule infeasible at alpha=%.2f; minimum harm %.3f; benefit capture %.3f\n",
    ALPHA_HARM, best_np$harm, best_np$bene
  ))
}

## ---- 4-save-outputs -------------------------------------------------------

tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 12),
  plot.tag.position = c(0, 1)
)

equalize_plot_grobs <- function(plots) {
  grobs <- lapply(plots, ggplotGrob)
  max_widths <- do.call(grid::unit.pmax, lapply(grobs, `[[`, "widths"))
  max_heights <- do.call(grid::unit.pmax, lapply(grobs, `[[`, "heights"))
  lapply(grobs, function(grob) {
    grob$widths <- max_widths
    grob$heights <- max_heights
    grob
  })
}

figure3_panels <- equalize_plot_grobs(list(
  p_cd4 + labs(tag = "A") + tag_theme,
  p_karnof + labs(tag = "B") + tag_theme,
  p_uplift + labs(tag = "C") + tag_theme,
  p_policy + labs(tag = "D") + tag_theme,
  p_np + labs(tag = "E") + tag_theme
))

figure3 <- gridExtra::arrangeGrob(
  grobs = figure3_panels,
  layout_matrix = rbind(
    c(1, 1, 2, 2, 3, 3),
    c(NA, NA, NA, NA, NA, NA),
    c(NA, 4, 4, 5, 5, NA)
  ),
  widths = rep(1, 6),
  heights = c(1, 0.12, 1)
)

ggsave(
  filename = "result/figures/actg175_figure3_ipcw.png",
  plot = figure3,
  width = 7.2,
  height = 5.65,
  dpi = 300,
  bg = "white"
)

summary_df <- data.frame(
  metric = c(
    "horizon_days",
    "raw_n",
    "analysis_n",
    "observed_horizon_outcomes",
    "early_censored_before_horizon",
    "events_by_horizon",
    "event_free_through_horizon",
    "covariate_count",
    "ipcw_min_observed",
    "ipcw_max_observed",
    "stage1_global_chisq",
    "stage1_global_df",
    "stage1_global_p",
    "stage1_proceed_alpha_0_05",
    "centered_auqc",
    "raw_auqc",
    "uplift_at_1",
    "best_threshold",
    "best_threshold_treat_fraction",
    "best_threshold_value",
    "treat_all_value",
    "treat_none_value",
    "fixed_comparator",
    "value_gain_vs_fixed",
    "np_alpha_harm",
    "np_feasible",
    "np_selected_or_best_threshold",
    "np_harm",
    "np_benefit_capture",
    "np_treat_fraction",
    "np_ppv_benefit",
    "tau_hat_min",
    "tau_hat_median",
    "tau_hat_max",
    "tau_dr_mean"
  ),
  value = c(
    HORIZON_DAYS,
    nrow(df_raw),
    nrow(dat_all),
    sum(dat_all$R_h),
    sum(dat_all$early_censored),
    sum(dat_all$event_by_h),
    sum(dat_all$R_h == 1L & dat_all$Y == 1),
    length(X_names),
    min(dat_all$ipcw[dat_all$R_h == 1L]),
    max(dat_all$ipcw[dat_all$R_h == 1L]),
    gate_summary$chisq,
    gate_summary$df,
    gate_summary$p_value,
    as.integer(gate_summary$p_value < 0.05),
    AUQC_centered,
    AUQC_raw,
    tail(U_q, 1),
    t_best,
    pv_df$treat_fraction[best_idx],
    v_best,
    v_all,
    v_none,
    fixed_rule,
    value_gain,
    ALPHA_HARM,
    as.integer(np_feasible),
    best_np$threshold,
    best_np$harm,
    best_np$bene,
    best_np$treat,
    best_np$ppv,
    min(eval_df$tau_hat, na.rm = TRUE),
    median(eval_df$tau_hat, na.rm = TRUE),
    max(eval_df$tau_hat, na.rm = TRUE),
    mean(eval_df$tau_dr, na.rm = TRUE)
  )
)

write.csv(summary_df, "result/actg175_ipcw_summary.csv", row.names = FALSE)
write.csv(option_c, "result/actg175_ipcw_optionC.csv", row.names = FALSE)
write.csv(stepp_cd4, "result/actg175_ipcw_stepp_cd4.csv", row.names = FALSE)
write.csv(stepp_karnof, "result/actg175_ipcw_stepp_karnofsky.csv", row.names = FALSE)
write.csv(uplift_df, "result/actg175_ipcw_uplift_curve.csv", row.names = FALSE)
write.csv(pv_df, "result/actg175_ipcw_policy_curve.csv", row.names = FALSE)
write.csv(np_df, "result/actg175_ipcw_np_curve.csv", row.names = FALSE)
write.csv(eval_df, "result/actg175_ipcw_eval.csv", row.names = FALSE)

capture.output(
  sessionInfo(),
  file = "result/actg175_ipcw_session_info.txt"
)

cat("\n[Saved]\n")
cat("- result/actg175_ipcw_summary.csv\n")
cat("- result/actg175_ipcw_optionC.csv\n")
cat("- result/figures/actg175_figure3_ipcw.png\n")
cat("- result/actg175_ipcw_session_info.txt\n")
