# library(speff2trial)
# data(ACTG175)
# str(ACTG175)
# 
# write.csv(ACTG175, 'data/ACTG175.csv', row.names = F)

set.seed(175)

need <- c("dplyr", "ggplot2", "data.table", "mgcv", "grf")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(dplyr)
library(ggplot2)
library(data.table)
library(mgcv)   # for optional GAM figure
library(grf)    # for causal forest (Stage 2)

theme_min <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
#df_raw <- fread("data/ACTG175.csv", na.strings = c("", "NA", ".", "NaN")) |> as.data.frame()

## Harmonize common column-name variants (CRAN vs UCI exports)
rename_if_present <- function(dat, from, to){
  for(i in seq_along(from)){
    if (from[i] %in% names(dat) && !(to[i] %in% names(dat))) {
      names(dat)[names(dat) == from[i]] <- to[i]
    }
  }
  dat
}

df_raw <- rename_if_present(
  df_raw,
  from = c("time","cid","trt","cd4.baseline","cd8.baseline"),
  to   = c("days","cens","arms","cd40","cd80")
)

## If a convenience binary treatment 'treat' is not present, create it:
##  treat = 0 for AZT alone (arms==0), 1 for "other" regimens (arms>0)
if (!("treat" %in% names(df_raw)) && ("arms" %in% names(df_raw))) {
  df_raw <- df_raw |> mutate(treat = ifelse(arms == 0, 0L, 1L))
}

## Sanity checks we need for the tutorial
stopifnot(all(c("days","cens") %in% names(df_raw)))
stopifnot("treat" %in% names(df_raw))

## ---- 1-define-data -------------------------------------------------------

## Analysis horizon: 96 weeks (days)
HORIZON_DAYS <- 96 * 7

## Binary event by horizon: 1 if event observed by horizon
event_by_h <- with(df_raw, ifelse(!is.na(days) & !is.na(cens) & days <= HORIZON_DAYS & cens == 1, 1L, 0L))
table(event_by_h)

## Define “good” outcome Y = 1 - event (higher is better)
df <- df_raw |>
  mutate(
    Y = 1L - event_by_h,
    A = as.integer(treat)  # binary treatment 0/1
  )
table(df$A)
table(df$Y)
table(df$arms)

df <- df[df$arms!=3,]
table(df$A)
table(df$Y)
table(df$arms)

## Candidate baseline covariates (keep only those that actually exist)
baseline_candidates <- c(
  "age","wtkg","hemo","homo","drugs","karnof",
  "oprior","z30","zprior","preanti","race","gender",
  "symptom","cd40","cd80","str2","strat"
)
X_names <- intersect(baseline_candidates, names(df))

## Keep complete cases on A, Y, and chosen X (no post-randomization vars!)
keep <- complete.cases(df[, c("A","Y", X_names)])
dat <- df[keep, c("A","Y", X_names)]

## Coerce binary/categorical appropriately
to_factor <- intersect(c("hemo","homo","drugs","race","gender","oprior","z30","zprior","symptom","str2","strat"), X_names)
for (nm in to_factor) dat[[nm]] <- factor(dat[[nm]])

## Quick glance
summary(dat[, c("A","Y")])
length(X_names); X_names
str(dat)

## ---- 2-1-optionA-LRT -----------------------------------------------------
## Build formulas: main effects and main+interactions (A:each X)
X_names <- setdiff(X_names, "zprior")
f_main <- as.formula(paste("Y ~ A +", paste(X_names, collapse = " + ")))
f_int  <- as.formula(paste("Y ~ A +", paste(X_names, collapse = " + "),
                           "+", paste(sprintf("A:%s", X_names), collapse = " + ")))

m_main <- glm(f_main, data = dat, family = binomial())
m_int  <- glm(f_int,  data = dat, family = binomial())

## Likelihood-ratio test for "no heterogeneity" (all A:Xi = 0)
lrt_stat <- 2 * (logLik(m_int) - logLik(m_main))
df_diff  <- attr(logLik(m_int), "df") - attr(logLik(m_main), "df")
p_lrt    <- pchisq(as.numeric(lrt_stat), df = df_diff, lower.tail = FALSE)

cat(sprintf("\n[Option A] Global heterogeneity LRT: chi^2 = %.2f (df=%d), p = %.4g\n",
            lrt_stat, df_diff, p_lrt))

## ---- 2-2-optionC-Wald-Holm ----------------------------------------------
## Choose clinically plausible prespecified biomarkers (edit if desired)
prespec <- intersect(c("cd40","karnof"), X_names)
if (length(prespec) == 0L) prespec <- X_names[1:min(2, length(X_names))]  # fallback

coefs <- summary(m_int)$coefficients
p_int <- sapply(prespec, function(x){
  nm <- paste0("A:", x)
  if (nm %in% rownames(coefs)) coefs[nm, "Pr(>|z|)"] else NA_real_
})

p_int_adj <- p.adjust(p_int, method = "holm")
cat("\n[Option C] Prespecified A×X tests with Holm adjustment:\n")
print(data.frame(term = prespec, p_raw = signif(p_int,3), p_holm = signif(p_int_adj,3)))


## ---- 2-3-stepp-like-plot -------------------------------------------------
## Simple sliding-window risk-difference plot along baseline CD4
dfp <- dat |> select(A, Y, cd40) |> arrange(cd40)

## Choose window parameters (n_win ~ 20 overlapping windows)
n <- nrow(dfp); win <- max(50, floor(n * 0.07)); step <- floor(win / 3)
idx_starts <- seq(1, n - win + 1, by = step)

stepp <- lapply(idx_starts, function(s){
  e <- s + win - 1
  sub <- dfp[s:e, ]
  if (length(unique(sub$A)) < 2) return(NULL)
  p1 <- mean(sub$Y[sub$A == 1], na.rm = TRUE)
  p0 <- mean(sub$Y[sub$A == 0], na.rm = TRUE)
  rd <- p1 - p0
  ## Wald SE for difference in proportions
  n1 <- sum(sub$A == 1); n0 <- sum(sub$A == 0)
  se <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
  data.frame(
    cd4_mid = median(sub$cd40, na.rm = TRUE),
    rd = rd, rd_lo = rd - 1.96*se, rd_hi = rd + 1.96*se,
    n_win = n1 + n0
  )
}) |> bind_rows()

ggplot(stepp, aes(cd4_mid, rd)) +
  geom_hline(yintercept = 0, linetype = 2) + 
  ylim(c(-0.3,0.3)) +
  #geom_ribbon(aes(ymin = rd_lo, ymax = rd_hi), alpha = 0.2) +
  geom_line(size = 0.8) +
  labs(x = "Baseline CD4 (cells/mm³)",
       y = "Risk difference (treated - control)") +
  theme_min

dfk <- dat |>
  dplyr::select(A, Y, karnof) |>
  dplyr::filter(!is.na(karnof)) |>
  dplyr::arrange(karnof)

stepp_k <- dfk |>
  dplyr::group_by(karnof) |>
  dplyr::summarise(
    n1 = sum(A == 1),
    n0 = sum(A == 0),
    p1 = mean(Y[A == 1], na.rm = TRUE),
    p0 = mean(Y[A == 0], na.rm = TRUE),
    rd = p1 - p0,
    se = sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0),
    rd_lo = rd - 1.96 * se,
    rd_hi = rd + 1.96 * se,
    n_win = n1 + n0,
    .groups = "drop"
  )

ggplot(stepp_k, aes(x = karnof, y = rd)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = rd_lo, ymax = rd_hi), width = 2) +
  geom_point(size = 2) +
  geom_line(size = 0.8) +
  labs(x = "Baseline Karnofsky score",
       y = "Risk difference (treated - control)",
       title = NULL) +
  theme_min

## ---- 3-1-crossfit-helpers ------------------------------------------------
K <- 5
n <- nrow(dat)
fold_id <- sample(rep(1:K, length.out = n))

## Build model matrices for convenience
X <- dat[, X_names, drop = FALSE]
Y <- dat$Y
A <- dat$A

## Helper: fit propensity e(X) and outcome models m_a(X) by fold
fit_nuisance_by_fold <- function(fold){
  idx_tr <- which(fold_id != fold)
  idx_te <- which(fold_id == fold)
  dtr <- dat[idx_tr, , drop = FALSE]
  dte <- dat[idx_te, , drop = FALSE]
  
  ## Propensity (logistic)
  f_e <- as.formula(paste("A ~", paste(X_names, collapse = " + ")))
  e_fit <- glm(f_e, data = dtr, family = binomial())
  e_hat <- plogis(predict(e_fit, newdata = dte))
  
  ## Outcome models (two logistic GLMs by arm)
  f_y <- as.formula(paste("Y ~", paste(X_names, collapse = " + ")))
  m0 <- glm(f_y, data = dtr[dtr$A == 0, ], family = binomial())
  m1 <- glm(f_y, data = dtr[dtr$A == 1, ], family = binomial())
  m0_hat <- plogis(predict(m0, newdata = dte))
  m1_hat <- plogis(predict(m1, newdata = dte))
  
  ## DR pseudo-outcome for ranking (CATE-type)
  ## tau_dr = (m1 - m0) + A/e*(Y - m1) - (1-A)/(1-e)*(Y - m0)
  Y_te <- dte$Y; A_te <- dte$A
  tau_dr <- (m1_hat - m0_hat) +
    (A_te / pmax(e_hat, 1e-4)) * (Y_te - m1_hat) -
    ((1 - A_te) / pmax(1 - e_hat, 1e-4)) * (Y_te - m0_hat)
  
  data.frame(
    idx = idx_te,
    e_hat = e_hat,
    m0_hat = m0_hat,
    m1_hat = m1_hat,
    tau_dr = tau_dr
  )
}

## Compute cross-fitted nuisance predictions & DR pseudo-outcomes
nus <- lapply(1:K, fit_nuisance_by_fold) |> bind_rows() |> arrange(idx)

## Cross-fitted causal forest CATE (trained out-of-fold)
tau_hat <- rep(NA_real_, n)

# One-hot encode factors, drop intercept
X.num <- model.matrix(~ . - 1, data = X)

for (fold in 1:K) {
  idx_tr <- which(fold_id != fold)
  idx_te <- which(fold_id == fold)
  cf <- causal_forest(X.num[idx_tr, , drop = FALSE], Y[idx_tr], A[idx_tr],
                      num.trees = 2000, seed = 175)
  tau_hat[idx_te] <- predict(cf, X.num[idx_te, , drop = FALSE])$predictions
}

## Gather evaluation frame
eval_df <- dat |>
  mutate(
    tau_hat = tau_hat,
    e_hat   = nus$e_hat,
    m0_hat  = nus$m0_hat,
    m1_hat  = nus$m1_hat,
    tau_dr  = nus$tau_dr
  )
summary(eval_df$tau_hat)


## ---- 3-2-uplift-AUQC -----------------------------------------------------
## Order by predicted benefit (descending) and accumulate DR pseudo-outcomes
ord <- order(eval_df$tau_hat, decreasing = TRUE, na.last = NA)
cum_dr <- cumsum(eval_df$tau_dr[ord]) / nrow(eval_df)

## Points at budgets q in {0.05, 0.10, ..., 1.00}
q_grid <- seq(0.05, 1, by = 0.05)
nq     <- floor(q_grid * nrow(eval_df))
U_q    <- cum_dr[nq]

uplift_df <- data.frame(q = q_grid, U = U_q)

## AUQC as Riemann sum (simple discrete integral over q)
AUQC <- sum(U_q) * 0.05

cat(sprintf("\n[Stage 2] AUQC (units of outcome gain) = %.3f\n", AUQC))

ggplot(uplift_df, aes(q, U)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  annotate("label",   x = .7,
    y = max(uplift_df$U, na.rm = TRUE),
    hjust = 1, vjust = 1,
    label = sprintf("AUQC = %.3f", AUQC)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = function(x) paste0(100 * x, "%")) +
  labs(x = "Top fraction by CATE",
    y = "Cumulative DR uplift") +
  theme_min +
  theme(plot.title = element_blank(),
    panel.grid.minor = element_blank())


## ---- 3-3-policy-value ----------------------------------------------------
## DR value for a candidate policy pi_t(x) = 1{ tau_hat >= t }
dr_value <- function(t){
  a_hat <- as.integer(eval_df$tau_hat >= t)
  ## m_hat(X, a_hat)
  m_hat <- ifelse(a_hat == 1, eval_df$m1_hat, eval_df$m0_hat)
  ## DR term
  w <- ifelse(eval_df$A == 1, 1 / pmax(eval_df$e_hat, 1e-4),
              1 / pmax(1 - eval_df$e_hat, 1e-4))
  aug <- ifelse(eval_df$A == a_hat, w * (eval_df$Y - ifelse(eval_df$A == 1, eval_df$m1_hat, eval_df$m0_hat)), 0)
  mean(m_hat + aug)
}

t_grid <- quantile(na.omit(eval_df$tau_hat), probs = seq(0.05, 0.95, by = 0.05))
val_grid <- sapply(t_grid, dr_value)
best_idx <- which.max(val_grid)
t_best   <- t_grid[best_idx]; v_best <- val_grid[best_idx]

cat(sprintf("[Stage 2] Best policy threshold t = %.4f; DR value = %.4f\n", t_best, v_best))

pv_df <- data.frame(threshold = as.numeric(t_grid), value = as.numeric(val_grid))

ggplot(pv_df, aes(threshold, value)) +
  geom_line(size = 1) + geom_point(size = 2) +
  geom_vline(xintercept = t_best, linetype = 2) +
  annotate("text", x = t_best, y = max(val_grid), vjust = -0.5,
           label = sprintf("best t = %.3f", t_best)) +
  labs(x = "Threshold on predicted benefit (tau_hat)",
       y = "Estimated policy value (DR)",
       title = "Policy value vs. threshold (cross-fitted, DR)") +
  theme_min

v_max <- max(pv_df$value, na.rm = TRUE)
ggplot(pv_df, aes(threshold, value)) +
  # smooth black line, no big points (same as simulation)
  geom_line(size = 0.6, colour = "black") +
  # dashed vertical at the chosen threshold t*
  geom_vline(xintercept = t_best, linetype = 2, colour = "grey30") +
  # label box in the upper-right corner, like AUQC label in sim plots
  annotate(
    "label",
    x = max(pv_df$threshold, na.rm = TRUE),
    y = v_max,
    hjust = 1, vjust = 1,
    label = sprintf("t* = %.3f", t_best),
    size = 5
  ) +
  labs(
    title = NULL,
    x = "Threshold on CATE",
    y = "Estimated policy value"
  ) +
  theme_min +
  theme(
    plot.title       = element_blank(),
    panel.grid.minor = element_blank()
  )


## ---- 3-4-NP-rule ---------------------------------------------------------
## Clinical margin delta for "benefit" on the DR scale (use 0 for "any benefit")
delta <- 0.0

## Surrogates: harm if DR pseudo-outcome < delta; benefit-capture if >= delta
eval_df <- eval_df |>
  mutate(harm_sur = as.integer(tau_dr <  delta),
         bene_sur = as.integer(tau_dr >= delta))

## For a policy pi_t, compute:
## - harm rate among treated: P(harm | pi_t(X)=1)
## - benefit-capture among all: P(bene & treated)
np_metrics <- function(t){
  a_hat <- as.integer(eval_df$tau_hat >= t)
  treated <- which(a_hat == 1)
  harm_rate <- if (length(treated) == 0) NA_real_ else mean(eval_df$harm_sur[treated])
  #bene_cap  <- mean(eval_df$bene_sur[a_hat == 1])
  bene_cap <- mean(eval_df$bene_sur * a_hat)
  treat_rate <- mean(a_hat)
  ppv_benefit <- if (length(treated) == 0) NA_real_ else mean(eval_df$bene_sur[treated]) 
  c(harm = harm_rate, bene = bene_cap, treat = treat_rate, ppv = ppv_benefit)
}

t_grid2 <- quantile(na.omit(eval_df$tau_hat), probs = seq(0.05, 0.95, by = 0.05))
np_mat <- t(sapply(t_grid2, np_metrics))
np_df <- data.frame(threshold = as.numeric(t_grid2),
                    harm = np_mat[, "harm"],
                    bene = np_mat[, "bene"],
                    treat = np_mat[, "treat"],
                    ppv = np_mat[, "ppv"])

## Choose smallest t whose harm <= alpha_harm; among ties pick largest bene
alpha_harm <- 0.10
ok <- which(!is.na(np_df$harm) & np_df$harm <= alpha_harm)
if (length(ok)){
  idx_np <- ok[ which.max(np_df$bene[ok]) ]
  t_np <- np_df$threshold[idx_np]; harm_np <- np_df$harm[idx_np]; bene_np <- np_df$bene[idx_np]
  cat(sprintf("[Stage 2] NP rule (alpha_harm=%.2f): threshold t = %.4f; harm = %.3f; benefit-capture = %.3f\n",
              alpha_harm, t_np, harm_np, bene_np))
} else {
  t_np <- NA_real_
  cat(sprintf("[Stage 2] NP rule: no threshold achieves harm <= %.2f; reporting frontier instead.\n", alpha_harm))
}

library(scales)  # for percent_format
## choose the "best attainable" point on the frontier:
## here: the point with the largest benefit capture
idx_best  <- which.min(np_df$harm)
harm_best <- np_df$harm[idx_best]
bene_best <- np_df$bene[idx_best]
best_df   <- data.frame(harm = harm_best, bene = bene_best)

## pretty percent labels
pct_lab <- function(x) paste0(round(100 * x), "%")

ggplot(np_df, aes(harm, bene)) +
  # NP constraint line at alpha_harm
  geom_vline(xintercept = alpha_harm, linetype = 2) +
  # frontier
  geom_path(size = 0.7) +
  geom_point(size = 1) +
  # highlight best-attainable point
  geom_point(data = best_df,
             aes(harm, bene),
             colour = "red", size = 2) +
  annotate("text",
           x = 0.21, y = bene_best,
           label = "Best attainable",
           colour = "red",
           hjust = -0.05, vjust = 1) +
  # axes in percent, like the simulation figure
  scale_x_continuous(
    name   = "Harm rate",
    labels = pct_lab
  ) +
  scale_y_continuous(
    name   = "Benefit-capture rate",
    labels = pct_lab
  ) +
  theme_min +
  theme(
    plot.title       = element_blank(),
    panel.grid.minor = element_blank()
  )
