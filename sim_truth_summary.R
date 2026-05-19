expit <- function(z) 1 / (1 + exp(-z))

baseline_risk <- function(X) {
  0.25 + 0.35 * expit(-0.4 + 0.45 * X$X1 - 0.25 * X$X2 + 0.35 * X$X3)
}

scenario_tau <- function(X, scenario) {
  switch(
    scenario,
    "No treatment effect" = rep(0, nrow(X)),
    "Constant benefit without HTE" = rep(0.04, nrow(X)),
    "Weak quantitative HTE" = 0.04 + 0.015 * tanh(X$X1 / 1.2),
    "Strong quantitative HTE" = 0.06 + 0.050 * tanh(X$X1),
    "Weak qualitative HTE" = 0.080 * tanh(X$X1 / 1.1),
    "Strong qualitative HTE" = 0.180 * tanh(X$X1),
    stop("Unknown simulation scenario: ", scenario)
  )
}

summarize_truth <- function(X, scenario, delta = 0.03) {
  p0 <- baseline_risk(X)
  tau <- scenario_tau(X, scenario)
  p1 <- p0 + tau
  treat_oracle <- tau > 0
  treat_margin <- tau > delta
  value_none <- mean(p0)
  value_all <- mean(p1)
  fixed_best <- max(value_none, value_all)
  value_oracle <- mean(ifelse(treat_oracle, p1, p0))
  value_margin <- mean(ifelse(treat_margin, p1, p0))
  data.frame(
    scenario = scenario,
    ate = mean(tau),
    tau_sd = sd(tau),
    tau_q05 = unname(quantile(tau, 0.05)),
    tau_q95 = unname(quantile(tau, 0.95)),
    benefit_rate_delta = mean(tau > delta),
    harm_rate = mean(tau < 0),
    value_treat_none = value_none,
    value_treat_all = value_all,
    value_oracle = value_oracle,
    oracle_value_gain = value_oracle - fixed_best,
    value_margin_policy = value_margin,
    margin_policy_value_gain = value_margin - fixed_best,
    oracle_treat_fraction = mean(treat_oracle),
    margin_policy_treat_fraction = mean(treat_margin)
  )
}

set.seed(20260517)
n <- as.integer(Sys.getenv("SIM_TRUTH_N", "1000000"))
X <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n),
  X3 = rbinom(n, 1, 0.5)
)

scenarios <- c(
  "No treatment effect",
  "Constant benefit without HTE",
  "Weak quantitative HTE",
  "Strong quantitative HTE",
  "Weak qualitative HTE",
  "Strong qualitative HTE"
)

truth <- do.call(rbind, lapply(scenarios, function(s) summarize_truth(X, s)))
truth$n_truth <- n

dir.create("result", showWarnings = FALSE)
write.csv(truth, "result/sim_truth_summary_step6.csv", row.names = FALSE)
print(truth)
