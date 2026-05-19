suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

setwd("~/Library/CloudStorage/Box-Box/Xi/HT_vs_ML")

source_prefix <- readLines("sim_1.R", warn = FALSE)
cut_at <- grep("^## ---- 5-mc-wrapper", source_prefix)
if (length(cut_at) != 1L) {
  stop("Could not locate Monte Carlo wrapper boundary in sim_1.R")
}
eval(parse(text = paste(source_prefix[seq_len(cut_at - 1L)], collapse = "\n")))

out_dir <- "result/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

make_single <- function(sim_fun, label, seed, output) {
  set.seed(seed)
  out <- plot_single_trial(
    sim_fun = sim_fun,
    n = 2000,
    delta = 0.03,
    alpha_harm = 0.10,
    main_label = "",
    plot_base_size = 9.2,
    annotation_size = 2.45
  )

  panels <- equalize_plot_grobs(list(out$pA, out$pB, out$pC, out$pD))
  fig <- gridExtra::arrangeGrob(grobs = panels, ncol = 2)
  ggsave(
    filename = file.path(out_dir, output),
    plot = fig,
    width = 7.2,
    height = 5.4,
    dpi = 450,
    bg = "white"
  )

  data.frame(
    figure = output,
    scenario = label,
    seed = seed,
    n = length(out$dat$Y),
    p_global = out$st1$p_global,
    proceed = out$st1$proceed,
    centered_AUQC = out$st2$AUQC,
    raw_AUQC = out$st2$AUQC_raw,
    value_gain_vs_selected_fixed = out$st2$value_gain_all,
    selected_fixed_treat = out$st2$selected_const_treat,
    learned_treat_fraction = out$st2$learned_treat_fraction,
    truth_value_gain_learned = out$st2$truth_value_gain_learned,
    NP_feasible = out$st2$np_out$feasible,
    NP_harm_upper_test = out$st2$np_out$harm_upper,
    NP_benefit_capture_test = out$st2$np_out$det
  )
}

summary <- rbind(
  make_single(
    sim_fun = sim_no_hte,
    label = "Constant benefit without HTE",
    seed = 20260516,
    output = "supplementary_figure_s1_constant_benefit_no_hte.png"
  ),
  make_single(
    sim_fun = sim_strong_qual_hte,
    label = "Strong qualitative HTE",
    seed = 20260517,
    output = "supplementary_figure_s2_strong_qualitative_hte.png"
  )
)

write.csv(summary, "result/supplementary_single_replicate_summary.csv", row.names = FALSE)
print(summary)
