suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

root <- normalizePath(getwd())
summary_path <- file.path(root, "result", "sim_main_SiM_summary.csv")
out_dir <- file.path(root, "result", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(summary_path, stringsAsFactors = FALSE)

# Accept the archived output names while using the revised terminology below.
if (!("stage1_rejection_rate" %in% names(dat)) && "proceed_rate" %in% names(dat)) {
  dat$stage1_rejection_rate <- dat$proceed_rate
}
if (!("se_stage1_rejection_rate" %in% names(dat)) && "se_proceed_rate" %in% names(dat)) {
  dat$se_stage1_rejection_rate <- dat$se_proceed_rate
}

scenario_order <- c(
  "No treatment effect",
  "No HTE",
  "Weak quantitative HTE",
  "Strong quantitative HTE",
  "Weak qualitative HTE",
  "Strong qualitative HTE"
)

scenario_labels <- c(
  "No\neffect",
  "Constant\nbenefit",
  "Weak\nquant.",
  "Strong\nquant.",
  "Weak\nqual.",
  "Strong\nqual."
)

dat <- dat |>
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    scenario_label = factor(
      scenario_labels[match(as.character(scenario), scenario_order)],
      levels = scenario_labels
    ),
    group = ifelse(
      grepl("qualitative", as.character(scenario)),
      "Qualitative HTE",
      ifelse(
        grepl("quantitative", as.character(scenario)),
        "Quantitative HTE",
        "Negative control"
      )
    )
  )

fill_values <- c(
  "Negative control" = "#d9d9d9",
  "Quantitative HTE" = "#9ecae1",
  "Qualitative HTE" = "#74c476"
)

theme_sim <- function() {
  theme_bw(base_size = 12.5, base_family = "Times") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.4),
      axis.text.x = element_text(
        angle = 42,
        hjust = 1,
        vjust = 1,
        color = "black",
        size = 9.2
      ),
      axis.text.y = element_text(color = "black", size = 10.5),
      axis.title = element_text(color = "black", size = 11.5),
      plot.title = element_text(face = "bold", hjust = 0, size = 14),
      legend.position = "none",
      plot.margin = margin(7, 10, 6, 7)
    )
}

p1 <- ggplot(dat, aes(x = scenario_label, y = stage1_rejection_rate, fill = group)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(
      ymin = pmax(0, stage1_rejection_rate - se_stage1_rejection_rate),
      ymax = pmin(1.05, stage1_rejection_rate + se_stage1_rejection_rate)
    ),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = fill_values) +
  labs(title = "A", x = NULL, y = "Stage 1 rejection rate") +
  theme_sim()

p2 <- ggplot(dat, aes(x = scenario_label, y = mean_valuegain, fill = group)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(
      ymin = mean_valuegain - se_valuegain,
      ymax = mean_valuegain + se_valuegain
    ),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  scale_y_continuous(
    limits = c(-0.012, 0.043),
    breaks = seq(-0.01, 0.04, 0.01)
  ) +
  scale_fill_manual(values = fill_values) +
  labs(title = "B", x = NULL, y = "Estimated incremental value") +
  theme_sim()

fig <- gridExtra::arrangeGrob(p1, p2, ncol = 2)

ggsave(
  filename = file.path(out_dir, "Figure_2_revised.png"),
  plot = fig,
  width = 7.2,
  height = 4.1,
  dpi = 450,
  bg = "white"
)
