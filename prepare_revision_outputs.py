#!/usr/bin/env python3
"""Prepare Supplementary Table S2 and an audit summary for a repository run."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SENS_DIR = Path(os.environ.get(
    "SBR_SENS_OUT_DIR",
    ROOT / "result" / "information_sensitivity",
))
ACTG_DIR = Path(os.environ.get(
    "ACTG_VALIDATION_OUT_DIR",
    ROOT / "result" / "actg_validation",
))
AUDIT_OUT = Path(os.environ.get(
    "SBR_REVISION_AUDIT_OUT",
    ROOT / "result" / "revision_analysis_audit.md",
))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fmt(value: str, digits: int = 3) -> str:
    cleaned = value.strip()
    if cleaned.upper() == "NA" or not cleaned:
        return "NA"
    numeric = float(cleaned)
    return f"{numeric:.{digits}f}" if math.isfinite(numeric) else "NA"


def build_supplementary_table() -> Path:
    rows = read_csv(SENS_DIR / "sim_information_sensitivity_summary.csv")
    by_key = {
        (r["scenario"], int(r["n"]), int(r["p"]), r["learner"]): r
        for r in rows
    }
    scenario_order = [
        "Constant benefit without HTE",
        "Strong quantitative HTE",
        "Strong qualitative HTE",
    ]
    output_rows: list[dict[str, str]] = []
    for scenario in scenario_order:
        for n in (500, 1000, 2000):
            for p in (3, 20):
                cf = by_key[(scenario, n, p, "causal_forest")]
                linear = by_key[(scenario, n, p, "linear_interaction")]
                output_rows.append(
                    {
                        "Scenario": scenario,
                        "n": str(n),
                        "p": str(p),
                        "Stage 1 rejection rate (SE)": (
                            f"{fmt(cf['stage1_rejection_rate'])} "
                            f"({fmt(cf['se_stage1_rejection_rate'])})"
                        ),
                        "Causal forest true gain vs best fixed (SE)": (
                            f"{fmt(cf['mean_true_incremental_value_best_fixed'])} "
                            f"({fmt(cf['se_true_incremental_value_best_fixed'])})"
                        ),
                        "Linear interaction true gain vs best fixed (SE)": (
                            f"{fmt(linear['mean_true_incremental_value_best_fixed'])} "
                            f"({fmt(linear['se_true_incremental_value_best_fixed'])})"
                        ),
                    }
                )

    output_path = SENS_DIR / "supplementary_table_s2.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    return output_path


def get_metric(rows: list[dict[str, str]], metric: str) -> dict[str, str]:
    for row in rows:
        if row["metric"] == metric:
            return row
    raise KeyError(metric)


def write_audit_summary(s2_path: Path) -> Path:
    sensitivity = read_csv(SENS_DIR / "sim_information_sensitivity_summary.csv")
    actg = read_csv(ACTG_DIR / "actg175_ipcw_validation_summary.csv")
    gain = get_metric(actg, "test_value_gain_vs_selected_fixed")
    value = get_metric(actg, "test_threshold_value")
    fixed = get_metric(actg, "test_treat_all_value")
    auqc = get_metric(actg, "test_centered_auqc")
    positive = get_metric(actg, "value_gain_positive_rate")
    feasible = get_metric(actg, "np_feasible_tune_rate")
    sensitivity_b_values = sorted({int(float(row["B"])) for row in sensitivity})
    sensitivity_b = ", ".join(str(value) for value in sensitivity_b_values)
    sensitivity_trees = int(os.environ.get("SBR_SENS_TREES", "1000"))
    sensitivity_cores = int(os.environ.get("SBR_SENS_CORES", "8"))
    actg_b = int(float(get_metric(actg, "B")["mean"]))
    actg_trees = int(float(get_metric(actg, "num_trees")["mean"]))
    actg_threads = int(float(get_metric(actg, "num_threads")["mean"]))
    actg_bootstrap_b = int(os.environ.get("ACTG_STAGE1_BOOT_B", "1000"))

    lines = [
        "# Final Revision Analysis Audit",
        "",
        "## Information sensitivity simulation",
        "",
        f"- Run design: {sensitivity_b} replicate(s) per setting, n = 500, 1000, or 2000; p = 3 or 20; three representative scenarios; causal forest and linear interaction learners evaluated on the same train, tune, and test splits.",
        f"- Computation: {sensitivity_trees} trees for each causal forest and {sensitivity_cores} worker process(es).",
        "- The p = 20 setting retains the three original covariates and adds 17 noise covariates. Stage 1 jointly tests all p treatment-by-covariate interactions.",
        f"- Manuscript-ready table: `{s2_path.relative_to(ROOT)}`.",
        "",
        "## ACTG 175 repeated split validation",
        "",
        f"- Run design: {actg_b} stratified 50/25/25 train, tune, and test split(s), {actg_bootstrap_b} Stage 1 bootstrap replicate(s), {actg_trees} causal-forest trees per split, and {actg_threads} worker thread(s).",
        "- Within each split, one model bundle is fitted on the training set and reused without refitting for tuning and test predictions.",
        f"- Learned policy mean test value: {float(value['mean']):.3f}; stability percentiles {float(value['p025']):.3f} to {float(value['p975']):.3f}.",
        f"- Treat-all mean test value: {float(fixed['mean']):.3f}; stability percentiles {float(fixed['p025']):.3f} to {float(fixed['p975']):.3f}.",
        f"- Mean gain over treat all: {float(gain['mean']):.4f}; stability percentiles {float(gain['p025']):.4f} to {float(gain['p975']):.4f}.",
        f"- Positive gain rate: {100 * float(positive['mean']):.1f}%.",
        f"- Mean cAUQC: {float(auqc['mean']):.3f}; stability percentiles {float(auqc['p025']):.3f} to {float(auqc['p975']):.3f}.",
        f"- The 0.10 surrogate harm constraint was feasible on the tuning split in {100 * float(feasible['mean']):.1f}% of splits.",
        "- Interpret the learned-policy comparison from the mean gain, stability percentiles, and positive-gain rate rather than from any single split.",
        "",
    ]
    output_path = AUDIT_OUT
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path


def main() -> None:
    s2_path = build_supplementary_table()
    audit_path = write_audit_summary(s2_path)
    print(s2_path)
    print(audit_path)


if __name__ == "__main__":
    main()
