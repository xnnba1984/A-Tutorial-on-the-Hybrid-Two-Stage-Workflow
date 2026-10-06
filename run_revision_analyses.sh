#!/bin/sh
set -eu

# Reproduce the analyses added or corrected for the SBR revision: the focused
# information sensitivity study, ACTG Stage 1 bootstrap inference, independent
# repeated-split policy validation, Supplementary Table S2, and an audit file.

REPO_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
R_BIN=${R_BIN:-Rscript}
PYTHON_BIN=${PYTHON_BIN:-python3}

SENS_OUT_DIR=${SBR_SENS_OUT_DIR:-"$REPO_ROOT/result/information_sensitivity"}
SENS_B=${SBR_SENS_B:-300}
SENS_TREES=${SBR_SENS_TREES:-1000}
SENS_CORES=${SBR_SENS_CORES:-8}

ACTG_OUT_DIR=${ACTG_VALIDATION_OUT_DIR:-"$REPO_ROOT/result/actg_validation"}
ACTG_B=${ACTG_VALIDATION_B:-200}
ACTG_BOOT_B=${ACTG_STAGE1_BOOT_B:-1000}
ACTG_TREES=${ACTG_VALIDATION_TREES:-1500}
ACTG_THREADS=${ACTG_VALIDATION_THREADS:-8}
AUDIT_OUT=${SBR_REVISION_AUDIT_OUT:-"$REPO_ROOT/result/revision_analysis_audit.md"}

cd "$REPO_ROOT"

SBR_SENS_OUT_DIR="$SENS_OUT_DIR" \
SBR_SENS_B="$SENS_B" \
SBR_SENS_TREES="$SENS_TREES" \
SBR_SENS_CORES="$SENS_CORES" \
"$R_BIN" sim_sensitivity.R

ACTG_VALIDATION_OUT_DIR="$ACTG_OUT_DIR" \
ACTG_VALIDATION_B="$ACTG_B" \
ACTG_STAGE1_BOOT_B="$ACTG_BOOT_B" \
ACTG_VALIDATION_TREES="$ACTG_TREES" \
ACTG_VALIDATION_THREADS="$ACTG_THREADS" \
"$R_BIN" actg_validation_uncertainty.R

SBR_SENS_OUT_DIR="$SENS_OUT_DIR" \
SBR_SENS_B="$SENS_B" \
SBR_SENS_TREES="$SENS_TREES" \
SBR_SENS_CORES="$SENS_CORES" \
ACTG_VALIDATION_OUT_DIR="$ACTG_OUT_DIR" \
ACTG_VALIDATION_B="$ACTG_B" \
ACTG_STAGE1_BOOT_B="$ACTG_BOOT_B" \
ACTG_VALIDATION_TREES="$ACTG_TREES" \
ACTG_VALIDATION_THREADS="$ACTG_THREADS" \
SBR_REVISION_AUDIT_OUT="$AUDIT_OUT" \
"$PYTHON_BIN" prepare_revision_outputs.py
