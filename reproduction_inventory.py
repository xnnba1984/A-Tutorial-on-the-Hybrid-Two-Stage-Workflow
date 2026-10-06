"""Shared, versioned artifact contract and fail-closed provenance checks."""

from contextlib import contextmanager
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import tempfile
import time

SCHEMA_VERSION = 1
STAGES = {
    "environment": ("R", "check_environment.R"),
    "data": ("R", "prepare_data.R"),
    "main": ("R", "sim_1.R"),
    "truth": ("R", "sim_truth_summary.R"),
    "sensitivity": ("R", "sim_sensitivity.R"),
    "actg_crossfit": ("R", "real_data.R"),
    "actg_validation": ("R", "actg_validation_uncertainty.R"),
    "tables": ("Python", "prepare_revision_outputs.py"),
    "figure2": ("R", "make_simulation_figure2.R"),
    "supplementary_figures": ("R", "make_supplementary_single_replicate_figures.R"),
    "figure1": ("Python", "make_framework_figure1.py"),
}
SETTINGS = {
    "full": {
        "SIM_B_MAIN": "500", "SIM_N_MAIN": "2000", "SIM_TRUTH_N": "1000000",
        "SBR_SENS_B": "300", "SBR_SENS_TREES": "1000", "SBR_SENS_CORES": "8",
        "ACTG_VALIDATION_B": "200", "ACTG_STAGE1_BOOT_B": "1000",
        "ACTG_VALIDATION_TREES": "1500", "ACTG_VALIDATION_THREADS": "8",
        "SUPP_SINGLE_N": "2000",
    },
    "quick": {
        "SIM_B_MAIN": "2", "SIM_N_MAIN": "500", "SIM_TRUTH_N": "10000",
        "SBR_SENS_B": "2", "SBR_SENS_TREES": "50", "SBR_SENS_CORES": "2",
        "ACTG_VALIDATION_B": "2", "ACTG_STAGE1_BOOT_B": "100",
        "ACTG_VALIDATION_TREES": "50", "ACTG_VALIDATION_THREADS": "2",
        "SUPP_SINGLE_N": "500",
    },
}
ALIASES = {
    "proceed": "stage1_reject",
    "proceed_rate": "stage1_rejection_rate",
    "se_proceed_rate": "se_stage1_rejection_rate",
}
OUTPUTS = {
    "environment": ("result/package_versions.csv", "result/session_info.txt"),
    "data": ("data/ACTG175.csv",),
    "main": ("result/sim_main_SiM.csv", "result/sim_main_SiM_summary.csv"),
    "truth": ("result/sim_truth_summary_step6.csv",),
    "sensitivity": tuple("result/information_sensitivity/" + name for name in (
        "sim_information_sensitivity_raw.csv", "sim_information_sensitivity_summary.csv",
        "sim_information_sensitivity_session_info.txt")),
    "actg_crossfit": tuple("result/actg_crossfit/" + name for name in (
        "actg175_ipcw_summary.csv", "actg175_ipcw_optionC.csv", "actg175_ipcw_stepp_cd4.csv",
        "actg175_ipcw_stepp_karnofsky.csv", "actg175_ipcw_uplift_curve.csv",
        "actg175_ipcw_policy_curve.csv", "actg175_ipcw_np_curve.csv", "actg175_ipcw_eval.csv",
        "actg175_ipcw_session_info.txt", "figures/actg175_figure3_ipcw.png")),
    "actg_validation": tuple("result/actg_validation/" + name for name in (
        "actg175_ipcw_stage1_effect_intervals.csv", "actg175_ipcw_stage1_global.csv",
        "actg175_ipcw_stage1_bootstrap_coefficients.csv",
        "actg175_ipcw_stage1_bootstrap_covariance.csv", "actg175_ipcw_validation_splits.csv",
        "actg175_ipcw_validation_summary.csv", "actg175_ipcw_validation_session_info.txt")),
    "tables": ("result/information_sensitivity/supplementary_table_s2.csv",
               "result/revision_analysis_audit.md"),
    "figure2": ("result/figures/Figure_2_revised.png",),
    "supplementary_figures": (
        "result/supplementary_single_replicate_summary.csv",
        "result/figures/supplementary_figure_s1_constant_benefit_no_hte.png",
        "result/figures/supplementary_figure_s2_strong_qualitative_hte.png"),
    "figure1": ("result/figures/Figure_1_revised.png",),
}
DEPENDENCIES = {
    "environment": (), "data": ("environment",), "main": ("environment",),
    "truth": ("environment",), "sensitivity": ("environment",),
    "actg_crossfit": ("environment", "data"),
    "actg_validation": ("environment", "data"),
    "tables": ("sensitivity", "actg_validation"), "figure2": ("main",),
    "supplementary_figures": ("environment",), "figure1": ("environment",),
}
REFERENCE_FILES = frozenset(name.removeprefix("result/")
    for stage, names in OUTPUTS.items() if stage not in ("environment", "data")
    for name in names if name.endswith(".csv"))
FINGERPRINT_ONLY = "actg_crossfit/actg175_ipcw_eval.csv"
SOURCE_FILES = tuple(script for _, script in STAGES.values()) + (
    "reproduce.py", "verify_results.py", "reproduction_inventory.py",
    "requirements.txt", "renv.lock", "environment_versions.csv",
)


class AuditError(Exception):
    """An actionable failed precondition, never a successful verification."""


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_path(root, relative):
    root = Path(root).resolve()
    path = root / relative
    if Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise AuditError(f"Unsafe path: {relative}")
    part = root
    for component in Path(relative).parts:
        part = part / component
        if part.is_symlink():
            raise AuditError(f"Symlinked artifact/source path is not supported: {path}")
    return path


def read_json(path):
    def reject_constant(value):
        raise AuditError(f"Non-finite JSON value: {value}")

    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise AuditError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    try:
        return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique,
                          parse_constant=reject_constant)
    except (OSError, ValueError) as exc:
        raise AuditError(f"Cannot read JSON {path}: {exc}") from exc


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix="." + path.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


@contextmanager
def result_lock(result):
    """Never steal a lock: a killed parent may have left an active R child."""
    result = Path(result)
    lock = result.parent / ("." + result.name + ".reproduction.lock")
    try:
        handle = lock.open("x", encoding="utf-8")
    except FileExistsError as exc:
        raise AuditError(f"Results are locked: {lock}. Another writer may be active. "
                         "Remove a leftover lock only after confirming its process and children stopped.") from exc
    try:
        with handle:
            json.dump({"pid": os.getpid(), "host": socket.gethostname(), "started": time.time()}, handle)
            handle.flush()
            os.fsync(handle.fileno())
        yield
    finally:
        lock.unlink(missing_ok=True)


def verdict(result, status, **details):
    # A rejected check on an absent result tree must not make a fresh run nonempty.
    if Path(result).is_dir():
        atomic_json(Path(result) / "verification.json",
                    {"schema_version": SCHEMA_VERSION, "status": status, **details})


def load_references(root):
    references = Path(root) / "reference_results"
    manifest = read_json(checked_path(root, "reference_results/manifest.json"))
    if not isinstance(manifest, dict) or not {"files", "legacy_column_aliases"} <= manifest.keys():
        raise AuditError("Reference manifest must contain files and legacy_column_aliases")
    if set(manifest) - {"files", "legacy_column_aliases", "description"}:
        raise AuditError("Unknown reference manifest fields")
    if "description" in manifest and not isinstance(manifest["description"], str):
        raise AuditError("Reference description must be text")
    if manifest["legacy_column_aliases"] != ALIASES:
        raise AuditError("Reference aliases differ from the supported three-column mapping")
    entries = manifest["files"]
    if not isinstance(entries, list) or not entries:
        raise AuditError("Reference manifest files must be a nonempty list")
    names = []
    for entry in entries:
        if (not isinstance(entry, dict) or not {"file", "sha256", "comparison"} <= entry.keys()
                or not isinstance(entry["file"], str)
                or not isinstance(entry["sha256"], str)
                or re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]) is None):
            raise AuditError("Malformed reference manifest entry")
        name = entry["file"]
        if name not in REFERENCE_FILES or name in names:
            raise AuditError(f"Unexpected or duplicate reference path: {name}")
        names.append(name)
        if name == FINGERPRINT_ONLY:
            if (entry["comparison"] != "sha256"
                    or set(entry) - {"file", "sha256", "comparison", "rows", "columns", "reason"}
                    or any(type(entry.get(key)) is not int or entry[key] < 1 for key in ("rows", "columns"))
                    or ("reason" in entry and not isinstance(entry["reason"], str))):
                raise AuditError("Participant-level reference requires sha256 mode and positive row/column counts")
            if (references / name).exists():
                raise AuditError("Participant-level reference CSV must not be distributed in reference_results")
        else:
            if entry["comparison"] != "cells" or set(entry) != {"file", "sha256", "comparison"}:
                raise AuditError(f"Only the participant-level evaluation may use fingerprint mode: {name}")
            if sha256(checked_path(root, "reference_results/" + name)) != entry["sha256"]:
                raise AuditError(f"Reference checksum mismatch: {name}")
            read_csv(references / name, ALIASES)
    if set(names) != REFERENCE_FILES:
        raise AuditError("Reference manifest must cover all 21 required CSVs; missing: "
                         + ", ".join(sorted(REFERENCE_FILES - set(names))))
    return manifest


def read_csv(path, aliases=None):
    aliases = aliases or {}
    try:
        with Path(path).open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle, strict=True)
            header = next(reader, None)
            if not header:
                raise AuditError(f"Missing CSV header: {path}")
            columns = [aliases.get(column, column) for column in header]
            if len(set(columns)) != len(columns):
                raise AuditError(f"Duplicate or alias-colliding CSV header: {path}")
            rows = []
            for row in reader:
                if len(row) != len(columns):
                    raise AuditError(f"Ragged CSV row at line {reader.line_num}: {path}")
                rows.append(dict(zip(columns, row)))
            if not rows:
                raise AuditError(f"CSV has no data rows: {path}")
            return columns, rows
    except (csv.Error, UnicodeError) as exc:
        raise AuditError(f"Malformed CSV {path}: {exc}") from exc


def artifact_hashes(root, stage):
    hashes = {}
    for name in OUTPUTS[stage]:
        path = checked_path(root, name)
        if not path.is_file() or path.stat().st_size == 0:
            raise AuditError(f"Missing or empty required artifact: {name}")
        if path.suffix == ".csv":
            read_csv(path)
        elif path.suffix == ".png":
            from PIL import Image
            with Image.open(path) as picture:
                if picture.format != "PNG":
                    raise AuditError(f"Required figure is not PNG: {name}")
                picture.load()
                if not any(low != high for low, high in picture.convert("RGB").getextrema()):
                    raise AuditError(f"Required figure is blank: {name}")
        hashes[name] = sha256(path)
    return hashes


def source_hashes(root):
    return {name: sha256(checked_path(root, name)) for name in SOURCE_FILES}


def affected_stages(selected):
    affected = set(selected)
    for stage in STAGES:
        if set(DEPENDENCIES[stage]) & affected:
            affected.add(stage)
    return affected


def stage_inputs(run, stage):
    inputs, attempts = {}, {}
    for dependency in DEPENDENCIES[stage]:
        record = run["stages"][dependency]
        if record["status"] != "complete":
            raise AuditError(f"{stage} requires completed same-run stage: {dependency}")
        inputs.update(record["outputs"])
        attempts[dependency] = record["attempt_id"]
    return inputs, attempts


def validate_run(root, run, complete=False, ignored=()):
    required = {"schema_version", "run_id", "root", "mode", "status", "settings", "sources",
                "reference_manifest_sha256", "stages", "runtime", "data_source_tarball"}
    if (not isinstance(run, dict) or set(run) != required
            or type(run.get("schema_version")) is not int or run.get("schema_version") != SCHEMA_VERSION
            or run.get("root") != str(Path(root).resolve())
            or re.fullmatch(r"[0-9a-f]{32}", str(run.get("run_id", ""))) is None
            or run.get("mode") not in SETTINGS
            or run.get("status") not in {"running", "failed", "incomplete", "complete"}):
        raise AuditError("Missing, old-format, or malformed run_manifest.json; start a fresh run")
    runtime = run["runtime"]
    if (not isinstance(runtime, dict) or set(runtime) != {
            "python", "python_executable", "pillow", "fonts", "r_launcher",
            "r_launcher_sha256", "r_script_first_wrapper", "r_environment"}
            or not isinstance(runtime.get("fonts"), dict)
            or set(runtime["fonts"]) != {"regular", "bold", "small"}
            or type(runtime.get("r_script_first_wrapper")) is not bool
            or not isinstance(runtime.get("r_environment"), dict)
            or any(not isinstance(runtime.get(key), str) or not runtime[key]
                   for key in ("python", "python_executable", "pillow", "r_launcher"))):
        raise AuditError("Malformed or missing Python/font/R preflight record")
    fingerprints = [*runtime["fonts"].values()]
    if run["data_source_tarball"] is not None:
        fingerprints.append(run["data_source_tarball"])
    for fingerprint in fingerprints:
        if (not isinstance(fingerprint, dict) or set(fingerprint) != {"path", "sha256"}
                or not isinstance(fingerprint["path"], str) or not fingerprint["path"]
                or re.fullmatch(r"[0-9a-f]{64}", str(fingerprint["sha256"])) is None):
            raise AuditError("Malformed font/data-source fingerprint")
    if re.fullmatch(r"[0-9a-f]{64}", str(runtime["r_launcher_sha256"])) is None:
        raise AuditError("Malformed R launcher fingerprint")
    if checked_path(root, "result/run_mode.txt").read_text().strip() != run["mode"]:
        raise AuditError("Run mode marker and run manifest disagree")
    settings = run.get("settings")
    if not isinstance(settings, dict) or set(settings) != set(SETTINGS[run["mode"]]):
        raise AuditError("Malformed run settings")
    for key, value in settings.items():
        if key in {"SBR_SENS_CORES", "ACTG_VALIDATION_THREADS"}:
            if not isinstance(value, str) or not value.isascii() or not value.isdecimal() or int(value) < 1:
                raise AuditError("Worker counts must be positive integers")
        elif value != SETTINGS[run["mode"]][key]:
            raise AuditError(f"Scientific setting differs from {run['mode']}: {key}")
    if settings["SBR_SENS_CORES"] != settings["ACTG_VALIDATION_THREADS"]:
        raise AuditError("Inconsistent worker settings")
    if run.get("sources") != source_hashes(root):
        raise AuditError("Source/script/dependency specification hashes changed; use a fresh run")
    if run.get("reference_manifest_sha256") != sha256(Path(root) / "reference_results/manifest.json"):
        raise AuditError("Reference manifest changed since this run started")
    records = run.get("stages")
    if not isinstance(records, dict) or set(records) != set(STAGES):
        raise AuditError("Run manifest must contain all required stages")
    for stage, record in records.items():
        if not isinstance(record, dict) or record.get("status") not in {
                "pending", "invalidated", "running", "failed", "complete"}:
            raise AuditError(f"Malformed stage record: {stage}")
        if stage in ignored:
            continue
        if record["status"] != "complete":
            if complete:
                raise AuditError(f"Required stage is not complete: {stage}")
            continue
        inputs, attempts = stage_inputs(run, stage)
        if (record.get("run_id") != run["run_id"] or record.get("exit_code") != 0
                or type(record.get("exit_code")) is not int
                or re.fullmatch(r"[0-9a-f]{32}", str(record.get("attempt_id", ""))) is None
                or record.get("script") != STAGES[stage][1]
                or record.get("script_sha256") != run["sources"][STAGES[stage][1]]
                or record.get("settings") != settings or record.get("inputs") != inputs
                or record.get("dependency_attempts") != attempts):
            raise AuditError(f"Invalid or mixed-run stage provenance: {stage}")
        if record.get("outputs") != artifact_hashes(root, stage):
            raise AuditError(f"Output hash mismatch: {stage}")
        if record.get("log_sha256") != sha256(checked_path(root, f"result/logs/{stage}.log")):
            raise AuditError(f"Stage log hash mismatch: {stage}")
        if read_json(checked_path(root, f"result/logs/{stage}.json")) != record:
            raise AuditError(f"Stage log record and run manifest disagree: {stage}")
    if complete and run["status"] != "complete":
        raise AuditError(f"Run is not complete: {run['status']}")
