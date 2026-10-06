#!/usr/bin/env python3
"""Verify complete same-run provenance, 20 cell comparisons and one exact fingerprint."""

import argparse
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from reproduction_inventory import (
    ALIASES, AuditError, FINGERPRINT_ONLY, OUTPUTS, checked_path, load_references,
    read_csv, read_json, result_lock, sha256, validate_run, verdict,
)

ROOT = Path(__file__).resolve().parent
MISSING = {"", "NA", "NaN", "nan"}
THREAD_SUMMARY = "actg_validation/actg175_ipcw_validation_summary.csv"


def equal_cell(left, right, atol, rtol):
    if left in MISSING or right in MISSING:
        return left in MISSING and right in MISSING, 0.0
    if left == right:
        return True, 0.0
    try:
        a, b = float(left), float(right)
    except (TypeError, ValueError):
        return False, None
    if not math.isfinite(a) or not math.isfinite(b):
        return a == b, None
    difference = abs(a - b)
    return math.isclose(a, b, abs_tol=atol, rel_tol=rtol), difference if math.isfinite(difference) else None


def compare_file(references, actual, entry, workers, atol, rtol):
    name = entry["file"]
    result = {"file": name, "comparison": entry["comparison"], "status": "PASS"}
    try:
        columns, observed = read_csv(actual, ALIASES)
        result.update(rows=len(observed), columns=len(columns), sha256=sha256(actual))
        if name == FINGERPRINT_ONLY:
            result.update(reference_sha256=entry["sha256"], reference_rows=entry["rows"],
                          reference_columns=entry["columns"])
            if (result["sha256"] != entry["sha256"] or len(observed) != entry["rows"]
                    or len(columns) != entry["columns"]):
                result.update(status="FAIL", error="Exact output fingerprint or row/column count mismatch")
            return result
        reference_columns, expected = read_csv(references / name, ALIASES)
        result.update(different_cells=0, max_absolute_difference=0.0, examples=[], metadata=[])
        if columns != reference_columns or len(observed) != len(expected):
            result.update(status="FAIL", error="Row count or ordered column schema differs")
            return result
        for i, (left, right) in enumerate(zip(expected, observed), start=1):
            for column in columns:
                if (name == THREAD_SUMMARY and column == "mean"
                        and left.get("metric") == right.get("metric") == "num_threads"):
                    try:
                        reference_workers, actual_workers = float(left[column]), float(right[column])
                        allowed = (math.isfinite(reference_workers) and reference_workers.is_integer()
                                   and reference_workers > 0 and math.isfinite(actual_workers)
                                   and actual_workers.is_integer() and actual_workers > 0
                                   and workers is not None and actual_workers == int(workers))
                    except (TypeError, ValueError, OverflowError):
                        allowed = False
                    result["metadata"].append({"row": i, "metric": "num_threads", "column": column,
                        "reference": left[column], "rerun": right[column], "matches_run_settings": allowed})
                    if not allowed:
                        result.update(status="FAIL", error="num_threads must be a positive integer matching the run settings")
                    continue
                equal, difference = equal_cell(left[column], right[column], atol, rtol)
                if difference is not None:
                    result["max_absolute_difference"] = max(result["max_absolute_difference"], difference)
                if not equal:
                    result["status"] = "FAIL"
                    result["different_cells"] += 1
                    if len(result["examples"]) < 12:
                        result["examples"].append({"row": i, "column": column,
                                                  "reference": left[column], "rerun": right[column]})
    except (AuditError, OSError, ValueError) as exc:
        result.update(status="FAIL", error=str(exc))
    return result


def verify(root, result, atol, rtol):
    if not all(math.isfinite(value) and value >= 0 for value in (atol, rtol)):
        raise AuditError("Tolerances must be finite and nonnegative")
    manifest = load_references(root)
    run = read_json(result / "run_manifest.json")
    errors = []
    if not isinstance(run, dict) or run.get("mode") != "full":
        raise AuditError("Verification requires a full run, not quick-test outputs")
    try:
        validate_run(root, run, complete=True)
    except (AuditError, OSError, ValueError) as exc:
        errors.append(str(exc))
    settings = run.get("settings")
    workers = settings.get("ACTG_VALIDATION_THREADS") if isinstance(settings, dict) else None
    comparisons = []
    for entry in manifest["files"]:
        actual = checked_path(root, "result/" + entry["file"])
        comparison = compare_file(root / "reference_results", actual, entry, workers, atol, rtol)
        comparisons.append(comparison)
        detail = ("exact fingerprint" if comparison["comparison"] == "sha256" else
                  f"differences={comparison.get('different_cells', 'unavailable')}; "
                  f"max_abs={comparison.get('max_absolute_difference', 'unavailable')}")
        print(f"{comparison['status']}: {entry['file']}; {detail}")
    # Recheck after comparisons so external edits cannot reuse a prior provenance check.
    try:
        validate_run(root, run, complete=True)
        if read_json(result / "run_manifest.json") != run:
            raise AuditError("Run manifest changed during verification")
        if load_references(root) != manifest:
            raise AuditError("Reference manifest changed during verification")
    except (AuditError, OSError, ValueError) as exc:
        if str(exc) not in errors:
            errors.append(str(exc))
    status = "PASS" if not errors and all(item["status"] == "PASS" for item in comparisons) else "FAIL"
    verdict(result, status, run_id=run.get("run_id"), absolute_tolerance=atol, relative_tolerance=rtol,
            errors=errors, files=comparisons,
            coverage={"cell_comparisons": 20, "fingerprint_comparisons": 1,
                      "required_artifacts": sum(map(len, OUTPUTS.values())), "required_figures": 5},
            figure_validation="Required nonblank PNGs, same-run hashes; not an R1 pixel-equivalence claim")
    print(f"Overall: {status}")
    return 0 if status == "PASS" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "result",
                        help="result/ directory in its intact originating repository (no detached/legacy outputs)")
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-10)
    # Locate the destination before parsing comparison arguments so invalid tolerances also replace old PASS.
    locator = argparse.ArgumentParser(add_help=False)
    locator.add_argument("--results", type=Path, default=ROOT / "result")
    location, _ = locator.parse_known_args()
    result = location.results.absolute()
    if result.name != "result":
        print("FAIL: --results must be the originating repository's result/ directory", file=sys.stderr)
        return 1
    root = result.parent.resolve()
    try:
        result = checked_path(root, "result")
        with result_lock(result):
            try:
                args = parser.parse_args()
                verdict(result, "RUNNING", error="Verification in progress")
                return verify(root, result, args.atol, args.rtol)
            except SystemExit as exc:
                if exc.code:
                    verdict(result, "FAIL", error="Invalid verifier arguments")
                raise
            except (Exception, KeyboardInterrupt) as exc:
                verdict(result, "FAIL", error=str(exc) or "Interrupted")
                print(f"FAIL: {exc}", file=sys.stderr)
                return 130 if isinstance(exc, KeyboardInterrupt) else 1
    except (AuditError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
