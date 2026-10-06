#!/usr/bin/env python3
"""Fresh reproduction by default; --stages resumes a checked, same-run result tree."""

import argparse
import ast
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

sys.dont_write_bytecode = True
from reproduction_inventory import (
    AuditError, OUTPUTS, SCHEMA_VERSION, SETTINGS, STAGES, affected_stages,
    artifact_hashes, atomic_json, checked_path, load_references, read_json,
    result_lock, sha256, source_hashes, stage_inputs, validate_run, verdict,
)

ROOT = Path(__file__).resolve().parent
OUTPUT_OVERRIDES = ("SBR_SENS_OUT_DIR", "ACTG_REALDATA_OUT_DIR", "ACTG_VALIDATION_OUT_DIR",
                    "SBR_REVISION_AUDIT_OUT", "SUPP_FIGURE_OUT_DIR", "SUPP_SUMMARY_OUT")
STARTUP_FILES = ("R_PROFILE", "R_PROFILE_USER", "R_ENVIRON", "R_ENVIRON_USER")


def execution_environment(root, settings):
    env = os.environ.copy()
    env.update(settings)
    for key in OUTPUT_OVERRIDES:
        env.pop(key, None)
    for key in STARTUP_FILES:
        env[key] = os.devnull
    env["R_DEFAULT_PACKAGES"] = "datasets,utils,grDevices,graphics,stats,methods"
    env["R_TESTS"] = ""
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if (root / "renv/library").is_dir():
        env["R_LIBS_USER"] = str(root / "renv/library")
    return env


def preflight(root, env):
    if sys.version_info < (3, 9):
        raise AuditError("Python 3.9 or newer is required")
    try:
        import PIL
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise AuditError("Pillow is required; install requirements.txt using this Python interpreter") from exc
    for name in ("reproduce.py", "verify_results.py", "reproduction_inventory.py",
                 "prepare_revision_outputs.py", "make_framework_figure1.py"):
        compile((root / name).read_text(encoding="utf-8"), name, "exec")
    # Execute only the existing font lookup function, never the figure's top-level drawing code.
    tree = ast.parse((root / "make_framework_figure1.py").read_text(encoding="utf-8"))
    lookup = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "load_font"]
    if len(lookup) != 1:
        raise AuditError("Cannot locate Figure 1 font lookup for preflight")
    namespace = {"os": os, "Path": Path, "ImageFont": ImageFont}
    exec(compile(ast.Module(body=lookup, type_ignores=[]), "font-preflight", "exec"), namespace)
    fonts = {}
    for label, size, bold in (("regular", 59, False), ("bold", 72, True), ("small", 52, False)):
        font = namespace["load_font"](size, bold=bold)
        ImageDraw.Draw(Image.new("RGB", (16, 16))).multiline_textbbox((0, 0), "Preflight", font=font)
        path = Path(font.path).resolve()
        fonts[label] = {"path": str(path), "sha256": sha256(path)}
    launcher = shutil.which(env.get("R_BIN", "Rscript"), path=env.get("PATH"))
    if launcher is None:
        raise AuditError("R launcher not found; set R_BIN to one executable or script-first wrapper path")
    launcher = str(Path(launcher).absolute())
    return {
        "python": sys.version, "python_executable": str(Path(sys.executable).absolute()),
        "pillow": PIL.__version__, "fonts": fonts, "r_launcher": launcher,
        "r_launcher_sha256": sha256(launcher), "r_script_first_wrapper": "R_BIN" in env,
        "r_environment": {key: env.get(key) for key in (
            *STARTUP_FILES, "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_DEFAULT_PACKAGES", "R_TESTS")},
    }


def run_stage(root, run, stage, env):
    inputs, attempts = stage_inputs(run, stage)
    if source_hashes(root) != run["sources"]:
        raise AuditError("Source hashes changed while the run was active")
    language, script = STAGES[stage]
    if language == "R":
        command = [run["runtime"]["r_launcher"]]
        if not run["runtime"]["r_script_first_wrapper"]:
            command.append("--vanilla")
        command.append(script)
    else:
        command = [sys.executable, "-B", script]
    if stage == "data" and run["data_source_tarball"] is not None:
        archive = run["data_source_tarball"]
        if sha256(archive["path"]) != archive["sha256"]:
            raise AuditError("Offline data source archive changed")
        command.extend(["--source-tarball", archive["path"]])
    record = {
        "run_id": run["run_id"], "attempt_id": uuid.uuid4().hex, "status": "running",
        "script": script, "script_sha256": run["sources"][script], "settings": run["settings"],
        "command": command, "inputs": inputs, "dependency_attempts": attempts, "outputs": {},
    }
    run["stages"][stage] = record
    result = root / "result"
    atomic_json(result / "run_manifest.json", run)
    atomic_json(result / f"logs/{stage}.json", record)
    started = time.monotonic()
    log = result / f"logs/{stage}.log"
    error = None
    try:
        # Only this stage's generated results are replaced. R revalidates the canonical input CSV.
        for name in OUTPUTS[stage]:
            if name.startswith("result/"):
                checked_path(root, name).unlink(missing_ok=True)
        print(f"[{run['mode']}] Starting {stage}", flush=True)
        with log.open("w", encoding="utf-8") as handle:
            process = subprocess.Popen(command, cwd=root, env=env, stdout=handle,
                                       stderr=subprocess.STDOUT, start_new_session=os.name == "posix")
            try:
                record["exit_code"] = process.wait()
            except BaseException:
                stop_process(process)
                raise
        if record["exit_code"] != 0:
            raise AuditError(f"{stage} exited {record['exit_code']}; inspect {log}")
        record["outputs"] = artifact_hashes(root, stage)
        if source_hashes(root) != run["sources"]:
            raise AuditError("Source hashes changed during stage execution")
        record["status"] = "complete"
    except (Exception, KeyboardInterrupt) as exc:
        error = exc
        record.update(status="failed", error=str(exc) or "Interrupted")
        run["status"] = "failed"
    finally:
        record["elapsed_seconds"] = round(time.monotonic() - started, 3)
        if log.is_file():
            record["log_sha256"] = sha256(log)
        atomic_json(result / f"logs/{stage}.json", record)
        atomic_json(result / "run_manifest.json", run)
    if error is not None:
        raise error
    print(f"[{run['mode']}] Finished {stage}: {record['elapsed_seconds']} seconds", flush=True)


def stop_process(process):
    # Forked R workers share the child's new process group. Stop them before releasing the run lock.
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    else:
        process.terminate()
    try:
        process.wait(timeout=30)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
    finally:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def execute(root, args):
    result = checked_path(root, "result")
    if args.workers is not None and args.workers < 1:
        raise AuditError("--workers must be positive")
    if args.stages is not None and len(set(args.stages)) != len(args.stages):
        raise AuditError("Duplicate --stages entries")
    load_references(root)
    if args.stages is None:
        if result.exists() and any(result.iterdir()):
            raise AuditError("Fresh runs require an empty/absent result/. Use a fresh copy, or --stages for a checked resume")
        mode = args.mode or "full"
        settings = dict(SETTINGS[mode])
        if args.workers is not None:
            settings.update(SBR_SENS_CORES=str(args.workers), ACTG_VALIDATION_THREADS=str(args.workers))
        run = None
        selected = list(STAGES)
    else:
        run = read_json(result / "run_manifest.json")
        affected = affected_stages(args.stages)
        validate_run(root, run, ignored=affected)
        mode, settings = run["mode"], run["settings"]
        if args.mode is not None and args.mode != mode:
            raise AuditError("Do not mix quick and full runs")
        if args.workers is not None and str(args.workers) != settings["ACTG_VALIDATION_THREADS"]:
            raise AuditError("Resume must retain this run's workers; start a fresh run to change settings")
        selected = [stage for stage in STAGES if stage in args.stages]
    env = execution_environment(root, settings)
    runtime = preflight(root, env)
    if run is not None and runtime != run.get("runtime"):
        raise AuditError("Python/font/R launcher or library configuration changed; use the original environment or a fresh run")
    archive = run.get("data_source_tarball") if run else None
    if args.data_source_tarball is not None:
        path = args.data_source_tarball.expanduser().resolve(strict=True)
        proposed = {"path": str(path), "sha256": sha256(path)}
        if archive is not None and proposed["sha256"] != archive["sha256"]:
            raise AuditError("Offline archive differs from this run's data source")
        archive = proposed
    if run is None:
        run = {"schema_version": SCHEMA_VERSION, "run_id": uuid.uuid4().hex, "root": str(root),
               "mode": mode, "settings": settings, "sources": source_hashes(root), "runtime": runtime,
               "reference_manifest_sha256": sha256(root / "reference_results/manifest.json"),
               "stages": {stage: {"status": "pending"} for stage in STAGES}}
        result.mkdir(parents=True, exist_ok=True)
        (result / "run_mode.txt").write_text(mode + "\n", encoding="utf-8")
    else:
        for stage in affected:
            run["stages"][stage]["status"] = "invalidated"
    run.update(status="running", data_source_tarball=archive)
    atomic_json(result / "run_manifest.json", run)
    verdict(result, "RUNNING", run_id=run["run_id"])
    try:
        for stage in selected:
            validate_run(root, run)
            run_stage(root, run, stage, env)
        run["status"] = "complete" if all(r["status"] == "complete" for r in run["stages"].values()) else "incomplete"
        validate_run(root, run, complete=run["status"] == "complete")
    except (Exception, KeyboardInterrupt):
        run["status"] = "failed"
        atomic_json(result / "run_manifest.json", run)
        raise
    atomic_json(result / "run_manifest.json", run)
    remaining = [stage for stage, record in run["stages"].items() if record["status"] != "complete"]
    verdict(result, "UNVERIFIED", run_id=run["run_id"], run_status=run["status"], remaining_stages=remaining)
    print(f"Run {run['run_id']}: {run['status']}. " + (
        "Resume required stages with --stages " + " ".join(remaining) if remaining else "Next: python3 verify_results.py"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=SETTINGS, default=None, help="Fresh default: full; resume default: recorded mode")
    parser.add_argument("--stages", nargs="+", choices=STAGES, help="Resume selected stages in dependency order")
    parser.add_argument("--workers", type=int, help="Fresh-run worker override; resumes retain recorded workers")
    parser.add_argument("--data-source-tarball", type=Path, help="Pinned speff2trial source archive for offline data preparation")
    root = ROOT.resolve()
    try:
        result = checked_path(root, "result")
        with result_lock(result):
            try:
                args = parser.parse_args()
                # Preserve emptiness until execute decides whether this is a fresh run.
                if (result / "verification.json").exists():
                    verdict(result, "UNVERIFIED", error="A new runner invocation is being checked")
                execute(root, args)
            except SystemExit as exc:
                if exc.code:
                    verdict(result, "FAIL", error="Invalid runner arguments")
                raise
            except (Exception, KeyboardInterrupt) as exc:
                verdict(result, "FAIL", error=str(exc) or "Interrupted")
                print(f"FAIL: {exc}", file=sys.stderr)
                return 130 if isinstance(exc, KeyboardInterrupt) else 1
    except (AuditError, OSError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
