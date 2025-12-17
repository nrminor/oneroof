#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "flask",
#     "werkzeug",
# ]
# ///
"""
OneRoof Globus Action Provider

This action provider implements the custom actions needed for the OneRoof
bioinformatics pipeline to work with Globus Flows.
"""

import logging
import os
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

from flask import Flask, jsonify, request
from werkzeug.exceptions import BadRequest

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)


def load_config() -> dict[str, Any]:
    """Load configuration from environment variables or config file."""
    config = {
        "NEXTFLOW_PATH": os.environ.get("NEXTFLOW_PATH", "nextflow"),
        "ONEROOF_PATH": os.environ.get(
            "ONEROOF_PATH",
            os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            ),
        ),
        "WORK_BASE_DIR": os.environ.get("ONEROOF_WORK_DIR", "/tmp/oneroof_runs"),
        "LOG_LEVEL": os.environ.get("LOG_LEVEL", "INFO"),
        "MAX_PARALLEL_RUNS": int(os.environ.get("MAX_PARALLEL_RUNS", "5")),
        "ENABLE_GPU": os.environ.get("ENABLE_GPU", "true").lower() == "true",
        "PIXI_PATH": os.environ.get("PIXI_PATH", "pixi"),
        "USE_PIXI": os.environ.get("USE_PIXI", "false").lower() == "true",
    }

    # Try to load from .env file if it exists
    env_file = Path(__file__).parent / ".env"
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    key, value = line.strip().split("=", 1)
                    if key not in os.environ:
                        config[key] = value

    return config


CONFIG = load_config()


@app.route("/", methods=["POST"])
def handle_action():
    """Main action handler - routes to specific action functions."""
    try:
        data = request.get_json()
        if not data:
            raise BadRequest("No JSON data provided")

        action = data.get("action")
        if not action:
            raise BadRequest("No action specified")

        # Route to appropriate handler
        handlers = {
            "validate": validate_input,
            "prepare_workspace": prepare_workspace,
            "run_pipeline": run_pipeline,
            "collect_logs": collect_logs,
            "cleanup": cleanup_workspace,
            "status": check_status,
        }

        handler = handlers.get(action)
        if not handler:
            raise BadRequest(f"Unknown action: {action}")

        result = handler(data)
        return jsonify(result)

    except Exception as e:
        logger.error(f"Action failed: {e!s}")
        return jsonify({"error": str(e)}), 500


def validate_input(data: dict[str, Any]) -> dict[str, Any]:
    """Validate input parameters for the pipeline."""
    platform = data.get("platform")
    if platform not in ["nanopore", "illumina"]:
        raise ValueError(f"Invalid platform: {platform}")

    required_files = ["primer_bed", "refseq", "ref_gbk"]
    missing = [f for f in required_files if not data.get(f)]
    if missing:
        raise ValueError(f"Missing required files: {missing}")

    input_data = data.get("input_data")
    if not input_data:
        raise ValueError("No input data specified")

    return {
        "status": "validated",
        "platform": platform,
        "message": f"Input validated for {platform} platform",
    }


def prepare_workspace(data: dict[str, Any]) -> dict[str, Any]:
    """Prepare a workspace directory for the pipeline run."""
    run_id = data.get("run_id", datetime.now().strftime("%Y%m%d_%H%M%S"))
    base_dir = Path(CONFIG["WORK_BASE_DIR"])
    work_dir = base_dir / run_id

    # Create directory structure
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "input").mkdir(exist_ok=True)
    (work_dir / "output").mkdir(exist_ok=True)
    (work_dir / "logs").mkdir(exist_ok=True)
    (work_dir / "work").mkdir(exist_ok=True)

    return {
        "run_id": run_id,
        "work_dir": str(work_dir),
        "input_path": str(work_dir / "input"),
        "output_path": str(work_dir / "output"),
        "primer_bed_path": str(work_dir / "primers.bed"),
        "refseq_path": str(work_dir / "reference.fasta"),
        "ref_gbk_path": str(work_dir / "reference.gbk"),
    }


def run_pipeline(data: dict[str, Any]) -> dict[str, Any]:
    """Execute the OneRoof pipeline."""
    platform = data["platform"]
    work_dir = Path(data["work_dir"])

    # Build nextflow command
    cmd = [CONFIG["NEXTFLOW_PATH"], "run", CONFIG["ONEROOF_PATH"]]

    # Add platform-specific parameters
    if platform == "nanopore":
        cmd.extend(["--pod5_dir", data["input_path"]])
    else:  # illumina
        cmd.extend(["--illumina_fastq_dir", data["input_path"]])

    # Add common parameters
    cmd.extend(
        [
            "--primer_bed",
            data["primer_bed_path"],
            "--refseq",
            data["refseq_path"],
            "--ref_gbk",
            data["ref_gbk_path"],
            "--outdir",
            str(work_dir / "output"),
            "-work-dir",
            str(work_dir / "work"),
        ],
    )

    # Add profile
    profile = data.get("profile", "docker")
    cmd.extend(["-profile", profile])

    # Add any additional nextflow parameters
    nf_params = data.get("nextflow_params", {})
    for key, value in nf_params.items():
        cmd.extend([f"--{key}", str(value)])

    # Set up environment
    env = os.environ.copy()
    if CONFIG["USE_PIXI"]:
        # Run command through pixi
        cmd = [CONFIG["PIXI_PATH"], "run", "-e", "oneroof"] + cmd

    # Execute pipeline
    log_file = work_dir / "logs" / "nextflow.log"
    try:
        with open(log_file, "w") as log:
            process = subprocess.run(
                cmd,
                cwd=CONFIG["ONEROOF_PATH"],
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )

        if process.returncode == 0:
            return {
                "status": "success",
                "output_path": str(work_dir / "output"),
                "log_file": str(log_file),
                "message": "Pipeline completed successfully",
            }
        return {
            "status": "failed",
            "output_path": str(work_dir / "output"),
            "log_file": str(log_file),
            "message": f"Pipeline failed with exit code {process.returncode}",
        }

    except Exception as e:
        logger.error(f"Pipeline execution failed: {e!s}")
        return {
            "status": "failed",
            "error": str(e),
            "log_file": str(log_file),
        }


def collect_logs(data: dict[str, Any]) -> dict[str, Any]:
    """Collect logs and debugging information from a failed run."""
    work_dir = Path(data["work_dir"])
    log_dir = work_dir / "logs"

    logs = {}

    # Collect nextflow log
    nf_log = log_dir / "nextflow.log"
    if nf_log.exists():
        with open(nf_log) as f:
            logs["nextflow"] = f.read()[-10000:]  # Last 10KB

    # Collect .nextflow.log from work directory
    nf_hidden_log = work_dir / "work" / ".nextflow.log"
    if nf_hidden_log.exists():
        with open(nf_hidden_log) as f:
            logs["nextflow_hidden"] = f.read()[-10000:]

    # Create summary
    message = "Pipeline execution failed.\n"
    if "nextflow" in logs:
        # Try to extract error from log
        lines = logs["nextflow"].split("\n")
        error_lines = [line for line in lines if "ERROR" in line or "WARN" in line]
        if error_lines:
            message += "\nRecent errors:\n" + "\n".join(error_lines[-10:])

    return {
        "message": message,
        "logs": logs,
        "run_id": data.get("run_id"),
        "work_dir": str(work_dir),
    }


def cleanup_workspace(data: dict[str, Any]) -> dict[str, Any]:
    """Clean up workspace after pipeline run."""
    work_dir = Path(data["work_dir"])
    keep_logs = data.get("keep_logs", False)

    if not work_dir.exists():
        return {"status": "success", "message": "Workspace already cleaned"}

    try:
        if keep_logs:
            # Keep logs directory, remove everything else
            for item in work_dir.iterdir():
                if item.name != "logs":
                    if item.is_dir():
                        shutil.rmtree(item)
                    else:
                        item.unlink()
        else:
            # Remove entire workspace
            shutil.rmtree(work_dir)

        return {
            "status": "success",
            "message": f"Workspace cleaned: {work_dir}",
        }

    except Exception as e:
        logger.error(f"Cleanup failed: {e!s}")
        return {
            "status": "failed",
            "error": str(e),
        }


def check_status(data: dict[str, Any]) -> dict[str, Any]:
    """Check the status of a pipeline run."""
    work_dir = Path(data.get("work_dir", ""))

    if not work_dir.exists():
        return {"status": "not_found", "message": "Run directory not found"}

    # Check for output files
    output_dir = work_dir / "output"
    if output_dir.exists() and any(output_dir.iterdir()):
        return {"status": "completed", "output_files": len(list(output_dir.iterdir()))}

    # Check if still running
    work_files = work_dir / "work"
    if work_files.exists() and any(work_files.glob("*/.command.begin")):
        return {"status": "running"}

    return {"status": "unknown"}


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "false").lower() == "true"
    app.run(host="0.0.0.0", port=port, debug=debug)
