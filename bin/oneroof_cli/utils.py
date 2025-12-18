"""
Utility functions for the OneRoof CLI.

Provides console output helpers, Nextflow command building, and shared constants.
"""

import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

from rich.console import Console

# Shared console instances
console = Console()
err_console = Console(stderr=True)

# Path to the pipeline root (parent of bin/)
PIPELINE_ROOT = Path(__file__).parent.parent.parent

# Resume file location
RESUME_FILE = Path(".nfresume")


# =============================================================================
# Console Output Helpers
# =============================================================================


def error(message: str, exit_code: int = 1) -> None:
    """Print an error message and optionally exit."""
    err_console.print(f"[bold red]Error:[/bold red] {message}")
    if exit_code:
        sys.exit(exit_code)


def success(message: str) -> None:
    """Print a success message."""
    console.print(f"[bold green]Success:[/bold green] {message}")


def info(message: str) -> None:
    """Print an informational message."""
    console.print(f"[cyan]Info:[/cyan] {message}")


def warning(message: str) -> None:
    """Print a warning message."""
    console.print(f"[yellow]Warning:[/yellow] {message}")


# =============================================================================
# Nextflow Command Building
# =============================================================================


def generate_nextflow_command(args: dict[str, Any]) -> str:
    """
    Generate a Nextflow command string from a dictionary of arguments.

    Args:
        args: Dictionary of argument names to values. Keys should match
              Nextflow parameter names. None values are skipped.

    Returns:
        A complete Nextflow command string ready for execution.

    Notes:
        - Boolean True values become flags (--flag)
        - Boolean False values are skipped
        - All other values are quoted appropriately
        - The 'profile' key is handled specially as a Nextflow option (-profile)
    """
    nextflow_params: list[str] = []

    for arg, value in args.items():
        if value is None:
            continue

        # Handle profile specially - it's a Nextflow option, not a pipeline param
        if arg == "profile":
            # profile can be a list of strings
            if isinstance(value, list):
                profile_str = ",".join(value)
            else:
                profile_str = str(value)
            nextflow_params.append(f"-profile {shlex.quote(profile_str)}")
            continue

        if isinstance(value, bool):
            if value:
                nextflow_params.append(f"--{arg}")
            # False booleans are simply omitted
        else:
            nextflow_params.append(f"--{arg} {shlex.quote(str(value))}")

    params_str = " ".join(nextflow_params)
    return f"nextflow run {PIPELINE_ROOT} {params_str}"


# =============================================================================
# Nextflow Execution
# =============================================================================


def run_nextflow(command: str) -> None:
    """
    Execute a Nextflow command and save it for potential resume.

    Args:
        command: The complete Nextflow command string to execute.

    Notes:
        - Appends -resume to the command if not already present
        - Saves the resume-enabled command to .nfresume for later use
        - Executes the command via subprocess
    """
    if command.endswith("-resume"):
        resume_command = command
    else:
        resume_command = f"{command} -resume"

    # Save for future resume
    RESUME_FILE.write_text(resume_command, encoding="utf-8")

    # Execute
    split_command = shlex.split(resume_command)
    subprocess.run(split_command, check=True)  # noqa: S603


def resume_nextflow() -> None:
    """
    Resume a previous Nextflow run from the cached .nfresume file.

    Raises:
        SystemExit: If no previous run is detected (.nfresume doesn't exist).
    """
    if not RESUME_FILE.exists():
        error(
            "No previous run detected. "
            "Make sure you start with 'oneroof run' before using 'oneroof resume'."
        )

    command = RESUME_FILE.read_text(encoding="utf-8").strip()
    info(f"Resuming with command: {command}")
    run_nextflow(command)
