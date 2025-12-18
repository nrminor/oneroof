"""
The 'env' command for the OneRoof CLI.

Checks that all dependencies are available in the environment.
Currently a stub - implementation coming soon.
"""

import sys

from oneroof_cli.app import app
from oneroof_cli.utils import info


@app.command("env")
def check_environment() -> None:
    """
    [bold cyan]Check[/bold cyan] that all dependencies are available.

    Validates that required tools (Nextflow, etc.) are installed and
    accessible in the current environment.

    [dim]Status:[/dim] This command is not yet fully implemented.
    """
    info("The 'env' subcommand is not yet implemented but will be soon!")
    sys.exit(0)
