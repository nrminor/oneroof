"""
The 'validate' command for the OneRoof CLI.

Validates provided inputs before running the pipeline.
Currently a stub - implementation coming soon.
"""

import sys

from oneroof_cli.app import app
from oneroof_cli.utils import info


@app.command("validate")
def validate_inputs() -> None:
    """
    [bold yellow]Validate[/bold yellow] provided inputs.

    Checks that input files exist, are readable, and are in the expected
    format before running the pipeline.

    [dim]Status:[/dim] This command is not yet fully implemented.
    """
    info("The 'validate' subcommand is not yet implemented but will be soon!")
    sys.exit(0)
