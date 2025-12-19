"""
Typer application instance for the OneRoof CLI.

This module defines the main Typer app and any shared configuration.
Commands are registered via the commands subpackage.
"""

import typer

# The main Typer application instance
app = typer.Typer(
    name="oneroof",
    help="OneRoof: Base-, Variant-, and Consensus-calling under One Proverbial Roof.",
    rich_markup_mode="rich",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
