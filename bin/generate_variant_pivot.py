#!/usr/bin/env python3


import argparse
from pathlib import Path

from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
    TODO
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_table",
        "-i",
        type=Path,
        required=True,
        help="Input table of VCF fields extracted from an annotated VCF with SnpSift",
    )

    return parser.parse_args()


def main() -> None:
    """
    TODO
    """
    _ = parse_command_line_args()
    logger.info("Hi mom!")


if __name__ == "__main__":
    main()
