#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
import requests
from loguru import logger


def parse_command_line_args() -> argparse.Namespace:
    """
    Strictly parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_tsv_dir",
        type=Path,
        required=True,
        help="Depth of coverage TSV file",
    )
    parser.add_argument(
        # depth of 20 for default
        "--depth",
        "-d",
        type=int,
        required=True,
        help="minimum depth of coverage",
    )
    # TODO(@akalweit5): Replace with run label str
    parser.add_argument(
        "--run_label",
        "-r",
        type=str,
        required=True,
        help="experiment number",
    )

    return parser.parse_args()


def passing_samples(df: pd.DataFrame, coverage_threshold) -> tuple[str, int]:
    # seeing what samples are above or equal to the coverage threshhold, also returning a count to find total number of samples passing
    passing_message = ""
    count = 0
    for i in range(len(df.iloc[:, 1])):
        proportion = df.iloc[i, 1]
        sample_id = df.iloc[i, 0]
        threshold = 1 / coverage_threshold
        if proportion >= threshold:
            count += 1
            passing_line = f"{sample_id}: {proportion}\n"
            passing_message += passing_line
    return (passing_message, count)


def failing_samples(df: pd.DataFrame, coverage_threshold) -> str:
    # seeing what samples are below the coverage threshhold
    failing_message = ""
    for i in range(len(df.iloc[:, 1])):
        proportion = df.iloc[i, 1]
        sample_id = df.iloc[i, 0]
        threshold = 1 / coverage_threshold
        if proportion < threshold:
            failing_line = f"{sample_id}: {proportion}\n"
            failing_message += failing_line
    return failing_message


def get_webhook_paths() -> list[str]:
    """Resolve list of webhook URLs from env or default file."""
    # Check ONEROOF_SLACK_HOOKS environment variable
    path_str = os.environ.get("ONEROOF_SLACK_HOOKS")

    path = (
        Path.home() / ".oneroof" / "slack.webhooks" if not path_str else Path(path_str)
    )

    # If the file doesn't exist or is empty, return empty list
    if not path.exists() or not path.is_file():
        return []

    # Read non-empty, non-comment lines
    with path.open() as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def send_slack_notification(
    run_label: str,
    stats_tsv: Path | str,
    coverage_threshold: int,
) -> None:
    # the webhook url
    # reading the tsv
    # print(stats_tsv)
    stats_df = pd.read_csv(stats_tsv, sep="\t")

    # getting the webhooks
    webhook_urls = get_webhook_paths()

    if not webhook_urls:
        logger.error("No webhook URLs found. Exiting.")
        return

    # finding passing and failing
    passing, count_passing = passing_samples(stats_df, coverage_threshold)
    failing = failing_samples(stats_df, coverage_threshold)

    # creating the return message
    message = (
        f"Oneroof has finished successfully for experiment {run_label}, "
        f"with {count_passing} samples passing. Below is a breakdown of "
        f"which samples had greater than or equal to {coverage_threshold}X coverage."
    )

    results = f"PASSING\n-------\n{passing}\n\nFAILING\n-------\n{failing}"

    complete_message = f"{message}\n```{results}```"

    payload = {"text": complete_message}
    for slack_webhook_url in webhook_urls:
        # TODO(@akalweit5): Add reasonable timeout and consider retry strategy
        r = requests.post(slack_webhook_url, json=payload)
        if (r.status_code) != 200:  # noqa: PLR2004
            msg = f"Error sending slack automation, response code: {r.status_code}"
            # TODO(@akalweit5): Find a better exception here. What actually is the error we're expecting? And could we just use
            # an assert somewhere to crash early?
            raise Exception(msg)


def main() -> None:
    args = parse_command_line_args()
    # get from launch dir
    run_label = args.run_label
    coverage_tsv_dir_path = args.input_tsv_dir
    coverage_depth = args.depth

    send_slack_notification(run_label, coverage_tsv_dir_path, coverage_depth)


if __name__ == "__main__":
    main()
