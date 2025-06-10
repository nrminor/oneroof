#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
import requests


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
        threshold = 1 - (1 / coverage_threshold)
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
        threshold = 1 - (1 / coverage_threshold)
        if proportion < threshold:
            failing_line = f"{sample_id}: {proportion}\n"
            failing_message += failing_line
    return failing_message


def get_user_ids() -> list[str]:
    """Resolve list of webhook URLs from env or default file."""
    # Check ONEROOF_SLACK_HOOKS environment variable
    path_str = os.environ.get("ONEROOF_SLACK_USER_IDS")

    path = (
        Path.home() / ".oneroof" / "slack.user_ids" if not path_str else Path(path_str)
    )

    # If the file doesn't exist or is empty, return empty list
    if not path.exists() or not path.is_file():
        return []

    # Read non-empty, non-comment lines
    with path.open() as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def get_slack_token() -> str:
    path_str = os.environ.get("ONEROOF_SLACK_TOKEN")

    path = Path.home() / ".oneroof" / "slack.token" if not path_str else Path(path_str)

    if not path.exists() or not path.is_file():
        return ""

    # Return the first non-empty line, stripped
    with path.open() as f:
        for line in f:
            if line.strip():
                return line.strip()
    return ""  # In case the file is empty


def send_slack_notification(
    run_label: str,
    stats_tsv: Path | str,
    coverage_threshold: int,
) -> None:
    stats_df = pd.read_csv(stats_tsv, sep="\t")

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

    user_id_list = get_user_ids()
    slack_token = get_slack_token()

    for user_id in user_id_list:
        resp = requests.post(
            "https://slack.com/api/conversations.open",
            headers={"Authorization": f"Bearer {slack_token}"},
            json={"users": user_id},
        )
        channel_id = resp.json().get("channel", {}).get("id")
        if not channel_id:
            raise RuntimeError("Failed to open conversation.")

        # Send the message
        msg_resp = requests.post(
            "https://slack.com/api/chat.postMessage",
            headers={"Authorization": f"Bearer {slack_token}"},
            json={"channel": channel_id, "text": complete_message},
        )


def main() -> None:
    args = parse_command_line_args()
    # get from launch dir
    run_label = args.run_label
    coverage_tsv_dir_path = args.input_tsv_dir
    coverage_depth = args.depth

    send_slack_notification(run_label, coverage_tsv_dir_path, coverage_depth)


if __name__ == "__main__":
    main()
