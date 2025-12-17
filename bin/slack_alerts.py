#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError


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
    if passing_message == "":
        passing_message = "No passing samples."
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
    if failing_message == "":
        failing_message = "No failing samples."
    return failing_message


def get_user_ids() -> list[str]:
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

    user_id_list = get_user_ids()
    slack_token = get_slack_token()

    block = [
        {
            "type": "section",
            "text": {
                "type": "plain_text",
                "text": f"Oneroof has finished successfully for experiment {run_label}, with {count_passing} samples passing. Below is a breakdown of which samples had greater than or equal to {coverage_threshold}X coverage.",
                "emoji": True,
            },
        },
        {
            "type": "section",
            "fields": [
                {
                    "type": "mrkdwn",
                    "text": f"*PASSING* \n ----------- \n {passing}  ",
                },
                {
                    "type": "mrkdwn",
                    "text": f"*FAILING* \n ----------- \n {failing}",
                },
            ],
        },
    ]

    for user_id in user_id_list:
        try:
            client = WebClient(token=slack_token)
            # Send the message to a channel
            client.chat_postMessage(
                channel=user_id,  # or a user ID like "U12345678"
                blocks=block,
            )
            print("Message sent successfully!")
        except SlackApiError as e:
            print(f"Error sending message: {e.response['error']}")


def main() -> None:
    args = parse_command_line_args()
    # get from launch dir
    run_label = args.run_label
    coverage_tsv_dir_path = args.input_tsv_dir
    coverage_depth = args.depth

    send_slack_notification(
        run_label,
        coverage_tsv_dir_path,
        coverage_depth,
    )


if __name__ == "__main__":
    main()
