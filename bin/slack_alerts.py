#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import pandas as pd
import requests

# SLACK_WEBHOOK_URL = (
#     "https://hooks.slack.com/services/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"
# )

# list of all users webhook url's
SLACK_WEBHOOK_URLS = [
    "https://hooks.slack.com/services/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",
]


def parse_command_line_args() -> argparse.Namespace:
    """
    Strictly parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_tsv",
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
    parser.add_argument(
        # depth of 20 for default
        "--exp_num",
        "-e",
        type=Path,
        required=True,
        help="experiment number",
    )
    return parser.parse_args()


def passing_samples(df, coverage_threshold):
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


def failing_samples(df, coverage_threshold):
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

    if path_str is None:
        # Fallback to default path
        path = Path.home() / ".oneroof" / "slack.webhooks"
    else:
        path = Path(path_str)

    # If the file doesn't exist or is empty, return empty list
    if not path.exists() or not path.is_file():
        return []

    # Read non-empty, non-comment lines
    with path.open() as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def send_slack_notification(exp_num, stats_tsv, coverage_threshold):
    # the webhook url
    # reading the tsv
    df = pd.read_csv(stats_tsv, sep="\t")

    # getting the webhooks
    webhook_urls = get_webhook_paths()

    if not webhook_urls:
        print("No webhook URLs found. Exiting.")
        return

    # finding passing and failing
    passing, count_passing = passing_samples(df, coverage_threshold)
    failing = failing_samples(df, coverage_threshold)

    # getting exp num
    exp_number = str(exp_num).split("/")[-1]

    # creating the return message
    message = (
        f"Oneroof has finished successfully for experiment {exp_number}, "
        f"with {count_passing} samples passing. Below is a breakdown of "
        f"which samples had greater than or equal to {coverage_threshold}X coverage."
    )

    results = f"PASSING\n-------\n{passing}\n\nFAILING\n-------\n{failing}"

    complete_message = f"{message}\n```{results}```"

    payload = {"text": complete_message}
    for SLACK_WEBHOOK_URL in webhook_urls:
        r = requests.post(SLACK_WEBHOOK_URL, json=payload)
        if (r.status_code) != 200:
            raise Exception(
                f"Error sending slack automation, response code: {r.status_code}",
            )


def main():
    args = parse_command_line_args()
    # get from launch dir
    experiment_number = args.exp_num
    coverage_tsv = args.input_tsv
    coverage_depth = args.depth

    send_slack_notification(experiment_number, coverage_tsv, coverage_depth)


if __name__ == "__main__":
    main()
