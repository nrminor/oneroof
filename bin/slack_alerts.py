import argparse
from pathlib import Path

import pandas as pd
import requests

SLACK_WEBHOOK_URL = (
    "https://hooks.slack.com/services/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"
)


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
        type=int,
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
        if proportion >= coverage_threshold:
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
        if proportion < coverage_threshold:
            failing_line = f"{sample_id}: {proportion}\n"
            failing_message += failing_line
    return failing_message


def send_slack_notification(exp_num, stats_tsv, coverage_threshold):
    # the webhook url
    SLACK_WEBHOOK_URL = ""
    # reading the tsv
    df = pd.read_csv(stats_tsv, sep="\t")

    # finding passing and failing
    passing, count_passing = passing_samples(df, coverage_threshold)
    failing = failing_samples(df, coverage_threshold)

    # creating the return message
    message = f"""
    Oneroof has finished successfully for experiment {exp_num}, with {count_passing} samples passing. Below is a breakdown of which samples had greater than or equal to 10X coverage.

    PASSING
    -------
    {passing}

    FAILING
    -------
    {failing}
    """

    payload = {"text": message}
    r = requests.post(SLACK_WEBHOOK_URL, json=payload)
    if (r.status_code) != 200:
        raise Exception(
            f"Error sending slack automation, response code: {r.status_code}"
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
