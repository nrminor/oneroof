import requests

SLACK_WEBHOOK_URL = "https://hooks.slack.com/services/YOUR/WEBHOOK/URL"


def slack_automated_message(experiment: str, samples_passing: int) -> None:
    message = f""" Oneroof has finished successfully for experiment {experiment}, with {samples_passing} samples passing. Below is a breakdown of which samples had greater than or equal to 10X coverage.

    PASSING
    -------
    sample1: 0.97
    sample4: 0.92

    FAILING
    -------
    sample2: 0.83
    sample4: 0.15'

    """

    payload = {"text": message}
    requests.post(SLACK_WEBHOOK_URL, json=payload)
