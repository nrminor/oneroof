#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "loguru",
#     "paramiko",
# ]
# ///
from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import yaml
from loguru import logger
from paramiko.client import SSHClient


@dataclass
class Credentials:
    watch_path: str
    pattern: str
    host: str
    username: str
    password: str
    watch_duration: int


@dataclass
class TransferRunner:
    client: SSHClient
    filename: str
    remote_path: str
    address: str
    username: str
    ready: bool = field(init=False)
    local_path: str | Path = field(default=Path.cwd())
    remote_hash: str = field(init=False)

    def __post_init__(self) -> TransferRunner:
        """
        Initializes the TransferRunner instance by checking if the remote file is ready.
        If the file is available and not already present in the local directory,
        the method returns a new instance of TransferRunner after a brief delay.

        Returns:
            TransferRunner: A new instance if the conditions are met; otherwise, returns self.
        """
        if self.is_remote_file_ready():
            object.__setattr__(self, "ready", True)
        else:
            object.__setattr__(self, "ready", False)

        if self.ready and self.filename not in os.listdir(self.local_path):
            remote_hash_bytes = str(self).encode("utf-8")
            self.remote_hash = hashlib.sha256(remote_hash_bytes).hexdigest()
            return self

        time.sleep(10)
        return TransferRunner(
            self.client,
            self.filename,
            self.remote_path,
            self.address,
            self.username,
        )

    def transfer_file(
        self,
        max_retries: int = 10,
        wait_time: int = 3,
        attempts: int = 0,
    ) -> None:
        """
        Transfers a file from the remote server to the local directory using SFTP.
        Retries the transfer if an error occurs, up to a maximum number of retries.

        Args:
            max_retries (int, optional): Maximum number of retry attempts. Defaults to 10.
            wait_time (int, optional): Time in seconds to wait between retries. Defaults to 3.
            attempts (int, optional): Current retry attempt count. Defaults to 0.

        Returns:
            None
        """
        sftp = self.client.open_sftp()
        remote_file_path = f"{self.remote_path}/{self.filename}"
        local_file_path = Path(self.local_path) / Path(self.filename)
        attempts += 1

        try:
            sftp.get(remote_file_path, local_file_path)
            logger.info(
                f"File {self.filename} successfully transferred to '{self.local_path}'.",
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(
                f"Error transferring file {self.filename}: {e}. {max_retries} retries remaining.",
            )
            time.sleep(wait_time)
            if attempts > max_retries:
                logger.error(
                    f"Transfer for file failed after {max_retries} retries. Skipping.",
                )
                return
            self.transfer_file(
                max_retries=max_retries,
                wait_time=wait_time,
                attempts=attempts,
            )

        finally:
            sftp.close()

    def is_remote_file_ready(
        self,
        wait_time: int = 3,
        max_checks: int = 10,
    ) -> bool:
        """
        Checks if the remote file is ready for transfer by monitoring its size.
        A file is considered ready if its size remains unchanged for a series of checks.

        Args:
            wait_time (int, optional): Time in seconds to wait between checks. Defaults to 3.
            max_checks (int, optional): Maximum number of checks to perform. Defaults to 10.

        Returns:
            bool: True if the file is ready, False otherwise.
        """
        sftp = self.client.open_sftp()
        previous_size = -1
        checks = 0

        while checks < max_checks:
            try:
                current_size = sftp.stat(
                    str(Path(self.remote_path) / Path(self.filename)),
                ).st_size
                if current_size == previous_size:
                    return True
                previous_size = current_size
                checks += 1
                time.sleep(wait_time)
            except OSError as e:  # noqa: PERF203
                # Handle file not found or other I/O errors
                logger.warning(f"Error checking file size: {e}")
                return False
        return False

    def verify_transfer(
        self,
        max_retries: int = 10,
        wait_time: int = 3,
        attempts: int = 0,
    ) -> bool:
        attempts += 1

        # encoding local file to bytes
        local_file_bytes = self.filename.encode("utf-8")

        # creating the hashes for the local file
        local_hash = hashlib.sha256(local_file_bytes).hexdigest()

        # checking if the two objects are the same
        transfer_check = local_hash == self.remote_hash

        if transfer_check:
            return transfer_check
        time.sleep(wait_time)
        if attempts > max_retries:
            return False
        return self.verify_transfer(attempts=attempts)


def parse_command_line_args() -> argparse.Namespace:
    """
    Parses command-line arguments for configuring the file watcher.

    Returns:
        argparse.Namespace: Parsed command-line arguments including:
            - watch_path (Path): Directory to watch on the remote host.
            - watch_pattern (str): Pattern for matching files on the remote host.
            - watch_duration (int): Number of hours to watch for new files.
            - host_config (Path): YAML configuration file for connecting to the remote host.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--watch_path",
        "-w",
        type=Path,
        required=False,
        default=None,
        help="Directory to watch on remote host.",
    )
    parser.add_argument(
        "--watch_pattern",
        "-p",
        type=str,
        required=False,
        default=None,
        help="Pattern to use when watching for files on the remote host.",
    )
    parser.add_argument(
        "--watch_duration",
        "-d",
        type=int,
        required=False,
        default=72,
        help="The number of hours to watch for new files.",
    )
    parser.add_argument(
        "--host_config",
        "-c",
        type=Path,
        required=True,
        default="file_watcher.yml",
        help="The configuration YAML file to use for connecting to the remote host with credentials.",
    )
    return parser.parse_args()


def runtime_config_check(config_dict: dict) -> None:
    """
    Validates the runtime configuration dictionary to ensure all required fields are present.

    Args:
        config_dict (dict): The configuration dictionary parsed from the YAML file.

    Raises:
        AssertionError: If any required field is missing or if the file pattern is unsupported.
    """
    entries = config_dict.keys()
    assert "watch_path" in entries, """
    A remote absolute file path, not including the IP address itself, must be
    provided in a `watch_path` field in the configuration YAML file.
    """
    assert "pattern" in entries, """
    A valid glob pattern to use for matching files, e.g., "*.pod5", must be
    provided in the `pattern` field of the configuration YAML file.
    """
    supported_patterns = ["*.fastq.gz", "*.bam", "*.pod5"]
    assert config_dict["pattern"] in supported_patterns, f"""
    The file watcher currently only supports searching for the following patterns:
    {supported_patterns}
    """
    assert "host" in entries, """
    A valid host IP address must be provided in the `host` field of the file
    watcher configuration YAML file.
    """
    assert "username" in entries, """
    A valid username at the remote host must be provided in the `username` field
    of the file watch configuration YAML file.
    """
    assert "password" in entries, """
    A valid password at the remote host must be provided in the `password` field
    of the file watch configuration YAML file.
    """


def make_credential_error(
    env_var: str,
    config_field: str,
    config_path: str | Path,
) -> str:
    """
    Generates an error message for missing environment variables or configuration fields.

    Args:
        env_var (str): The name of the required environment variable.
        config_field (str): The corresponding field in the configuration file.
        config_path (str | Path): The path to the configuration file.

    Returns:
        str: The formatted error message.
    """
    return f"The ${env_var} environment variable is unset and the '{config_field}' field in the provided config file, '{config_path}', is missing."


def try_access_env_setting(
    env_var: str,
    config_field: str,
    config_path: str | Path,
) -> str:
    """
    Attempts to access a required setting from an environment variable.

    Args:
        env_var (str): The name of the required environment variable.
        config_field (str): The corresponding field in the configuration file.
        config_path (str | Path): The path to the configuration file.

    Returns:
        str: The retrieved setting value.

    Raises:
        OSError: If the environment variable is not set or empty.
    """
    message = make_credential_error(
        env_var,
        config_field,
        config_path,
    )
    if os.getenv(env_var):
        setting = os.getenv(env_var)
        if not setting:
            raise OSError(message)
    else:
        raise OSError(message)

    return setting


def try_access_config(
    env_var: str,
    config_field: str,
    config_path: str | Path,
    config_dict: dict[str, str],
) -> str:
    """
    Attempts to retrieve a required setting from the configuration dictionary.

    Args:
        env_var (str): The corresponding environment variable name (for logging purposes).
        config_field (str): The field name in the configuration file.
        config_path (str | Path): The path to the configuration file.
        config_dict (dict[str, str]): The parsed configuration dictionary.

    Returns:
        str: The retrieved setting value.

    Raises:
        OSError: If the required setting is missing in the configuration file.
    """
    message = make_credential_error(
        env_var,
        config_field,
        config_path,
    )
    setting = config_dict.get(config_field)
    if not setting:
        raise OSError(message)

    return setting


def try_access_setting(
    env_var: str,
    config_field: str,
    config_path: str | Path,
    config_dict: dict[str, str],
) -> str:
    """
    Retrieves a required setting from either an environment variable or the configuration file.

    Args:
        env_var (str): The corresponding environment variable name.
        config_field (str): The field name in the configuration file.
        config_path (str | Path): The path to the configuration file.
        config_dict (dict[str, str]): The parsed configuration dictionary.

    Returns:
        str: The retrieved setting value.

    Raises:
        AssertionError: If the provided configuration file does not exist.
    """
    if os.getenv(env_var):
        setting = try_access_env_setting(env_var, config_field, config_path)
    else:
        assert Path(config_path).exists(), (
            f"The provided config path {config_path} does not point to a file that exists."
        )
        setting = try_access_config(env_var, config_field, config_path, config_dict)

    return setting


def find_credentials(config_path: str | Path) -> Credentials:
    """
    Parses the configuration file and retrieves credentials required for file watching.

    Args:
        config_path (str | Path): The path to the configuration file.

    Returns:
        Credentials: A Credentials object containing the extracted configuration values.

    Raises:
        OSError: If any required field is missing from the configuration or environment variables.
    """
    config_check = Path(config_path).exists()
    config_dict: dict[str, str] = {}
    if config_check:
        with Path(config_path).open(encoding="utf8") as config_handle:
            config_dict = yaml.safe_load(config_handle)
        runtime_config_check(config_dict)

    watch_path = try_access_setting(
        "ONEROOF_WATCH_PATH",
        "watch_path",
        config_path,
        config_dict,
    )

    pattern = try_access_setting(
        "ONEROOF_WATCH_PATTERN",
        "pattern",
        config_path,
        config_dict,
    )

    host = try_access_setting(
        "ONEROOF_WATCH_HOST",
        "host",
        config_path,
        config_dict,
    )

    username = try_access_setting(
        "ONEROOF_WATCH_USERNAME",
        "username",
        config_path,
        config_dict,
    )

    password = try_access_setting(
        "ONEROOF_WATCH_PASSWORD",
        "password",
        config_path,
        config_dict,
    )

    watch_duration = int(
        try_access_setting(
            "ONEROOF_WATCH_DURATION",
            "watch_duration",
            config_path,
            config_dict,
        ),
    )

    return Credentials(
        watch_path,
        pattern,
        host,
        username,
        password,
        watch_duration,
    )


def main() -> None:
    # parse any command line args
    args = parse_command_line_args()

    # parse the config file
    creds = find_credentials(args.config_path)

    # set up ssh client and connection
    client = SSHClient()
    # This may or may not be necessary: client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(creds.host, username=creds.username, password=creds.password)

    # Record the start time
    start_time = time.time()
    duration = creds.watch_duration * 60 * 60

    # run the watcher until the watch duration, checking for files every 30
    # seconds
    try:
        while True:
            # Check if the elapsed time exceeds 72 hours
            if time.time() - start_time > duration:
                logger.info("Time limit reached. Stopping the file transfer process.")
                break

            _, stdout, _ = client.exec_command(f"ls {creds.watch_path}")
            all_files = stdout.read().decode(encoding="utf8").splitlines()
            file_queue = [file for file in all_files if file.endswith(creds.pattern)]
            for file in file_queue:
                runner = TransferRunner(
                    client=client,
                    filename=file,
                    remote_path=creds.watch_path,
                    address=creds.host,
                    username=creds.username,
                )
                runner.transfer_file()
                # TODO: check local hash against remote hash here
                if not runner.verify_transfer():
                    logger.error(
                        f"The file, {runner.filename}, was corrupted during the transfer process.",
                    )
            time.sleep(30)
    except KeyboardInterrupt:
        logger.debug("Process interrupted. Exiting...")
        sys.exit(1)
    finally:
        # Close the SSH connection
        client.close()
        logger.info("SSH connection closed.")


if __name__ == "__main__":
    main()
