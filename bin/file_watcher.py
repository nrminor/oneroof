#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "loguru",
#     "paramiko",
#     "pyyaml",
# ]
# ///
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import yaml
from loguru import logger
import paramiko
from paramiko.client import SSHClient

# Default timing constants (in seconds unless noted)
DEFAULT_WATCH_INTERVAL_SECONDS = 30
DEFAULT_FILE_READY_CHECK_INTERVAL = 3
DEFAULT_FILE_READY_MAX_CHECKS = 10
DEFAULT_TRANSFER_RETRY_WAIT = 3
DEFAULT_TRANSFER_MAX_RETRIES = 10
DEFAULT_VERIFY_RETRY_WAIT = 3
DEFAULT_VERIFY_MAX_RETRIES = 10
SECONDS_PER_HOUR = 3600


@dataclass(frozen=True)
class Credentials:
    watch_path: str
    pattern: str
    host: str
    username: str
    password: str
    watch_duration: int


class TransferRunner:
    """
    Handles file transfers from a remote server via SFTP.

    Attributes:
        client: SSH client connection to the remote server
        filename: Name of the file to transfer
        remote_path: Directory path on the remote server
        address: Remote server address
        username: Username for the remote connection
        local_path: Local directory to save transferred files
        ready: Whether the file is ready for transfer
        remote_hash: SHA256 hash of the remote file (set after prepare())
    """

    def __init__(
        self,
        client: SSHClient,
        filename: str,
        remote_path: str,
        address: str,
        username: str,
        local_path: str | Path = Path.cwd(),
    ) -> None:
        self.client = client
        self.filename = filename
        self.remote_path = remote_path
        self.address = address
        self.username = username
        self.local_path = local_path
        self.ready = False
        self.remote_hash: str | None = None
        self._temp_path: Path | None = None

    def prepare(self) -> bool:
        """
        Prepare for file transfer by checking readiness and computing remote hash.

        Checks if:
        1. The remote file is ready (size has stabilized)
        2. The file doesn't already exist locally

        If both conditions are met, computes the remote file hash for later
        verification.

        Returns:
            True if ready to transfer, False otherwise
        """
        # Check if file already exists locally
        if self.filename in os.listdir(self.local_path):
            logger.debug(f"File {self.filename} already exists locally, skipping")
            self.ready = False
            return False

        # Check if remote file is ready
        if not self.is_remote_file_ready():
            logger.debug(f"Remote file {self.filename} is not ready yet")
            self.ready = False
            return False

        # Compute remote hash before transfer
        self.remote_hash = self._compute_remote_hash()
        self.ready = True
        return True

    def transfer_file(
        self,
        max_retries: int = DEFAULT_TRANSFER_MAX_RETRIES,
        wait_time: int = DEFAULT_TRANSFER_RETRY_WAIT,
        attempts: int = 0,
    ) -> None:
        """
        Transfers a file from the remote server to a temporary local file.
        Retries the transfer if an error occurs, up to a maximum number of retries.

        The file is downloaded to a .tmp file first. Call verify_transfer() to
        verify the hash and rename to the final filename.

        Args:
            max_retries (int, optional): Maximum number of retry attempts. Defaults to 10.
            wait_time (int, optional): Time in seconds to wait between retries. Defaults to 3.
            attempts (int, optional): Current retry attempt count. Defaults to 0.

        Returns:
            None
        """
        sftp = self.client.open_sftp()
        remote_file_path = f"{self.remote_path}/{self.filename}"
        self._temp_path = Path(self.local_path) / f"{self.filename}.tmp"
        attempts += 1

        try:
            sftp.get(remote_file_path, self._temp_path)
            logger.info(
                f"File {self.filename} successfully transferred to temp file.",
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(
                f"Error transferring file {self.filename}: {e}. {max_retries} retries remaining.",
            )
            # Clean up partial temp file if it exists
            if self._temp_path and self._temp_path.exists():
                self._temp_path.unlink()
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
        wait_time: int = DEFAULT_FILE_READY_CHECK_INTERVAL,
        max_checks: int = DEFAULT_FILE_READY_MAX_CHECKS,
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

        try:
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
        finally:
            sftp.close()

    def _compute_file_hash(self, file_path: str | Path, chunk_size: int = 65536) -> str:
        """
        Compute SHA256 hash of a local file by reading it in chunks.

        Args:
            file_path: Path to the local file
            chunk_size: Size of chunks to read (default 64KB)

        Returns:
            Hexadecimal SHA256 hash string
        """
        sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                sha256.update(chunk)
        return sha256.hexdigest()

    def _compute_remote_hash(self, chunk_size: int = 65536) -> str:
        """
        Compute SHA256 hash of the remote file by reading it via SFTP in chunks.

        Args:
            chunk_size: Size of chunks to read (default 64KB)

        Returns:
            Hexadecimal SHA256 hash string
        """
        sftp = self.client.open_sftp()
        try:
            remote_file_path = f"{self.remote_path}/{self.filename}"
            sha256 = hashlib.sha256()
            with sftp.open(remote_file_path, "rb") as f:
                for chunk in iter(lambda: f.read(chunk_size), b""):
                    sha256.update(chunk)
            return sha256.hexdigest()
        finally:
            sftp.close()

    def verify_transfer(self) -> bool:
        """
        Verify the transferred file by comparing temp file hash against remote hash.

        If verification succeeds, renames the temp file to the final filename.
        If verification fails, removes the temp file.

        Returns:
            True if hashes match and file was renamed, False otherwise
        """
        if not self._temp_path or not self._temp_path.exists():
            logger.error(f"Temp file for {self.filename} does not exist")
            return False

        local_hash = self._compute_file_hash(self._temp_path)
        transfer_check = local_hash == self.remote_hash

        if transfer_check:
            # Rename temp file to final filename
            final_path = Path(self.local_path) / self.filename
            self._temp_path.rename(final_path)
            logger.info(f"File {self.filename} verified and saved.")
            return True

        # Verification failed - clean up temp file
        logger.warning(
            f"Hash mismatch for {self.filename}: local={local_hash}, remote={self.remote_hash}"
        )
        self._temp_path.unlink()
        return False


def parse_command_line_args() -> argparse.Namespace:
    """
    Parses command-line arguments for configuring the file watcher.

    Returns:
        argparse.Namespace: Parsed command-line arguments including:
            - host_config (Path): YAML configuration file for connecting to the remote host.
    """
    parser = argparse.ArgumentParser(
        description="Watch a remote directory for files and transfer them locally."
    )
    parser.add_argument(
        "--host_config",
        "-c",
        type=Path,
        required=True,
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
    assert "watch_duration" in entries, """
    A watch duration in hours must be provided in the `watch_duration` field
    of the file watcher configuration YAML file.
    """
    watch_duration = config_dict["watch_duration"]
    assert isinstance(watch_duration, int) and watch_duration > 0, """
    The `watch_duration` field must be a positive integer (hours to watch).
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
    creds = find_credentials(args.host_config)

    # set up ssh client and connection
    client = SSHClient()
    client.load_system_host_keys()
    logger.warning(
        "Unknown host keys will be automatically added. "
        "Verify the host key fingerprint for production use."
    )
    client.set_missing_host_key_policy(paramiko.WarningPolicy())
    client.connect(creds.host, username=creds.username, password=creds.password)

    # Record the start time
    start_time = time.time()
    duration = creds.watch_duration * SECONDS_PER_HOUR

    # run the watcher until the watch duration
    try:
        while True:
            # Check if the elapsed time exceeds the configured duration
            if time.time() - start_time > duration:
                logger.info("Time limit reached. Stopping the file transfer process.")
                break

            _, stdout, _ = client.exec_command(f"ls {creds.watch_path}")
            all_files = stdout.read().decode(encoding="utf8").splitlines()
            file_queue = [
                file for file in all_files if fnmatch.fnmatch(file, creds.pattern)
            ]
            for file in file_queue:
                runner = TransferRunner(
                    client=client,
                    filename=file,
                    remote_path=creds.watch_path,
                    address=creds.host,
                    username=creds.username,
                )
                try:
                    if not runner.prepare():
                        continue
                    runner.transfer_file()
                    if not runner.verify_transfer():
                        logger.error(
                            f"The file, {runner.filename}, was corrupted during the transfer process.",
                        )
                except FileNotFoundError:
                    logger.warning(f"File {file} no longer exists on remote, skipping")
                except OSError as e:
                    logger.warning(f"Error processing file {file}: {e}, skipping")
            time.sleep(DEFAULT_WATCH_INTERVAL_SECONDS)
    except KeyboardInterrupt:
        logger.debug("Process interrupted. Exiting...")
        sys.exit(1)
    finally:
        # Close the SSH connection
        client.close()
        logger.info("SSH connection closed.")


if __name__ == "__main__":
    main()
