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
import os
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

    def __post_init__(self) -> TransferRunner:
        if self.is_file_done_writing():
            object.__setattr__(self, "ready", True)
        else:
            object.__setattr__(self, "ready", False)

        if self.ready and self.filename not in os.listdir(self.local_path):
            return self

        time.sleep(10)
        return TransferRunner(
            self.client,
            self.filename,
            self.remote_path,
            self.address,
            self.username,
        )

    def transfer_file(self) -> None:
        sftp = self.client.open_sftp()
        remote_file_path = f"{self.remote_path}/{self.filename}"
        local_file_path = Path(self.local_path) / Path(self.filename)

        try:
            sftp.get(remote_file_path, local_file_path)
            logger.info(
                f"File {self.filename} successfully transferred to '{self.local_path}'.",
            )
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Error transferring file {self.filename}: {e}")
        finally:
            sftp.close()

    def is_file_done_writing(
        self,
        wait_time: int = 3,
        max_checks: int = 10,
    ) -> bool:
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


def parse_command_line_args() -> argparse.Namespace:
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


def parse_credential_config(args: argparse.Namespace) -> Credentials:
    with Path(args.config_path).open(encoding="utf8") as config_handle:
        config_dict = yaml.safe_load(config_handle)

    runtime_config_check(config_dict)

    config_dict["watch_path"] = (
        args.watch_path if args.watch_path else config_dict["watch_path"]
    )
    config_dict["pattern"] = (
        args.watch_pattern if args.watch_path else config_dict["pattern"]
    )

    return Credentials(
        watch_path=config_dict["watch_path"],
        pattern=config_dict["pattern"],
        host=config_dict["host"],
        username=config_dict["username"],
        password=config_dict["password"],
        watch_duration=args.watch_duration
        if args.watch_duration
        else config_dict["watch_duration"],
    )


def main() -> None:
    # parse any command line args
    args = parse_command_line_args()

    # parse the config file
    creds = parse_credential_config(args)

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
            time.sleep(30)
    except KeyboardInterrupt:
        logger.debug("Process interrupted. Exiting...")
    finally:
        # Close the SSH connection
        client.close()
        logger.info("SSH connection closed.")


if __name__ == "__main__":
    main()
