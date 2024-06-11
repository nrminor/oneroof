#!/usr/bin/env python3

import argparse
import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import paramiko
import yaml
from paramiko import SSHClient


@dataclass
class Credentials:
    watch_path: str
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
    local_path: str | Path = field(default=os.getcwd())

    def __post_init__(self):
        if self.is_file_done_writing():
            object.__setattr__(self, "ready", True)
        else:
            object.__setattr__(self, "ready", False)

        if self.ready and self.filename not in os.listdir(self.local_path):
            return self
        else:
            time.sleep(10)
            return TransferRunner(
                self.client,
                self.filename,
                self.remote_path,
                self.address,
                self.username,
            )

    def transfer_file(self):
        sftp = self.client.open_sftp()
        remote_file_path = f"{self.remote_path}/{self.filename}"
        local_file_path = os.path.join(self.local_path, self.filename)

        try:
            sftp.get(remote_file_path, local_file_path)
            print(
                f"File {self.filename} successfully transferred to '{self.local_path}'."
            )
        except Exception as e:
            print(f"Error transferring file {self.filename}: {e}")
        finally:
            sftp.close()

    def is_file_done_writing(
        self,
        wait_time: int = 3,
        max_checks: int = 10,
    ):
        sftp = self.client.open_sftp()
        previous_size = -1
        checks = 0

        while checks < max_checks:
            try:
                current_size = sftp.stat(
                    os.path.join(self.remote_path, self.filename)
                ).st_size
                if current_size == previous_size:
                    return True
                previous_size = current_size
                checks += 1
                time.sleep(wait_time)
            except IOError as e:
                # Handle file not found or other I/O errors
                print(f"Error checking file size: {e}")
                return False
        return False


def parse_command_line_args() -> argparse.Namespace:
    """
    TODO
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--watchpath",
        "-w",
        type=Path,
        required=False,
        help="Directory to watch on remote host.",
    )
    parser.add_argument(
        "--watch_duration",
        "-d",
        type=int,
        required=False,
        default=72,
        help="The number of hours to watch for new pod5 files.",
    )
    parser.add_argument(
        "--host_config",
        "-c",
        type=Path,
        required=True,
        default="pod5_watcher.yml",
        help="The configuration YAML file to use for connecting to the remote host with credentials.",
    )

    args = parser.parse_args()
    return args


def parse_credential_config(args: argparse.Namespace) -> Credentials:
    with open(args.config_path, "r", encoding="utf8") as config_handle:
        config_dict = yaml.safe_load(config_handle)

    return Credentials(
        watch_path=args.watch_path if args.watch_path else config_dict["watch_path"],
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
    client = paramiko.client.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(creds.host, username=creds.username, password=creds.password)

    # Record the start time
    start_time = time.time()
    duration = creds.watch_duration * 60 * 60

    # run the watcher until the watch duration, checking fo pod5 files every 30
    # seconds
    try:
        while True:
            # Check if the elapsed time exceeds 72 hours
            if time.time() - start_time > duration:
                print("Time limit reached. Stopping the file transfer process.")
                break

            _stdin, _stdout, _stderr = client.exec_command(f"ls {creds.watch_path}")
            all_files = _stdout.read().decode(encoding="utf8").splitlines()
            pod5_queue = [file for file in all_files if file.endswith(".pod5")]
            for pod5 in pod5_queue:
                runner = TransferRunner(
                    client=client,
                    filename=pod5,
                    remote_path=creds.watch_path,
                    address=creds.host,
                    username=creds.username,
                )
                runner.transfer_file()
            time.sleep(30)
    except KeyboardInterrupt:
        print("Process interrupted. Exiting...")
    finally:
        # Close the SSH connection
        client.close()
        print("SSH connection closed.")


if __name__ == "__main__":
    main()
