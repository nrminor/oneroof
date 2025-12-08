#!/usr/bin/env python3
"""
Comprehensive pytest test module for file_watcher.py
Tests all major components including file watching, remote connections,
configuration parsing, file transfers, and error handling.
"""

from __future__ import annotations

import hashlib
import os
import sys
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import Mock, mock_open, patch

import pytest
import yaml
from paramiko.client import SSHClient
from paramiko.sftp_attr import SFTPAttributes

# Import the modules we're testing
from file_watcher import (
    Credentials,
    TransferRunner,
    find_credentials,
    make_credential_error,
    parse_command_line_args,
    runtime_config_check,
    try_access_config,
    try_access_env_setting,
    try_access_setting,
)


# Fixtures
@pytest.fixture
def mock_ssh_client():
    """Mock SSH client for testing"""
    client = Mock(spec=SSHClient)
    sftp = Mock()
    client.open_sftp.return_value = sftp
    return client


@pytest.fixture
def mock_credentials():
    """Mock credentials for testing"""
    return Credentials(
        watch_path="/remote/path",
        pattern="*.pod5",
        host="192.168.1.100",
        username="testuser",
        password="testpass",
        watch_duration=72,
    )


@pytest.fixture
def valid_config_dict():
    """Valid configuration dictionary"""
    return {
        "watch_path": "/remote/path",
        "pattern": "*.pod5",
        "host": "192.168.1.100",
        "username": "testuser",
        "password": "testpass",
        "watch_duration": "72",
    }


@pytest.fixture
def mock_sftp_attr():
    """Mock SFTP file attributes"""
    attr = Mock(spec=SFTPAttributes)
    attr.st_size = 1024
    return attr


class TestCredentials:
    """Test the Credentials dataclass"""

    def test_credentials_creation(self):
        """Test creating a Credentials instance"""
        creds = Credentials(
            watch_path="/path",
            pattern="*.pod5",
            host="localhost",
            username="user",
            password="pass",
            watch_duration=24,
        )
        assert creds.watch_path == "/path"
        assert creds.pattern == "*.pod5"
        assert creds.host == "localhost"
        assert creds.username == "user"
        assert creds.password == "pass"
        assert creds.watch_duration == 24

    def test_credentials_frozen(self):
        """Test that Credentials is frozen (immutable)"""
        creds = Credentials(
            watch_path="/path",
            pattern="*.pod5",
            host="localhost",
            username="user",
            password="pass",
            watch_duration=24,
        )
        with pytest.raises(FrozenInstanceError):
            creds.watch_path = "/new/path"


class TestTransferRunner:
    """Test the TransferRunner class"""

    @patch("file_watcher.time.sleep")
    @patch("os.listdir")
    def test_transfer_runner_init_file_ready(
        self, mock_listdir, mock_sleep, mock_ssh_client, mock_sftp_attr
    ):
        """Test TransferRunner initialization when file is ready"""
        mock_listdir.return_value = []  # File not in local directory
        mock_ssh_client.open_sftp.return_value.stat.return_value = mock_sftp_attr

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        assert runner.ready is True
        assert hasattr(runner, "remote_hash")
        assert isinstance(runner.remote_hash, str)

    @patch("file_watcher.time.sleep")
    @patch("os.listdir")
    def test_transfer_runner_init_file_not_ready(
        self, mock_listdir, mock_sleep, mock_ssh_client
    ):
        """Test TransferRunner initialization when file is not ready"""
        mock_listdir.return_value = []
        mock_ssh_client.open_sftp.return_value.stat.side_effect = [
            Mock(st_size=100),
            Mock(st_size=200),  # Size changing
        ]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        assert runner.ready is False
        mock_sleep.assert_called_once_with(10)

    @patch("file_watcher.time.sleep")
    def test_transfer_file_success(self, mock_sleep, mock_ssh_client):
        """Test successful file transfer"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )
        runner.ready = True

        runner.transfer_file()

        sftp_mock.get.assert_called_once()
        sftp_mock.close.assert_called_once()

    @patch("file_watcher.time.sleep")
    def test_transfer_file_retry_on_error(self, mock_sleep, mock_ssh_client):
        """Test file transfer retry mechanism"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        sftp_mock.get.side_effect = [Exception("Connection error"), None]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )
        runner.ready = True

        runner.transfer_file(max_retries=2, wait_time=1)

        assert sftp_mock.get.call_count == 2
        mock_sleep.assert_called_with(1)

    @patch("file_watcher.time.sleep")
    def test_transfer_file_max_retries_exceeded(self, mock_sleep, mock_ssh_client):
        """Test file transfer when max retries is exceeded"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        sftp_mock.get.side_effect = Exception("Connection error")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )
        runner.ready = True

        runner.transfer_file(max_retries=2, wait_time=1)

        # Should attempt initial + 2 retries = 3 total
        assert sftp_mock.get.call_count == 3

    @patch("file_watcher.time.sleep")
    def test_is_remote_file_ready_stable_size(self, mock_sleep, mock_ssh_client):
        """Test checking if remote file is ready with stable size"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        # Same size twice = file is ready
        sftp_mock.stat.side_effect = [
            Mock(st_size=1024),
            Mock(st_size=1024),
        ]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready(wait_time=1, max_checks=2)
        assert result is True

    @patch("file_watcher.time.sleep")
    def test_is_remote_file_ready_changing_size(self, mock_sleep, mock_ssh_client):
        """Test checking if remote file is ready with changing size"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        # Different sizes each time
        sftp_mock.stat.side_effect = [
            Mock(st_size=1024),
            Mock(st_size=2048),
            Mock(st_size=3072),
        ]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready(wait_time=1, max_checks=3)
        assert result is False

    def test_is_remote_file_ready_os_error(self, mock_ssh_client):
        """Test handling of OSError when checking remote file"""
        sftp_mock = mock_ssh_client.open_sftp.return_value
        sftp_mock.stat.side_effect = OSError("File not found")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready()
        assert result is False

    @patch("file_watcher.time.sleep")
    def test_verify_transfer_success(self, mock_sleep, mock_ssh_client):
        """Test successful transfer verification"""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        # Set up matching hashes
        filename_bytes = "test.pod5".encode("utf-8")
        expected_hash = hashlib.sha256(filename_bytes).hexdigest()
        runner.remote_hash = expected_hash

        result = runner.verify_transfer()
        assert result is True

    @patch("file_watcher.time.sleep")
    def test_verify_transfer_failure_with_retry(self, mock_sleep, mock_ssh_client):
        """Test transfer verification failure with retry"""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        # Set up non-matching hash
        runner.remote_hash = "different_hash"

        result = runner.verify_transfer(max_retries=2, wait_time=1)
        assert result is False
        assert mock_sleep.call_count == 3  # Called for each retry


class TestParseCommandLineArgs:
    """Test command line argument parsing"""

    def test_parse_args_all_provided(self):
        """Test parsing all command line arguments"""
        test_args = [
            "--watch_path",
            "/test/path",
            "--watch_pattern",
            "*.pod5",
            "--watch_duration",
            "48",
            "--host_config",
            "config.yml",
        ]

        with patch.object(sys, "argv", ["file_watcher.py"] + test_args):
            args = parse_command_line_args()

        assert args.watch_path == Path("/test/path")
        assert args.watch_pattern == "*.pod5"
        assert args.watch_duration == 48
        assert args.host_config == Path("config.yml")

    def test_parse_args_minimal(self):
        """Test parsing with only required arguments"""
        test_args = ["--host_config", "config.yml"]

        with patch.object(sys, "argv", ["file_watcher.py"] + test_args):
            args = parse_command_line_args()

        assert args.watch_path is None
        assert args.watch_pattern is None
        assert args.watch_duration == 72  # Default value
        assert args.host_config == Path("config.yml")

    def test_parse_args_short_options(self):
        """Test parsing with short option names"""
        test_args = ["-w", "/test/path", "-p", "*.bam", "-d", "24", "-c", "config.yml"]

        with patch.object(sys, "argv", ["file_watcher.py"] + test_args):
            args = parse_command_line_args()

        assert args.watch_path == Path("/test/path")
        assert args.watch_pattern == "*.bam"
        assert args.watch_duration == 24
        assert args.host_config == Path("config.yml")


class TestRuntimeConfigCheck:
    """Test runtime configuration validation"""

    def test_valid_config(self, valid_config_dict):
        """Test validation of valid configuration"""
        # Should not raise any exceptions
        runtime_config_check(valid_config_dict)

    def test_missing_watch_path(self, valid_config_dict):
        """Test validation with missing watch_path"""
        del valid_config_dict["watch_path"]
        with pytest.raises(AssertionError, match="watch_path"):
            runtime_config_check(valid_config_dict)

    def test_missing_pattern(self, valid_config_dict):
        """Test validation with missing pattern"""
        del valid_config_dict["pattern"]
        with pytest.raises(AssertionError, match="pattern"):
            runtime_config_check(valid_config_dict)

    def test_unsupported_pattern(self, valid_config_dict):
        """Test validation with unsupported file pattern"""
        valid_config_dict["pattern"] = "*.txt"
        with pytest.raises(
            AssertionError, match="file watcher currently only supports"
        ):
            runtime_config_check(valid_config_dict)

    def test_missing_host(self, valid_config_dict):
        """Test validation with missing host"""
        del valid_config_dict["host"]
        with pytest.raises(AssertionError, match="host"):
            runtime_config_check(valid_config_dict)

    def test_missing_username(self, valid_config_dict):
        """Test validation with missing username"""
        del valid_config_dict["username"]
        with pytest.raises(AssertionError, match="username"):
            runtime_config_check(valid_config_dict)

    def test_missing_password(self, valid_config_dict):
        """Test validation with missing password"""
        del valid_config_dict["password"]
        with pytest.raises(AssertionError, match="password"):
            runtime_config_check(valid_config_dict)

    def test_supported_patterns(self, valid_config_dict):
        """Test all supported file patterns"""
        for pattern in ["*.fastq.gz", "*.bam", "*.pod5"]:
            valid_config_dict["pattern"] = pattern
            runtime_config_check(valid_config_dict)  # Should not raise


class TestCredentialHelpers:
    """Test credential helper functions"""

    def test_make_credential_error(self):
        """Test error message generation"""
        error = make_credential_error("MY_VAR", "my_field", "/path/to/config.yml")
        assert "MY_VAR" in error
        assert "my_field" in error
        assert "/path/to/config.yml" in error

    @patch.dict(os.environ, {"TEST_VAR": "test_value"})
    def test_try_access_env_setting_success(self):
        """Test successful environment variable access"""
        result = try_access_env_setting("TEST_VAR", "test_field", "config.yml")
        assert result == "test_value"

    @patch.dict(os.environ, {"TEST_VAR": ""})
    def test_try_access_env_setting_empty(self):
        """Test accessing empty environment variable"""
        with pytest.raises(OSError):
            try_access_env_setting("TEST_VAR", "test_field", "config.yml")

    @patch.dict(os.environ, {}, clear=True)
    def test_try_access_env_setting_missing(self):
        """Test accessing missing environment variable"""
        with pytest.raises(OSError):
            try_access_env_setting("MISSING_VAR", "test_field", "config.yml")

    def test_try_access_config_success(self):
        """Test successful config dictionary access"""
        config_dict = {"test_field": "test_value"}
        result = try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "test_value"

    def test_try_access_config_missing(self):
        """Test accessing missing config field"""
        config_dict = {}
        with pytest.raises(OSError):
            try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)

    def test_try_access_config_empty(self):
        """Test accessing empty config field"""
        config_dict = {"test_field": ""}
        with pytest.raises(OSError):
            try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)

    @patch.dict(os.environ, {"TEST_VAR": "env_value"})
    def test_try_access_setting_env_priority(self):
        """Test that environment variable takes priority over config"""
        config_dict = {"test_field": "config_value"}
        result = try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "env_value"

    @patch.dict(os.environ, {}, clear=True)
    @patch("pathlib.Path.exists")
    def test_try_access_setting_config_fallback(self, mock_exists):
        """Test falling back to config when env var is missing"""
        mock_exists.return_value = True
        config_dict = {"test_field": "config_value"}
        result = try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "config_value"

    @patch.dict(os.environ, {}, clear=True)
    @patch("pathlib.Path.exists")
    def test_try_access_setting_config_not_exists(self, mock_exists):
        """Test error when config file doesn't exist"""
        mock_exists.return_value = False
        config_dict = {}
        with pytest.raises(
            AssertionError, match="does not point to a file that exists"
        ):
            try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)


class TestFindCredentials:
    """Test the find_credentials function"""

    @patch("builtins.open", new_callable=mock_open)
    @patch("pathlib.Path.exists")
    @patch.dict(os.environ, {}, clear=True)
    def test_find_credentials_from_config(
        self, mock_exists, mock_file, valid_config_dict
    ):
        """Test finding credentials from config file"""
        mock_exists.return_value = True
        mock_file.return_value.read.return_value = yaml.dump(valid_config_dict)

        with patch("file_watcher.yaml.safe_load", return_value=valid_config_dict):
            creds = find_credentials("config.yml")

        assert creds.watch_path == "/remote/path"
        assert creds.pattern == "*.pod5"
        assert creds.host == "192.168.1.100"
        assert creds.username == "testuser"
        assert creds.password == "testpass"
        assert creds.watch_duration == 72

    @patch.dict(
        os.environ,
        {
            "ONEROOF_WATCH_PATH": "/env/path",
            "ONEROOF_WATCH_PATTERN": "*.bam",
            "ONEROOF_WATCH_HOST": "10.0.0.1",
            "ONEROOF_WATCH_USERNAME": "envuser",
            "ONEROOF_WATCH_PASSWORD": "envpass",
            "ONEROOF_WATCH_DURATION": "48",
        },
    )
    def test_find_credentials_from_env(self):
        """Test finding credentials from environment variables"""
        creds = find_credentials("config.yml")

        assert creds.watch_path == "/env/path"
        assert creds.pattern == "*.bam"
        assert creds.host == "10.0.0.1"
        assert creds.username == "envuser"
        assert creds.password == "envpass"
        assert creds.watch_duration == 48

    @patch("builtins.open", new_callable=mock_open)
    @patch("pathlib.Path.exists")
    @patch.dict(os.environ, {"ONEROOF_WATCH_PATH": "/env/path"}, clear=True)
    def test_find_credentials_mixed(self, mock_exists, mock_file, valid_config_dict):
        """Test finding credentials from both env and config"""
        mock_exists.return_value = True
        valid_config_dict["watch_path"] = "/config/path"  # Should be overridden by env

        with patch("file_watcher.yaml.safe_load", return_value=valid_config_dict):
            creds = find_credentials("config.yml")

        assert creds.watch_path == "/env/path"  # From env
        assert creds.pattern == "*.pod5"  # From config

    @patch("pathlib.Path.exists")
    @patch.dict(os.environ, {}, clear=True)
    def test_find_credentials_missing_config(self, mock_exists):
        """Test error when config file is missing and no env vars"""
        mock_exists.return_value = False

        with pytest.raises(OSError):
            find_credentials("config.yml")


class TestMainFunction:
    """Test the main function and overall workflow"""

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    def test_main_time_limit_reached(
        self, mock_sleep, mock_time, mock_parse_args, mock_find_creds, mock_ssh_class
    ):
        """Test main function stops after time limit"""
        # Setup mocks
        mock_args = Mock(config_path="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,  # 1 hour
            watch_path="/path",
            pattern=".pod5",
        )
        mock_find_creds.return_value = mock_creds

        # Mock time to exceed duration
        mock_time.side_effect = [0, 3700]  # Start at 0, then jump past 1 hour

        mock_ssh = mock_ssh_class.return_value
        mock_ssh.exec_command.return_value = (
            None,
            Mock(read=Mock(return_value=b"")),
            None,
        )

        # Import and run main
        from file_watcher import main

        main()

        # Verify SSH connection was closed
        mock_ssh.close.assert_called_once()

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    @patch("file_watcher.TransferRunner")
    def test_main_file_transfer(
        self,
        mock_transfer_class,
        mock_sleep,
        mock_time,
        mock_parse_args,
        mock_find_creds,
        mock_ssh_class,
    ):
        """Test main function with file transfer"""
        # Setup mocks
        mock_args = Mock(config_path="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="pod5",  # Note: without the dot
        )
        mock_find_creds.return_value = mock_creds

        # Mock time to run once then exceed duration
        mock_time.side_effect = [0, 10, 3700]

        mock_ssh = mock_ssh_class.return_value
        stdout_mock = Mock()
        stdout_mock.read.return_value = b"test1.pod5\ntest2.pod5\nother.txt\n"
        mock_ssh.exec_command.return_value = (None, stdout_mock, None)

        # Mock transfer runner
        mock_runner = Mock()
        mock_runner.verify_transfer.return_value = True
        mock_transfer_class.return_value = mock_runner

        # Import and run main
        from file_watcher import main

        main()

        # Verify transfers were attempted for pod5 files
        assert mock_transfer_class.call_count == 2
        mock_runner.transfer_file.assert_called()

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.sys.exit")
    def test_main_keyboard_interrupt(
        self, mock_exit, mock_parse_args, mock_find_creds, mock_ssh_class
    ):
        """Test handling of keyboard interrupt"""
        # Setup mocks
        mock_args = Mock(config_path="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="pod5",
        )
        mock_find_creds.return_value = mock_creds

        mock_ssh = mock_ssh_class.return_value
        mock_ssh.exec_command.side_effect = KeyboardInterrupt()

        # Import and run main
        from file_watcher import main

        main()

        # Verify clean exit
        mock_ssh.close.assert_called_once()
        mock_exit.assert_called_once_with(1)

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    @patch("file_watcher.TransferRunner")
    def test_main_corrupted_transfer(
        self,
        mock_transfer_class,
        mock_sleep,
        mock_time,
        mock_parse_args,
        mock_find_creds,
        mock_ssh_class,
    ):
        """Test handling of corrupted file transfer"""
        # Setup mocks
        mock_args = Mock(config_path="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="pod5",
        )
        mock_find_creds.return_value = mock_creds

        # Mock time to run once then exceed duration
        mock_time.side_effect = [0, 10, 3700]

        mock_ssh = mock_ssh_class.return_value
        stdout_mock = Mock()
        stdout_mock.read.return_value = b"test.pod5\n"
        mock_ssh.exec_command.return_value = (None, stdout_mock, None)

        # Mock transfer runner with failed verification
        mock_runner = Mock()
        mock_runner.filename = "test.pod5"
        mock_runner.verify_transfer.return_value = False
        mock_transfer_class.return_value = mock_runner

        # Import and run main with captured logs
        from file_watcher import main

        with patch("file_watcher.logger") as mock_logger:
            main()

            # Verify error was logged for corrupted transfer
            mock_logger.error.assert_called_with(
                "The file, test.pod5, was corrupted during the transfer process."
            )


# Additional edge case tests
class TestEdgeCases:
    """Test edge cases and error conditions"""

    @patch("file_watcher.time.sleep")
    def test_transfer_runner_str_representation(self, mock_sleep, mock_ssh_client):
        """Test string representation of TransferRunner"""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        # The string representation is used for hashing
        str_repr = str(runner)
        assert isinstance(str_repr, str)

    def test_empty_file_list(self, mock_ssh_client):
        """Test handling of empty file list from remote"""
        stdout_mock = Mock()
        stdout_mock.read.return_value = b""
        mock_ssh_client.exec_command.return_value = (None, stdout_mock, None)

        # This should not raise any exceptions
        all_files = stdout_mock.read().decode("utf8").splitlines()
        assert all_files == []

    @patch("builtins.open", side_effect=IOError("Permission denied"))
    @patch("pathlib.Path.exists", return_value=True)
    def test_config_file_permission_error(self, mock_exists, mock_open):
        """Test handling of permission error when reading config"""
        with pytest.raises(IOError, match="Permission denied"):
            find_credentials("config.yml")

    @patch("file_watcher.yaml.safe_load", side_effect=yaml.YAMLError("Invalid YAML"))
    @patch("builtins.open", new_callable=mock_open)
    @patch("pathlib.Path.exists", return_value=True)
    def test_invalid_yaml_config(self, mock_exists, mock_file, mock_yaml):
        """Test handling of invalid YAML in config file"""
        with pytest.raises(yaml.YAMLError, match="Invalid YAML"):
            find_credentials("config.yml")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
