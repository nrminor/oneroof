#!/usr/bin/env python3
"""
Tests for file_watcher.py

Tests the SFTP file transfer client including credentials handling,
transfer operations, and configuration parsing.
"""

from __future__ import annotations

import os
import sys
import tempfile
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import Mock, mock_open, patch

import pytest
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
from paramiko.client import SSHClient

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def mock_ssh_client():
    """Mock SSH client with SFTP support."""
    client = Mock(spec=SSHClient)
    sftp = Mock()
    client.open_sftp.return_value = sftp
    return client


@pytest.fixture
def valid_config_dict():
    """Valid configuration dictionary with all required fields."""
    return {
        "watch_path": "/remote/data",
        "pattern": "*.pod5",
        "host": "192.168.1.100",
        "username": "testuser",
        "password": "testpass",
        "watch_duration": 72,
    }


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


# =============================================================================
# Test Credentials Dataclass
# =============================================================================


class TestCredentials:
    """Tests for the Credentials dataclass."""

    def test_creation(self):
        """Credentials can be created with all required fields."""
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

    def test_frozen(self):
        """Credentials is immutable (frozen)."""
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


# =============================================================================
# Test TransferRunner
# =============================================================================


class TestTransferRunnerInit:
    """Tests for TransferRunner initialization."""

    def test_init_sets_attributes(self, mock_ssh_client):
        """__init__ sets all attributes correctly."""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path="/local/path",
        )

        assert runner.client is mock_ssh_client
        assert runner.filename == "test.pod5"
        assert runner.remote_path == "/remote/path"
        assert runner.address == "192.168.1.100"
        assert runner.username == "testuser"
        assert runner.local_path == "/local/path"

    def test_init_defaults(self, mock_ssh_client):
        """__init__ sets correct default values."""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        assert runner.ready is False
        assert runner.remote_hash is None
        assert runner._temp_path is None
        assert runner.local_path == Path.cwd()


class TestTransferRunnerPrepare:
    """Tests for TransferRunner.prepare()."""

    @patch("file_watcher.time.sleep")
    def test_prepare_returns_false_if_file_exists_locally(
        self, mock_sleep, mock_ssh_client, temp_dir
    ):
        """prepare() returns False if file already exists locally."""
        # Create local file
        local_file = temp_dir / "test.pod5"
        local_file.touch()

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        result = runner.prepare()

        assert result is False
        assert runner.ready is False

    @patch("file_watcher.time.sleep")
    def test_prepare_returns_false_if_remote_not_ready(
        self, mock_sleep, mock_ssh_client, temp_dir
    ):
        """prepare() returns False if remote file size is still changing."""
        sftp = mock_ssh_client.open_sftp.return_value
        # Size keeps changing - file not ready (need enough values for max_checks)
        sftp.stat.side_effect = [Mock(st_size=i * 100) for i in range(15)]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        result = runner.prepare()

        assert result is False
        assert runner.ready is False

    @patch("file_watcher.time.sleep")
    def test_prepare_returns_true_when_ready(
        self, mock_sleep, mock_ssh_client, temp_dir
    ):
        """prepare() returns True when file is ready and not local."""
        sftp = mock_ssh_client.open_sftp.return_value
        # Size stable - file ready
        sftp.stat.side_effect = [Mock(st_size=1024), Mock(st_size=1024)]
        # Mock remote file for hash computation
        mock_file = Mock()
        mock_file.read.side_effect = [b"test content", b""]
        sftp.open.return_value.__enter__ = Mock(return_value=mock_file)
        sftp.open.return_value.__exit__ = Mock(return_value=False)

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        result = runner.prepare()

        assert result is True
        assert runner.ready is True
        assert runner.remote_hash is not None


class TestTransferRunnerIsRemoteFileReady:
    """Tests for TransferRunner.is_remote_file_ready()."""

    @patch("file_watcher.time.sleep")
    def test_returns_true_when_size_stable(self, mock_sleep, mock_ssh_client):
        """Returns True when file size remains unchanged."""
        sftp = mock_ssh_client.open_sftp.return_value
        sftp.stat.side_effect = [Mock(st_size=1024), Mock(st_size=1024)]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready(wait_time=0, max_checks=2)

        assert result is True
        sftp.close.assert_called_once()

    @patch("file_watcher.time.sleep")
    def test_returns_false_when_size_changing(self, mock_sleep, mock_ssh_client):
        """Returns False when file size keeps changing."""
        sftp = mock_ssh_client.open_sftp.return_value
        sftp.stat.side_effect = [
            Mock(st_size=100),
            Mock(st_size=200),
            Mock(st_size=300),
        ]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready(wait_time=0, max_checks=3)

        assert result is False
        sftp.close.assert_called_once()

    def test_returns_false_on_os_error(self, mock_ssh_client):
        """Returns False when OSError occurs."""
        sftp = mock_ssh_client.open_sftp.return_value
        sftp.stat.side_effect = OSError("File not found")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner.is_remote_file_ready()

        assert result is False
        sftp.close.assert_called_once()


class TestTransferRunnerTransferFile:
    """Tests for TransferRunner.transfer_file()."""

    @patch("file_watcher.time.sleep")
    def test_downloads_to_temp_file(self, mock_sleep, mock_ssh_client, temp_dir):
        """transfer_file() downloads to a .tmp file."""
        sftp = mock_ssh_client.open_sftp.return_value

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        runner.transfer_file()

        # Check sftp.get was called with temp path
        call_args = sftp.get.call_args
        assert call_args is not None
        remote_path, local_path = call_args[0]
        assert remote_path == "/remote/path/test.pod5"
        assert str(local_path).endswith(".tmp")
        assert runner._temp_path == local_path
        sftp.close.assert_called()

    @patch("file_watcher.time.sleep")
    def test_retries_on_failure(self, mock_sleep, mock_ssh_client, temp_dir):
        """transfer_file() retries on failure."""
        sftp = mock_ssh_client.open_sftp.return_value
        sftp.get.side_effect = [Exception("Connection error"), None]

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        runner.transfer_file(max_retries=2, wait_time=0)

        assert sftp.get.call_count == 2

    @patch("file_watcher.time.sleep")
    def test_stops_after_max_retries(self, mock_sleep, mock_ssh_client, temp_dir):
        """transfer_file() stops retrying after max_retries."""
        sftp = mock_ssh_client.open_sftp.return_value
        sftp.get.side_effect = Exception("Connection error")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )

        runner.transfer_file(max_retries=2, wait_time=0)

        # Initial attempt + 2 retries = 3 total
        assert sftp.get.call_count == 3


class TestTransferRunnerVerifyTransfer:
    """Tests for TransferRunner.verify_transfer()."""

    def test_returns_false_if_no_temp_file(self, mock_ssh_client, temp_dir):
        """verify_transfer() returns False if temp file doesn't exist."""
        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )
        runner._temp_path = None

        result = runner.verify_transfer()

        assert result is False

    def test_returns_true_and_renames_on_match(self, mock_ssh_client, temp_dir):
        """verify_transfer() renames temp file when hashes match."""
        # Create temp file with known content
        temp_file = temp_dir / "test.pod5.tmp"
        temp_file.write_bytes(b"test content")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )
        runner._temp_path = temp_file
        # Set remote_hash to match the temp file content
        runner.remote_hash = runner._compute_file_hash(temp_file)

        result = runner.verify_transfer()

        assert result is True
        assert not temp_file.exists()
        assert (temp_dir / "test.pod5").exists()

    def test_returns_false_and_deletes_on_mismatch(self, mock_ssh_client, temp_dir):
        """verify_transfer() deletes temp file when hashes don't match."""
        # Create temp file
        temp_file = temp_dir / "test.pod5.tmp"
        temp_file.write_bytes(b"test content")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
            local_path=temp_dir,
        )
        runner._temp_path = temp_file
        runner.remote_hash = "different_hash"

        result = runner.verify_transfer()

        assert result is False
        assert not temp_file.exists()
        assert not (temp_dir / "test.pod5").exists()


class TestTransferRunnerComputeFileHash:
    """Tests for TransferRunner._compute_file_hash()."""

    def test_computes_correct_hash(self, mock_ssh_client, temp_dir):
        """_compute_file_hash() computes correct SHA256 hash."""
        test_file = temp_dir / "test.txt"
        test_file.write_bytes(b"hello world")

        runner = TransferRunner(
            client=mock_ssh_client,
            filename="test.pod5",
            remote_path="/remote/path",
            address="192.168.1.100",
            username="testuser",
        )

        result = runner._compute_file_hash(test_file)

        # Known SHA256 of "hello world"
        expected = "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        assert result == expected


# =============================================================================
# Test parse_command_line_args
# =============================================================================


class TestParseCommandLineArgs:
    """Tests for parse_command_line_args()."""

    def test_parses_host_config(self):
        """Parses --host_config argument."""
        with patch.object(
            sys, "argv", ["file_watcher.py", "--host_config", "config.yml"]
        ):
            args = parse_command_line_args()

        assert args.host_config == Path("config.yml")

    def test_parses_short_option(self):
        """Parses -c short option."""
        with patch.object(sys, "argv", ["file_watcher.py", "-c", "config.yml"]):
            args = parse_command_line_args()

        assert args.host_config == Path("config.yml")

    def test_requires_host_config(self):
        """Exits if --host_config is not provided."""
        with patch.object(sys, "argv", ["file_watcher.py"]):
            with pytest.raises(SystemExit):
                parse_command_line_args()


# =============================================================================
# Test runtime_config_check
# =============================================================================


class TestRuntimeConfigCheck:
    """Tests for runtime_config_check()."""

    def test_valid_config_passes(self, valid_config_dict):
        """Valid config passes without error."""
        runtime_config_check(valid_config_dict)  # Should not raise

    def test_missing_watch_path(self, valid_config_dict):
        """Raises on missing watch_path."""
        del valid_config_dict["watch_path"]
        with pytest.raises(AssertionError, match="watch_path"):
            runtime_config_check(valid_config_dict)

    def test_missing_pattern(self, valid_config_dict):
        """Raises on missing pattern."""
        del valid_config_dict["pattern"]
        with pytest.raises(AssertionError, match="pattern"):
            runtime_config_check(valid_config_dict)

    def test_unsupported_pattern(self, valid_config_dict):
        """Raises on unsupported pattern."""
        valid_config_dict["pattern"] = "*.txt"
        with pytest.raises(AssertionError, match="currently only supports"):
            runtime_config_check(valid_config_dict)

    def test_missing_host(self, valid_config_dict):
        """Raises on missing host."""
        del valid_config_dict["host"]
        with pytest.raises(AssertionError, match="host"):
            runtime_config_check(valid_config_dict)

    def test_missing_username(self, valid_config_dict):
        """Raises on missing username."""
        del valid_config_dict["username"]
        with pytest.raises(AssertionError, match="username"):
            runtime_config_check(valid_config_dict)

    def test_missing_password(self, valid_config_dict):
        """Raises on missing password."""
        del valid_config_dict["password"]
        with pytest.raises(AssertionError, match="password"):
            runtime_config_check(valid_config_dict)

    def test_missing_watch_duration(self, valid_config_dict):
        """Raises on missing watch_duration."""
        del valid_config_dict["watch_duration"]
        with pytest.raises(AssertionError, match="watch_duration"):
            runtime_config_check(valid_config_dict)

    def test_invalid_watch_duration_type(self, valid_config_dict):
        """Raises on non-integer watch_duration."""
        valid_config_dict["watch_duration"] = "72"  # String instead of int
        with pytest.raises(AssertionError, match="positive integer"):
            runtime_config_check(valid_config_dict)

    def test_invalid_watch_duration_value(self, valid_config_dict):
        """Raises on non-positive watch_duration."""
        valid_config_dict["watch_duration"] = 0
        with pytest.raises(AssertionError, match="positive integer"):
            runtime_config_check(valid_config_dict)

    def test_all_supported_patterns(self, valid_config_dict):
        """All supported patterns pass validation."""
        for pattern in ["*.fastq.gz", "*.bam", "*.pod5"]:
            valid_config_dict["pattern"] = pattern
            runtime_config_check(valid_config_dict)  # Should not raise


# =============================================================================
# Test credential helper functions
# =============================================================================


class TestMakeCredentialError:
    """Tests for make_credential_error()."""

    def test_generates_error_message(self):
        """Generates error message with all components."""
        error = make_credential_error("MY_VAR", "my_field", "/path/to/config.yml")

        assert "MY_VAR" in error
        assert "my_field" in error
        assert "/path/to/config.yml" in error


class TestTryAccessEnvSetting:
    """Tests for try_access_env_setting()."""

    @patch.dict(os.environ, {"TEST_VAR": "test_value"})
    def test_returns_env_value(self):
        """Returns environment variable value."""
        result = try_access_env_setting("TEST_VAR", "test_field", "config.yml")
        assert result == "test_value"

    @patch.dict(os.environ, {"TEST_VAR": ""})
    def test_raises_on_empty(self):
        """Raises OSError on empty environment variable."""
        with pytest.raises(OSError):
            try_access_env_setting("TEST_VAR", "test_field", "config.yml")

    @patch.dict(os.environ, {}, clear=True)
    def test_raises_on_missing(self):
        """Raises OSError on missing environment variable."""
        with pytest.raises(OSError):
            try_access_env_setting("MISSING_VAR", "test_field", "config.yml")


class TestTryAccessConfig:
    """Tests for try_access_config()."""

    def test_returns_config_value(self):
        """Returns config dictionary value."""
        config_dict = {"test_field": "test_value"}
        result = try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "test_value"

    def test_raises_on_missing(self):
        """Raises OSError on missing config field."""
        config_dict = {}
        with pytest.raises(OSError):
            try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)

    def test_raises_on_empty(self):
        """Raises OSError on empty config field."""
        config_dict = {"test_field": ""}
        with pytest.raises(OSError):
            try_access_config("TEST_VAR", "test_field", "config.yml", config_dict)


class TestTryAccessSetting:
    """Tests for try_access_setting()."""

    @patch.dict(os.environ, {"TEST_VAR": "env_value"})
    def test_env_takes_priority(self):
        """Environment variable takes priority over config."""
        config_dict = {"test_field": "config_value"}
        result = try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "env_value"

    @patch.dict(os.environ, {}, clear=True)
    @patch("pathlib.Path.exists")
    def test_falls_back_to_config(self, mock_exists):
        """Falls back to config when env var missing."""
        mock_exists.return_value = True
        config_dict = {"test_field": "config_value"}
        result = try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)
        assert result == "config_value"

    @patch.dict(os.environ, {}, clear=True)
    @patch("pathlib.Path.exists")
    def test_raises_if_config_not_exists(self, mock_exists):
        """Raises AssertionError if config file doesn't exist."""
        mock_exists.return_value = False
        config_dict = {}
        with pytest.raises(AssertionError, match="does not point to a file"):
            try_access_setting("TEST_VAR", "test_field", "config.yml", config_dict)


# =============================================================================
# Test find_credentials
# =============================================================================


class TestFindCredentials:
    """Tests for find_credentials()."""

    @patch("pathlib.Path.open", new_callable=mock_open)
    @patch("pathlib.Path.exists")
    @patch.dict(os.environ, {}, clear=True)
    def test_loads_from_config_file(
        self, mock_exists, mock_path_open, valid_config_dict
    ):
        """Loads credentials from config file."""
        mock_exists.return_value = True

        with patch("file_watcher.yaml.safe_load", return_value=valid_config_dict):
            creds = find_credentials("config.yml")

        assert creds.watch_path == "/remote/data"
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
    def test_loads_from_env_vars(self):
        """Loads credentials from environment variables."""
        creds = find_credentials("config.yml")

        assert creds.watch_path == "/env/path"
        assert creds.pattern == "*.bam"
        assert creds.host == "10.0.0.1"
        assert creds.username == "envuser"
        assert creds.password == "envpass"
        assert creds.watch_duration == 48


# =============================================================================
# Test main function
# =============================================================================


class TestMain:
    """Tests for main()."""

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    def test_stops_after_duration(
        self, mock_sleep, mock_time, mock_parse_args, mock_find_creds, mock_ssh_class
    ):
        """main() stops after watch duration expires."""
        mock_args = Mock(host_config="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,  # 1 hour
            watch_path="/path",
            pattern="*.pod5",
        )
        mock_find_creds.return_value = mock_creds

        # Time jumps past duration on second call
        mock_time.side_effect = [0, 3700]

        mock_ssh = mock_ssh_class.return_value
        mock_ssh.exec_command.return_value = (
            None,
            Mock(read=Mock(return_value=b"")),
            None,
        )

        from file_watcher import main

        main()

        mock_ssh.close.assert_called_once()

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.sys.exit")
    def test_handles_keyboard_interrupt(
        self, mock_exit, mock_parse_args, mock_find_creds, mock_ssh_class
    ):
        """main() handles KeyboardInterrupt gracefully."""
        mock_args = Mock(host_config="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="*.pod5",
        )
        mock_find_creds.return_value = mock_creds

        mock_ssh = mock_ssh_class.return_value
        mock_ssh.exec_command.side_effect = KeyboardInterrupt()

        from file_watcher import main

        main()

        mock_ssh.close.assert_called_once()
        mock_exit.assert_called_once_with(1)

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    @patch("file_watcher.TransferRunner")
    def test_calls_prepare_before_transfer(
        self,
        mock_runner_class,
        mock_sleep,
        mock_time,
        mock_parse_args,
        mock_find_creds,
        mock_ssh_class,
    ):
        """main() calls prepare() before transfer_file()."""
        mock_args = Mock(host_config="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="*.pod5",
        )
        mock_find_creds.return_value = mock_creds

        mock_time.side_effect = [0, 10, 3700]

        mock_ssh = mock_ssh_class.return_value
        stdout_mock = Mock()
        stdout_mock.read.return_value = b"test.pod5\n"
        mock_ssh.exec_command.return_value = (None, stdout_mock, None)

        mock_runner = Mock()
        mock_runner.prepare.return_value = True
        mock_runner.verify_transfer.return_value = True
        mock_runner_class.return_value = mock_runner

        from file_watcher import main

        main()

        mock_runner.prepare.assert_called()
        mock_runner.transfer_file.assert_called()

    @patch("file_watcher.SSHClient")
    @patch("file_watcher.find_credentials")
    @patch("file_watcher.parse_command_line_args")
    @patch("file_watcher.time.time")
    @patch("file_watcher.time.sleep")
    @patch("file_watcher.TransferRunner")
    def test_skips_transfer_if_prepare_fails(
        self,
        mock_runner_class,
        mock_sleep,
        mock_time,
        mock_parse_args,
        mock_find_creds,
        mock_ssh_class,
    ):
        """main() skips transfer if prepare() returns False."""
        mock_args = Mock(host_config="config.yml")
        mock_parse_args.return_value = mock_args

        mock_creds = Mock(
            host="localhost",
            username="user",
            password="pass",
            watch_duration=1,
            watch_path="/path",
            pattern="*.pod5",
        )
        mock_find_creds.return_value = mock_creds

        mock_time.side_effect = [0, 10, 3700]

        mock_ssh = mock_ssh_class.return_value
        stdout_mock = Mock()
        stdout_mock.read.return_value = b"test.pod5\n"
        mock_ssh.exec_command.return_value = (None, stdout_mock, None)

        mock_runner = Mock()
        mock_runner.prepare.return_value = False  # Not ready
        mock_runner_class.return_value = mock_runner

        from file_watcher import main

        main()

        mock_runner.prepare.assert_called()
        mock_runner.transfer_file.assert_not_called()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
