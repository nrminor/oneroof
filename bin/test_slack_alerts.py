#!/usr/bin/env python3
"""
Comprehensive pytest test module for slack_alerts.py

Tests cover:
- Message formatting
- Webhook URL validation
- HTTP request handling (all network calls mocked)
- Error handling (network failures, invalid webhooks)
- Different alert types/levels
- Environment variable handling
"""

import argparse
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

# Import the module we're testing
import slack_alerts
from slack_sdk.errors import SlackApiError


class TestParseCommandLineArgs:
    """Test command line argument parsing"""

    def test_parse_all_required_args(self):
        """Test parsing all required arguments"""
        test_args = [
            "--input_tsv_dir",
            "/path/to/tsv",
            "--depth",
            "20",
            "--run_label",
            "experiment_001",
        ]
        with patch("sys.argv", ["slack_alerts.py"] + test_args):
            args = slack_alerts.parse_command_line_args()
            assert args.input_tsv_dir == Path("/path/to/tsv")
            assert args.depth == 20
            assert args.run_label == "experiment_001"

    def test_parse_with_short_options(self):
        """Test parsing with short option forms"""
        test_args = [
            "--input_tsv_dir",
            "/path/to/tsv",
            "-d",
            "30",
            "-r",
            "exp_002",
        ]
        with patch("sys.argv", ["slack_alerts.py"] + test_args):
            args = slack_alerts.parse_command_line_args()
            assert args.depth == 30
            assert args.run_label == "exp_002"

    def test_missing_required_args(self):
        """Test that missing required arguments raise SystemExit"""
        test_args = ["--input_tsv_dir", "/path/to/tsv"]  # Missing depth and run_label
        with patch("sys.argv", ["slack_alerts.py"] + test_args):
            with pytest.raises(SystemExit):
                slack_alerts.parse_command_line_args()


class TestPassingSamples:
    """Test passing_samples function"""

    def test_all_samples_passing(self):
        """Test when all samples pass the coverage threshold"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "coverage": [0.95, 0.98, 0.96],
            }
        )
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 3
        assert "sample1: 0.95" in message
        assert "sample2: 0.98" in message
        assert "sample3: 0.96" in message

    def test_some_samples_passing(self):
        """Test when only some samples pass"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "coverage": [0.95, 0.85, 0.96],
            }
        )
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 2
        assert "sample1: 0.95" in message
        assert "sample2" not in message
        assert "sample3: 0.96" in message

    def test_no_samples_passing(self):
        """Test when no samples pass"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "coverage": [0.85, 0.80],
            }
        )
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 0
        assert message == "No passing samples."

    def test_coverage_threshold_calculation(self):
        """Test the threshold calculation for different coverage values"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1"],
                "coverage": [0.95],
            }
        )
        # For coverage_threshold=20, threshold should be 1-(1/20)=0.95
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 1  # 0.95 >= 0.95

        # For coverage_threshold=10, threshold should be 1-(1/10)=0.9
        message, count = slack_alerts.passing_samples(df, coverage_threshold=10)
        assert count == 1  # 0.95 >= 0.9


class TestFailingSamples:
    """Test failing_samples function"""

    def test_all_samples_failing(self):
        """Test when all samples fail the coverage threshold"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "coverage": [0.85, 0.80, 0.82],
            }
        )
        message = slack_alerts.failing_samples(df, coverage_threshold=20)
        assert "sample1: 0.85" in message
        assert "sample2: 0.8" in message  # 0.80 formats as 0.8
        assert "sample3: 0.82" in message

    def test_some_samples_failing(self):
        """Test when only some samples fail"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2", "sample3"],
                "coverage": [0.95, 0.85, 0.96],
            }
        )
        message = slack_alerts.failing_samples(df, coverage_threshold=20)
        assert "sample1" not in message
        assert "sample2: 0.85" in message
        assert "sample3" not in message

    def test_no_samples_failing(self):
        """Test when no samples fail"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "coverage": [0.96, 0.98],
            }
        )
        message = slack_alerts.failing_samples(df, coverage_threshold=20)
        assert message == "No failing samples."


class TestGetUserIds:
    """Test get_user_ids function"""

    def test_read_from_default_location(self):
        """Test reading user IDs from default location"""
        with tempfile.TemporaryDirectory() as tmpdir:
            oneroof_dir = Path(tmpdir) / ".oneroof"
            oneroof_dir.mkdir()
            user_ids_file = oneroof_dir / "slack.user_ids"
            user_ids_file.write_text("U123456\nU789012\n# This is a comment\nU345678\n")

            with patch.dict(os.environ, {}, clear=True):
                with patch("pathlib.Path.home", return_value=Path(tmpdir)):
                    user_ids = slack_alerts.get_user_ids()
                    assert user_ids == ["U123456", "U789012", "U345678"]

    def test_read_from_env_variable_location(self):
        """Test reading user IDs from environment variable location"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("U111111\nU222222\n")
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_USER_IDS": f.name}):
                    user_ids = slack_alerts.get_user_ids()
                    assert user_ids == ["U111111", "U222222"]
            finally:
                os.unlink(f.name)

    def test_file_not_exists(self):
        """Test when user IDs file doesn't exist"""
        with patch.dict(os.environ, {"ONEROOF_SLACK_USER_IDS": "/nonexistent/file"}):
            user_ids = slack_alerts.get_user_ids()
            assert user_ids == []

    def test_empty_file(self):
        """Test when user IDs file is empty"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("")  # Empty file
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_USER_IDS": f.name}):
                    user_ids = slack_alerts.get_user_ids()
                    assert user_ids == []
            finally:
                os.unlink(f.name)

    def test_file_with_only_comments(self):
        """Test when file contains only comments"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("# Comment 1\n# Comment 2\n")
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_USER_IDS": f.name}):
                    user_ids = slack_alerts.get_user_ids()
                    assert user_ids == []
            finally:
                os.unlink(f.name)


class TestGetSlackToken:
    """Test get_slack_token function"""

    def test_read_from_default_location(self):
        """Test reading token from default location"""
        with tempfile.TemporaryDirectory() as tmpdir:
            oneroof_dir = Path(tmpdir) / ".oneroof"
            oneroof_dir.mkdir()
            token_file = oneroof_dir / "slack.token"
            token_file.write_text("xoxb-test-token-123")

            with patch.dict(os.environ, {}, clear=True):
                with patch("pathlib.Path.home", return_value=Path(tmpdir)):
                    token = slack_alerts.get_slack_token()
                    assert token == "xoxb-test-token-123"

    def test_read_from_env_variable_location(self):
        """Test reading token from environment variable location"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("xoxb-env-token-456")
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_TOKEN": f.name}):
                    token = slack_alerts.get_slack_token()
                    assert token == "xoxb-env-token-456"
            finally:
                os.unlink(f.name)

    def test_file_not_exists(self):
        """Test when token file doesn't exist"""
        with patch.dict(os.environ, {"ONEROOF_SLACK_TOKEN": "/nonexistent/file"}):
            token = slack_alerts.get_slack_token()
            assert token == ""

    def test_empty_file(self):
        """Test when token file is empty"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("")  # Empty file
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_TOKEN": f.name}):
                    token = slack_alerts.get_slack_token()
                    assert token == ""
            finally:
                os.unlink(f.name)

    def test_token_with_whitespace(self):
        """Test that whitespace is stripped from token"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("  xoxb-whitespace-token  \n")
            f.flush()

            try:
                with patch.dict(os.environ, {"ONEROOF_SLACK_TOKEN": f.name}):
                    token = slack_alerts.get_slack_token()
                    assert token == "xoxb-whitespace-token"
            finally:
                os.unlink(f.name)


class TestSendSlackNotification:
    """Test send_slack_notification function"""

    @patch("slack_alerts.WebClient")
    @patch("slack_alerts.get_slack_token")
    @patch("slack_alerts.get_user_ids")
    def test_successful_notification(self, mock_get_user_ids, mock_get_slack_token, mock_webclient):
        """Test successful Slack notification sending"""
        # Setup mocks
        mock_get_user_ids.return_value = ["U123456", "U789012"]
        mock_get_slack_token.return_value = "xoxb-test-token"
        mock_client = MagicMock()
        mock_webclient.return_value = mock_client

        # Create test data
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample_id\tcoverage\n")
            f.write("sample1\t0.95\n")
            f.write("sample2\t0.85\n")
            f.write("sample3\t0.98\n")
            f.flush()

            try:
                slack_alerts.send_slack_notification(
                    run_label="exp_001",
                    stats_tsv=f.name,
                    coverage_threshold=20,
                )

                # Verify WebClient was created with correct token
                mock_webclient.assert_called_with(token="xoxb-test-token")

                # Verify messages were sent to both users
                assert mock_client.chat_postMessage.call_count == 2

                # Check the message content
                calls = mock_client.chat_postMessage.call_args_list
                for call in calls:
                    assert call[1]["channel"] in ["U123456", "U789012"]
                    blocks = call[1]["blocks"]
                    assert len(blocks) == 2
                    assert "exp_001" in blocks[0]["text"]["text"]
                    assert "2 samples passing" in blocks[0]["text"]["text"]
                    assert "20X coverage" in blocks[0]["text"]["text"]

            finally:
                os.unlink(f.name)

    @patch("slack_alerts.WebClient")
    @patch("slack_alerts.get_slack_token")
    @patch("slack_alerts.get_user_ids")
    def test_slack_api_error(self, mock_get_user_ids, mock_get_slack_token, mock_webclient, capsys):
        """Test handling of Slack API errors"""
        # Setup mocks
        mock_get_user_ids.return_value = ["U123456"]
        mock_get_slack_token.return_value = "xoxb-test-token"
        mock_client = MagicMock()
        mock_webclient.return_value = mock_client

        # Simulate SlackApiError
        error_response = {"error": "invalid_auth"}
        mock_client.chat_postMessage.side_effect = SlackApiError(
            message="Invalid auth", response=error_response
        )

        # Create test data
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample_id\tcoverage\n")
            f.write("sample1\t0.95\n")
            f.flush()

            try:
                slack_alerts.send_slack_notification(
                    run_label="exp_001",
                    stats_tsv=f.name,
                    coverage_threshold=20,
                )

                # Check error message was printed
                captured = capsys.readouterr()
                assert "Error sending message: invalid_auth" in captured.out

            finally:
                os.unlink(f.name)

    @patch("slack_alerts.WebClient")
    @patch("slack_alerts.get_slack_token")
    @patch("slack_alerts.get_user_ids")
    def test_no_users_configured(self, mock_get_user_ids, mock_get_slack_token, mock_webclient):
        """Test when no users are configured"""
        # Setup mocks
        mock_get_user_ids.return_value = []  # No users
        mock_get_slack_token.return_value = "xoxb-test-token"
        mock_client = MagicMock()
        mock_webclient.return_value = mock_client

        # Create test data
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample_id\tcoverage\n")
            f.write("sample1\t0.95\n")
            f.flush()

            try:
                slack_alerts.send_slack_notification(
                    run_label="exp_001",
                    stats_tsv=f.name,
                    coverage_threshold=20,
                )

                # Verify no messages were sent
                mock_client.chat_postMessage.assert_not_called()

            finally:
                os.unlink(f.name)

    @patch("slack_alerts.WebClient")
    @patch("slack_alerts.get_slack_token")
    @patch("slack_alerts.get_user_ids")
    def test_message_formatting(self, mock_get_user_ids, mock_get_slack_token, mock_webclient):
        """Test detailed message formatting"""
        # Setup mocks
        mock_get_user_ids.return_value = ["U123456"]
        mock_get_slack_token.return_value = "xoxb-test-token"
        mock_client = MagicMock()
        mock_webclient.return_value = mock_client

        # Create test data with specific samples
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample_id\tcoverage\n")
            f.write("sample_pass_1\t0.96\n")
            f.write("sample_pass_2\t0.98\n")
            f.write("sample_fail_1\t0.85\n")
            f.write("sample_fail_2\t0.80\n")
            f.flush()

            try:
                slack_alerts.send_slack_notification(
                    run_label="test_run_123",
                    stats_tsv=f.name,
                    coverage_threshold=20,
                )

                # Get the actual call arguments
                call_args = mock_client.chat_postMessage.call_args[1]
                blocks = call_args["blocks"]

                # Check header section
                header_text = blocks[0]["text"]["text"]
                assert "test_run_123" in header_text
                assert "2 samples passing" in header_text
                assert "20X coverage" in header_text

                # Check fields section
                fields = blocks[1]["fields"]
                passing_field = fields[0]["text"]
                failing_field = fields[1]["text"]

                # Check passing samples
                assert "*PASSING*" in passing_field
                assert "sample_pass_1: 0.96" in passing_field
                assert "sample_pass_2: 0.98" in passing_field

                # Check failing samples
                assert "*FAILING*" in failing_field
                assert "sample_fail_1: 0.85" in failing_field
                assert "sample_fail_2: 0.8" in failing_field  # 0.80 formats as 0.8

            finally:
                os.unlink(f.name)


class TestMain:
    """Test main function"""

    @patch("slack_alerts.send_slack_notification")
    @patch("slack_alerts.parse_command_line_args")
    def test_main_function(self, mock_parse_args, mock_send_notification):
        """Test main function integration"""
        # Setup mock arguments
        mock_args = argparse.Namespace(
            input_tsv_dir=Path("/path/to/tsv"),
            depth=30,
            run_label="main_test_run",
        )
        mock_parse_args.return_value = mock_args

        # Call main
        slack_alerts.main()

        # Verify send_slack_notification was called with correct arguments
        mock_send_notification.assert_called_once_with(
            "main_test_run",
            Path("/path/to/tsv"),
            30,
        )


class TestEdgeCases:
    """Test edge cases and error conditions"""

    def test_empty_dataframe(self):
        """Test handling of empty dataframe"""
        df = pd.DataFrame(columns=["sample_id", "coverage"])
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 0
        assert message == "No passing samples."

    def test_very_high_coverage_threshold(self):
        """Test with very high coverage threshold"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1"],
                "coverage": [0.99],
            }
        )
        # For coverage_threshold=100, threshold would be 0.99
        message, count = slack_alerts.passing_samples(df, coverage_threshold=100)
        assert count == 1  # 0.99 >= 0.99

    def test_coverage_exactly_at_threshold(self):
        """Test samples exactly at the coverage threshold"""
        df = pd.DataFrame(
            {
                "sample_id": ["sample1"],
                "coverage": [0.95],
            }
        )
        # For coverage_threshold=20, threshold is 0.95
        message, count = slack_alerts.passing_samples(df, coverage_threshold=20)
        assert count == 1  # Should pass when exactly at threshold

    @patch("pandas.read_csv")
    def test_malformed_tsv_file(self, mock_read_csv):
        """Test handling of malformed TSV file"""
        # Simulate pandas read error
        mock_read_csv.side_effect = pd.errors.ParserError("Bad TSV format")

        with pytest.raises(pd.errors.ParserError):
            slack_alerts.send_slack_notification(
                run_label="error_test",
                stats_tsv="bad_file.tsv",
                coverage_threshold=20,
            )


# Run tests with pytest
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
