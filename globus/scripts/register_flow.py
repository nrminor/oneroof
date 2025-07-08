#!/usr/bin/env python3
"""
Register OneRoof Flow with Globus Flows service.

This script registers the flow definition and sets up necessary permissions.
"""

import json
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
GLOBUS_DIR = SCRIPT_DIR.parent
FLOW_FILE = GLOBUS_DIR / "flows" / "oneroof_flow.json"
CONFIG_FILE = GLOBUS_DIR / "config" / ".env"


def load_env_config() -> dict[str, str]:
    """Load configuration from .env file."""
    config = {}

    if not CONFIG_FILE.exists():
        print(f"Error: {CONFIG_FILE} not found!")
        print(f"Please copy {CONFIG_FILE}.template to {CONFIG_FILE} and configure it.")
        sys.exit(1)

    with open(CONFIG_FILE) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                if "=" in line:
                    key, value = line.split("=", 1)
                    config[key.strip()] = value.strip()

    return config


def check_globus_cli():
    """Check if Globus CLI is installed and user is logged in."""
    try:
        result = subprocess.run(
            ["globus", "whoami"],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"Logged in as: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError:
        print("Error: Not logged in to Globus CLI")
        print("Please run: globus login")
        return False
    except FileNotFoundError:
        print("Error: Globus CLI not found")
        print("Please install with: pip install globus-cli")
        return False


def process_flow_definition(config: dict[str, str]) -> str:
    """Process flow definition with environment variables."""
    with open(FLOW_FILE) as f:
        flow_content = f.read()

    # Replace environment variables
    replacements = {
        "${ONEROOF_ACTION_URL}": config.get("ONEROOF_ACTION_URL", ""),
        "${ONEROOF_ACTION_SCOPE}": config.get("ONEROOF_ACTION_SCOPE", ""),
    }

    for key, value in replacements.items():
        if not value:
            print(f"Warning: {key} not set in config")
        flow_content = flow_content.replace(key, value)

    return flow_content


def register_flow(flow_definition: str, config: dict[str, str]) -> str | None:
    """Register the flow with Globus."""
    # Create temporary file with processed flow
    temp_flow = SCRIPT_DIR / "temp_flow.json"

    try:
        # Validate JSON
        flow_dict = json.loads(flow_definition)

        # Add metadata
        flow_dict["title"] = "OneRoof Bioinformatics Pipeline"
        flow_dict["description"] = (
            "Automated bioinformatics pipeline for base-calling, variant-calling, "
            "and consensus-calling of amplicon sequencing data. Supports both "
            "Nanopore and Illumina platforms."
        )

        with open(temp_flow, "w") as f:
            json.dump(flow_dict, f, indent=2)

        # Register flow
        print("Registering flow with Globus...")
        result = subprocess.run(
            ["globus", "flows", "create", "--definition", str(temp_flow)],
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            print(f"Error registering flow: {result.stderr}")
            return None

        # Extract flow ID from output
        output_lines = result.stdout.strip().split("\n")
        flow_id = None

        for line in output_lines:
            if "Flow ID:" in line:
                flow_id = line.split("Flow ID:")[1].strip()
                break

        if flow_id:
            print("\nFlow registered successfully!")
            print(f"Flow ID: {flow_id}")

            # Save flow ID to config
            with open(GLOBUS_DIR / "config" / "flow_id.txt", "w") as f:
                f.write(flow_id)

            return flow_id
        print("Warning: Could not extract flow ID from output")
        print(f"Output: {result.stdout}")

    finally:
        # Clean up temp file
        if temp_flow.exists():
            temp_flow.unlink()

    return None


def update_flow_permissions(flow_id: str, config: dict[str, str]):
    """Update flow permissions if needed."""
    print("\nUpdating flow permissions...")

    # Make flow runnable by all authenticated users
    subprocess.run(
        [
            "globus",
            "flows",
            "update",
            flow_id,
            "--runnable-by",
            "all_authenticated_users",
        ],
        check=False,
    )

    print("Permissions updated")


def print_usage_instructions(flow_id: str, config: dict[str, str]):
    """Print instructions for using the registered flow."""
    print("\n" + "=" * 60)
    print("Flow Registration Complete!")
    print("=" * 60)
    print(f"\nFlow ID: {flow_id}")
    print("\nNext steps:")
    print("\n1. Start the action provider on your compute endpoint:")
    print(f"   cd {GLOBUS_DIR}/scripts")
    print("   ./deploy.sh")

    print("\n2. Test the flow:")
    print(f"   python test_flow.py --flow-id {flow_id}")

    print("\n3. Run the flow with your data:")
    print(f"   globus flows run {flow_id} --input your_input.json")

    print("\nExample input file:")
    example_input = {
        "run_id": "test_run_001",
        "platform": "nanopore",
        "input_data": "/path/on/source/endpoint/pod5_files",
        "primer_bed": "/path/to/primers.bed",
        "refseq": "/path/to/reference.fasta",
        "ref_gbk": "/path/to/reference.gbk",
        "source_endpoint_id": config.get(
            "STORAGE_ENDPOINT_ID", "YOUR_SOURCE_ENDPOINT_ID"
        ),
        "compute_endpoint_id": config.get(
            "COMPUTE_ENDPOINT_ID", "YOUR_COMPUTE_ENDPOINT_ID"
        ),
        "destination_endpoint_id": config.get(
            "STORAGE_ENDPOINT_ID", "YOUR_DEST_ENDPOINT_ID"
        ),
        "output_path": "/path/on/destination/endpoint/results",
        "nextflow_params": {},
        "profile": "docker",
        "notification_recipients": [],
        "keep_logs": False,
    }

    print("\n" + json.dumps(example_input, indent=2))


def main():
    """Main registration process."""
    print("OneRoof Flow Registration")
    print("========================\n")

    # Load configuration
    config = load_env_config()

    # Check prerequisites
    if not check_globus_cli():
        sys.exit(1)

    # Verify required configuration
    required = ["ONEROOF_ACTION_URL", "ONEROOF_ACTION_SCOPE"]
    missing = [k for k in required if not config.get(k)]

    if missing:
        print(f"Error: Missing required configuration: {missing}")
        print(f"Please configure these in {CONFIG_FILE}")
        sys.exit(1)

    # Process and register flow
    flow_definition = process_flow_definition(config)
    flow_id = register_flow(flow_definition, config)

    if flow_id:
        update_flow_permissions(flow_id, config)
        print_usage_instructions(flow_id, config)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
