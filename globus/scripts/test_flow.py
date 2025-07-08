#!/usr/bin/env python3
"""
Test OneRoof Globus Flow with sample data.

This script helps test the registered flow with various configurations.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).parent
GLOBUS_DIR = SCRIPT_DIR.parent
CONFIG_FILE = GLOBUS_DIR / "config" / ".env"
FLOW_ID_FILE = GLOBUS_DIR / "config" / "flow_id.txt"


def load_config() -> dict[str, str]:
    """Load configuration from .env file."""
    config = {}

    if CONFIG_FILE.exists():
        with open(CONFIG_FILE) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    if "=" in line:
                        key, value = line.split("=", 1)
                        config[key.strip()] = value.strip()

    return config


def get_flow_id(args_flow_id: str | None = None) -> str:
    """Get flow ID from arguments or saved file."""
    if args_flow_id:
        return args_flow_id

    if FLOW_ID_FILE.exists():
        with open(FLOW_ID_FILE) as f:
            return f.read().strip()

    print("Error: No flow ID provided and no saved flow ID found")
    print("Please provide --flow-id or run register_flow.py first")
    sys.exit(1)


def create_test_input(
    platform: str,
    input_data: str,
    primer_bed: str,
    refseq: str,
    ref_gbk: str,
    output_path: str,
    config: dict[str, str],
    extra_params: dict[str, Any],
) -> dict[str, Any]:
    """Create test input for the flow."""
    run_id = f"test_{platform}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    flow_input = {
        "run_id": run_id,
        "platform": platform,
        "input_data": input_data,
        "primer_bed": primer_bed,
        "refseq": refseq,
        "ref_gbk": ref_gbk,
        "source_endpoint_id": config.get("STORAGE_ENDPOINT_ID", ""),
        "compute_endpoint_id": config.get("COMPUTE_ENDPOINT_ID", ""),
        "destination_endpoint_id": config.get("STORAGE_ENDPOINT_ID", ""),
        "work_dir": f"{config.get('ONEROOF_WORK_DIR', '/tmp/oneroof_runs')}/{run_id}",
        "output_path": output_path,
        "nextflow_params": {},
        "profile": "docker",
        "notification_recipients": [],
        "keep_logs": True,
    }

    # Add extra parameters
    if extra_params.get("profile"):
        flow_input["profile"] = extra_params["profile"]

    if extra_params.get("notifications"):
        flow_input["notification_recipients"] = extra_params["notifications"]

    if extra_params.get("nextflow_params"):
        flow_input["nextflow_params"] = extra_params["nextflow_params"]

    # Override with any endpoint IDs from command line
    if extra_params.get("source_endpoint"):
        flow_input["source_endpoint_id"] = extra_params["source_endpoint"]
    if extra_params.get("compute_endpoint"):
        flow_input["compute_endpoint_id"] = extra_params["compute_endpoint"]
    if extra_params.get("dest_endpoint"):
        flow_input["destination_endpoint_id"] = extra_params["dest_endpoint"]

    return flow_input


def run_flow(flow_id: str, flow_input: dict[str, Any]) -> str | None:
    """Submit the flow run."""
    # Save input to temporary file
    temp_input = SCRIPT_DIR / "temp_input.json"

    try:
        with open(temp_input, "w") as f:
            json.dump(flow_input, f, indent=2)

        print("\nSubmitting flow run...")
        print(f"Flow ID: {flow_id}")
        print(f"Run ID: {flow_input['run_id']}")

        # Submit flow
        result = subprocess.run(
            ["globus", "flows", "run", flow_id, "--input", str(temp_input)],
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            print(f"Error submitting flow: {result.stderr}")
            return None

        # Extract run ID from output
        output_lines = result.stdout.strip().split("\n")
        run_id = None

        for line in output_lines:
            if "Run ID:" in line:
                run_id = line.split("Run ID:")[1].strip()
                break

        if run_id:
            print("\nFlow submitted successfully!")
            print(f"Flow Run ID: {run_id}")
            return run_id
        print("Warning: Could not extract run ID from output")
        print(f"Output: {result.stdout}")

    finally:
        # Clean up temp file
        if temp_input.exists():
            temp_input.unlink()

    return None


def monitor_flow(flow_id: str, run_id: str):
    """Monitor the flow run status."""
    print("\nMonitoring flow run...")
    print("To check status manually:")
    print(f"  globus flows run-status {flow_id} {run_id}")
    print("\nTo view logs:")
    print(f"  globus flows run-logs {flow_id} {run_id}")

    # Initial status check
    subprocess.run(
        [
            "globus",
            "flows",
            "run-status",
            flow_id,
            run_id,
        ],
        check=False,
    )


def main():
    """Main test function."""
    parser = argparse.ArgumentParser(
        description="Test OneRoof Globus Flow",
    )

    parser.add_argument(
        "--flow-id",
        help="Flow ID (if not provided, will read from saved file)",
    )

    parser.add_argument(
        "--platform",
        choices=["nanopore", "illumina"],
        required=True,
        help="Sequencing platform",
    )

    parser.add_argument(
        "--input-data",
        required=True,
        help="Path to input data on source endpoint",
    )

    parser.add_argument(
        "--primer-bed",
        required=True,
        help="Path to primer BED file",
    )

    parser.add_argument(
        "--refseq",
        required=True,
        help="Path to reference FASTA",
    )

    parser.add_argument(
        "--ref-gbk",
        required=True,
        help="Path to reference GenBank file",
    )

    parser.add_argument(
        "--output-path",
        required=True,
        help="Output path on destination endpoint",
    )

    parser.add_argument(
        "--profile",
        choices=["docker", "singularity", "containerless"],
        default="docker",
        help="Nextflow profile to use",
    )

    parser.add_argument(
        "--source-endpoint",
        help="Override source endpoint ID",
    )

    parser.add_argument(
        "--compute-endpoint",
        help="Override compute endpoint ID",
    )

    parser.add_argument(
        "--dest-endpoint",
        help="Override destination endpoint ID",
    )

    parser.add_argument(
        "--notifications",
        nargs="+",
        help="Email addresses for notifications",
    )

    parser.add_argument(
        "--nextflow-param",
        action="append",
        help="Additional Nextflow parameters (format: key=value)",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print input JSON without submitting",
    )

    args = parser.parse_args()

    # Load configuration
    config = load_config()

    # Get flow ID
    flow_id = get_flow_id(args.flow_id)

    # Process extra parameters
    extra_params = {
        "profile": args.profile,
        "notifications": args.notifications or [],
        "nextflow_params": {},
    }

    if args.source_endpoint:
        extra_params["source_endpoint"] = args.source_endpoint
    if args.compute_endpoint:
        extra_params["compute_endpoint"] = args.compute_endpoint
    if args.dest_endpoint:
        extra_params["dest_endpoint"] = args.dest_endpoint

    # Parse nextflow parameters
    if args.nextflow_param:
        for param in args.nextflow_param:
            if "=" in param:
                key, value = param.split("=", 1)
                # Try to parse as number or boolean
                if value.lower() in ["true", "false"]:
                    value = value.lower() == "true"
                else:
                    try:
                        value = int(value)
                    except ValueError:
                        try:
                            value = float(value)
                        except ValueError:
                            pass  # Keep as string

                extra_params["nextflow_params"][key] = value

    # Create flow input
    flow_input = create_test_input(
        args.platform,
        args.input_data,
        args.primer_bed,
        args.refseq,
        args.ref_gbk,
        args.output_path,
        config,
        extra_params,
    )

    # Validate endpoint IDs
    endpoint_ids = [
        flow_input["source_endpoint_id"],
        flow_input["compute_endpoint_id"],
        flow_input["destination_endpoint_id"],
    ]

    if not all(endpoint_ids):
        print("Error: Missing endpoint IDs")
        print("Please configure endpoint IDs in .env or provide via command line")
        print(f"  Source: {flow_input['source_endpoint_id'] or 'NOT SET'}")
        print(f"  Compute: {flow_input['compute_endpoint_id'] or 'NOT SET'}")
        print(f"  Destination: {flow_input['destination_endpoint_id'] or 'NOT SET'}")
        sys.exit(1)

    # Print input
    print("\nFlow Input:")
    print(json.dumps(flow_input, indent=2))

    if args.dry_run:
        print("\nDry run - not submitting flow")
        return

    # Run flow
    run_id = run_flow(flow_id, flow_input)

    if run_id:
        monitor_flow(flow_id, run_id)


if __name__ == "__main__":
    main()
