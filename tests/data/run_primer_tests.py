#!/usr/bin/env python3
"""
Simple test runner for primer-related tests.
This script helps verify that all test modules can be imported and discovered by pytest.
"""

import subprocess
import sys
from pathlib import Path


def main():
    """Run all primer-related tests"""
    test_dir = Path(__file__).parent

    # List of test modules
    test_modules = [
        "test_make_primer_patterns.py",
        "test_resplice_primers.py",
        "test_split_primer_combos.py",
        "test_validate_primer_bed.py",
    ]

    print("Running primer-related tests...")
    print("=" * 60)

    # Run pytest on each module
    for module in test_modules:
        test_file = test_dir / module
        if test_file.exists():
            print(f"\nRunning tests in {module}...")
            cmd = [sys.executable, "-m", "pytest", str(test_file), "-v", "--tb=short"]
            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                print(f"✓ {module} passed")
            else:
                print(f"✗ {module} failed")
                print(result.stdout)
                print(result.stderr)
        else:
            print(f"✗ {module} not found")

    print("\n" + "=" * 60)
    print("Test run complete!")


if __name__ == "__main__":
    main()
