#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# ///

"""
fix-index-paths.py
Post-render script to fix paths in the generated site
This script:
1. Copies docs/index.html to the site root and fixes relative paths
2. Copies markdown files from _site/docs back to docs directory
"""

import re
import shutil
from pathlib import Path


def main():
    # Define paths relative to project root
    site_dir = Path("_site")
    docs_dir = Path("docs")
    source_index = site_dir / "docs" / "index.html"
    dest_index = site_dir / "index.html"

    # Part 1: Copy and fix index.html
    if source_index.exists():
        try:
            # Read the content
            content = source_index.read_text(encoding="utf-8")

            # Fix relative paths
            content = re.sub(r'href="../docs/', 'href="docs/', content)
            content = re.sub(r'src="../site_libs/', 'src="site_libs/', content)

            # Write the fixed content
            dest_index.write_text(content, encoding="utf-8")
            print("✓ Fixed paths in root index.html")

        except Exception as e:
            print(f"⚠️  Failed to process index.html: {e}")
    else:
        print("⚠️  No _site directory found - using default homepage")

    # Part 2: Copy markdown files back to docs directory
    site_docs_dir = site_dir / "docs"
    if site_docs_dir.is_dir():
        copied = False
        try:
            for md_file in site_docs_dir.glob("*.md"):
                dest_file = docs_dir / md_file.name
                shutil.copy2(md_file, dest_file)
                copied = True

            if copied:
                print("✓ Markdown files copied to docs/")
            else:
                print("⚠️  No markdown files found")

        except Exception as e:
            print(f"⚠️  Failed to copy markdown files: {e}")


if __name__ == "__main__":
    main()
