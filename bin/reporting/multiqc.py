"""
MultiQC custom content file generators for OneRoof reporting.

MultiQC supports custom content via specially formatted files with embedded
YAML configuration in comments. This module generates TSV files that MultiQC
will automatically detect and include in reports.

Reference: https://multiqc.info/docs/development/modules/#custom-content

File format:
    # id: 'oneroof_coverage'
    # section_name: 'Coverage Summary'
    # plot_type: 'table'
    # pconfig:
    #     namespace: 'OneRoof'
    # headers:
    #     col1:
    #         title: 'Column 1'
    #         format: '{:,.0f}'
    Sample	col1	col2
    sample1	100	200
"""

from pathlib import Path
from typing import Any


def _write_tsv_with_header(
    output_path: Path,
    yaml_config: dict[str, Any],
    headers: list[str],
    rows: list[list[Any]],
) -> None:
    """
    Write a TSV file with embedded YAML header for MultiQC.

    Args:
        output_path: Path to write the TSV file
        yaml_config: YAML configuration dict to embed in comments
        headers: Column headers
        rows: Data rows (list of lists)
    """
    lines = []

    # Write YAML config as comments
    for key, value in yaml_config.items():
        if isinstance(value, dict):
            lines.append(f"# {key}:")
            for subkey, subvalue in value.items():
                if isinstance(subvalue, dict):
                    lines.append(f"#     {subkey}:")
                    for subsubkey, subsubvalue in subvalue.items():
                        lines.append(f"#         {subsubkey}: '{subsubvalue}'")
                else:
                    lines.append(f"#     {subkey}: '{subvalue}'")
        else:
            lines.append(f"# {key}: '{value}'")

    # Write header row
    lines.append("\t".join(headers))

    # Write data rows
    for row in rows:
        lines.append("\t".join(str(v) for v in row))

    output_path.write_text("\n".join(lines) + "\n")


def generate_general_stats_tsv(
    samples: dict[str, dict],
    output_path: Path,
) -> Path:
    """
    Generate a TSV file for MultiQC General Statistics table.

    This adds OneRoof metrics to the main General Statistics table that
    appears at the top of every MultiQC report.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Path to write the TSV file

    Returns:
        Path to the written file
    """
    yaml_config = {
        "id": "oneroof_general_stats",
        "plot_type": "generalstats",
        "pconfig": {
            "namespace": "OneRoof",
        },
        "headers": {
            "mean_coverage": {
                "title": "Mean Cov",
                "description": "Mean coverage depth",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "genome_coverage_pct": {
                "title": "Genome %",
                "description": "Percentage of genome with ≥10x coverage",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
                "min": 0,
                "max": 100,
            },
            "completeness_pct": {
                "title": "Complete %",
                "description": "Consensus completeness (non-N bases)",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
                "min": 0,
                "max": 100,
            },
            "n_pct": {
                "title": "N %",
                "description": "Percentage of N bases in consensus",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "OrRd",
                "min": 0,
                "max": 100,
            },
        },
    }

    headers = [
        "Sample",
        "mean_coverage",
        "genome_coverage_pct",
        "completeness_pct",
        "n_pct",
    ]
    rows = []

    for sample_id, metrics in samples.items():
        # Extract coverage metrics
        coverage = metrics.get("coverage", {})
        consensus = metrics.get("consensus", {})

        mean_cov = coverage.get("mean_coverage", 0)
        genome_pct = coverage.get("genome_coverage_at_10x", 0) * 100
        completeness = consensus.get("completeness", 0) * 100 if consensus else 0
        n_pct = consensus.get("n_percentage", 0) if consensus else 0

        rows.append([sample_id, mean_cov, genome_pct, completeness, n_pct])

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


def generate_coverage_table_tsv(
    samples: dict[str, dict],
    output_path: Path,
) -> Path:
    """
    Generate a TSV file for a detailed coverage table section.

    Args:
        samples: Dict mapping sample_id to metrics dict
        output_path: Path to write the TSV file

    Returns:
        Path to the written file
    """
    yaml_config = {
        "id": "oneroof_coverage_table",
        "section_name": "Coverage Summary",
        "description": "Detailed coverage statistics from OneRoof analysis",
        "plot_type": "table",
        "pconfig": {
            "namespace": "OneRoof",
            "id": "oneroof_coverage_table",
            "title": "Coverage Statistics",
        },
        "headers": {
            "mean_coverage": {
                "title": "Mean Coverage",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "median_coverage": {
                "title": "Median Coverage",
                "format": "{:,.1f}",
                "scale": "Blues",
            },
            "genome_1x": {
                "title": "≥1x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "genome_10x": {
                "title": "≥10x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "genome_100x": {
                "title": "≥100x %",
                "format": "{:,.1f}",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "min_coverage": {
                "title": "Min",
                "format": "{:,.0f}",
            },
            "max_coverage": {
                "title": "Max",
                "format": "{:,.0f}",
            },
        },
    }

    headers = [
        "Sample",
        "mean_coverage",
        "median_coverage",
        "genome_1x",
        "genome_10x",
        "genome_100x",
        "min_coverage",
        "max_coverage",
    ]
    rows = []

    for sample_id, metrics in samples.items():
        coverage = metrics.get("coverage", {})

        rows.append(
            [
                sample_id,
                coverage.get("mean_coverage", 0),
                coverage.get("median_coverage", 0),
                coverage.get("genome_coverage_at_1x", 0) * 100,
                coverage.get("genome_coverage_at_10x", 0) * 100,
                coverage.get("genome_coverage_at_100x", 0) * 100,
                coverage.get("min_coverage", 0),
                coverage.get("max_coverage", 0),
            ]
        )

    _write_tsv_with_header(output_path, yaml_config, headers, rows)
    return output_path


# Default relative path to the MultiQC config template from project root
DEFAULT_MULTIQC_TEMPLATE = Path("assets/multiqc_config.yaml")


def generate_multiqc_config(
    run_metadata: dict,
    output_path: Path,
    template_path: Path | None = None,
) -> Path:
    """
    Generate a MultiQC configuration file from template.

    Loads the base template and substitutes dynamic values like platform
    and reference name.

    Args:
        run_metadata: Run metadata dict with platform, reference info, etc.
        output_path: Path to write the YAML config file
        template_path: Path to template file. If None, uses DEFAULT_MULTIQC_TEMPLATE
                      relative to current working directory.

    Returns:
        Path to the written file

    Raises:
        FileNotFoundError: If template_path does not exist
    """
    if template_path is None:
        template_path = DEFAULT_MULTIQC_TEMPLATE

    if not template_path.exists():
        msg = f"MultiQC config template not found: {template_path}"
        raise FileNotFoundError(msg)

    platform = run_metadata.get("platform", "unknown").upper()
    reference_name = run_metadata.get("reference", {}).get("name", "unknown")

    # Load template and substitute placeholders
    template_content = template_path.read_text()
    config_content = template_content.format(
        platform=platform,
        reference=reference_name,
    )

    output_path.write_text(config_content)
    return output_path
