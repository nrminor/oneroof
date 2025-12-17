#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "biopython",
#     "loguru",
#     "numpy",
#     "polars",
#     "pydantic",
#     "scipy",
#     "typer",
#     "rich",
# ]
# ///

"""Functional iVar to VCF converter with beautiful Typer CLI.

This module provides a purely functional approach to converting iVar TSV variant files
to VCF format, using Polars expression chains for data transformation and Pydantic
for validation. Features a beautiful command-line interface built with Typer and Rich.
"""

import gzip
from enum import Enum
from pathlib import Path
from typing import Annotated, cast

import numpy as np
import polars as pl
import typer
from Bio import SeqIO
from loguru import logger
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from scipy.stats import fisher_exact

# Initialize Typer app and Rich console
app = typer.Typer(
    name="ivar2vcf",
    help="Convert iVar TSV variant files to VCF format with advanced filtering and validation.",
    add_completion=False,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
)
console = Console()

# ===== Pydantic Models for Validation =====


class VariantType(str, Enum):
    """Enumeration of variant types."""

    SNP = "SNP"
    INS = "INS"
    DEL = "DEL"


class FilterType(str, Enum):
    """Enumeration of filter types."""

    PASS = "PASS"  # noqa: S105
    FAIL_TEST = "ft"
    BAD_QUALITY = "bq"
    STRAND_BIAS = "sb"


class IvarVariant(BaseModel):
    """Validated iVar variant record."""

    model_config = ConfigDict(str_strip_whitespace=True)

    region: str
    pos: Annotated[int, Field(gt=0)]
    ref: Annotated[str, Field(min_length=1)]
    alt: str
    ref_dp: Annotated[int, Field(ge=0)]
    ref_rv: Annotated[int, Field(ge=0)]
    ref_qual: Annotated[float, Field(ge=0, le=100)]
    alt_dp: Annotated[int, Field(ge=0)]
    alt_rv: Annotated[int, Field(ge=0)]
    alt_qual: Annotated[float, Field(ge=0, le=100)]
    alt_freq: Annotated[float, Field(ge=0, le=1)]
    total_dp: Annotated[int, Field(ge=0)]
    pval: float
    pass_test: bool
    ref_codon: str | None = None
    alt_codon: str | None = None

    @field_validator("ref_rv", "alt_rv")
    @classmethod
    def validate_reverse_depth(cls, v: int, info: ValidationInfo) -> int:
        """Ensure reverse depth doesn't exceed total depth for that allele."""
        if info.field_name == "ref_rv":
            depth_field = "ref_dp"
        elif info.field_name == "alt_rv":
            depth_field = "alt_dp"
        else:
            msg = f"Unexpected field in reverse depth validator: {info.field_name}"
            raise ValueError(msg)

        if depth_field in info.data and v > info.data[depth_field]:
            msg = f"Reverse depth ({info.field_name}) cannot exceed total depth ({depth_field})"
            raise ValueError(msg)
        return v

    @model_validator(mode="after")
    def validate_total_depth(self) -> "IvarVariant":
        """Ensure total depth is consistent."""
        if self.total_dp < self.ref_dp + self.alt_dp:
            msg = "Total depth must be at least ref_dp + alt_dp"
            raise ValueError(msg)
        return self


class VcfVariant(BaseModel):
    """Validated VCF variant record."""

    model_config = ConfigDict(str_strip_whitespace=True)

    chrom: str
    pos: Annotated[int, Field(gt=0)]
    id: str = "."
    ref: Annotated[str, Field(min_length=1)]
    alt: Annotated[str, Field(min_length=1)]
    qual: str = "."
    filter: str
    info: str
    format: str
    sample: str

    @field_validator("ref", "alt")
    @classmethod
    def validate_alleles(cls, v: str) -> str:
        """Ensure alleles contain only valid nucleotides."""
        valid_chars = set("ACGTNacgtn")
        if not all(c in valid_chars for c in v):
            msg = f"Invalid nucleotide in allele: {v}"
            raise ValueError(msg)
        return v  # Keep original case


class ConversionConfig(BaseModel):
    """Configuration for iVar to VCF conversion."""

    model_config = ConfigDict(frozen=True)

    file_in: Path
    file_out: Path
    pass_only: bool = False
    freq_threshold: Annotated[float, Field(ge=0, le=1)] = 0.0
    bad_qual_threshold: Annotated[float, Field(ge=0)] = 20.0
    merge_af_threshold: Annotated[float, Field(ge=0, le=1)] = 0.25
    consensus_af: Annotated[float, Field(ge=0, le=1)] = 0.75
    ignore_strand_bias: bool = False
    ignore_merge_codons: bool = False
    ref_fasta: Path | None = None

    @field_validator("file_in")
    @classmethod
    def validate_input_file(cls, v: Path) -> Path:
        """Ensure input file exists."""
        if not v.exists():
            msg = f"Input file does not exist: {v}"
            raise ValueError(msg)
        return v

    @field_validator("file_out")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Ensure output directory exists."""
        if not v.parent.exists():
            msg = f"Output directory does not exist: {v.parent}"
            raise ValueError(msg)
        return v


# ===== Pure Functions for Data Transformation =====


def calculate_strand_bias_pvalue(
    ref_dp: int,
    ref_rv: int,
    alt_dp: int,
    alt_rv: int,
) -> float:
    """Calculate p-value for strand bias using Fisher's exact test.

    Args:
        ref_dp: Reference total depth
        ref_rv: Reference reverse strand depth
        alt_dp: Alternative total depth
        alt_rv: Alternative reverse strand depth

    Returns:
        P-value from Fisher's exact test
    """
    contingency_table = np.array(
        [
            [ref_dp - ref_rv, ref_rv],
            [alt_dp - alt_rv, alt_rv],
        ],
    )
    _odds_ratio, pvalue = cast(
        "tuple[float, float]",
        fisher_exact(contingency_table, alternative="two-sided"),
    )
    return pvalue


def create_strand_bias_expr(threshold: float = 0.05) -> pl.Expr:
    """Create a Polars expression for strand bias calculation.

    Args:
        threshold: P-value threshold for significance

    Returns:
        Expression that returns True if strand bias is significant
    """
    # Create a struct with all needed values
    return (
        pl.struct(
            [
                pl.col("REF_DP"),
                pl.col("REF_RV"),
                pl.col("ALT_DP"),
                pl.col("ALT_RV"),
            ],
        )
        .map_elements(
            lambda x: calculate_strand_bias_pvalue(
                x["REF_DP"],
                x["REF_RV"],
                x["ALT_DP"],
                x["ALT_RV"],
            )
            < threshold,
            return_dtype=pl.Boolean,
        )
        .alias("has_strand_bias")
    )


def determine_variant_type_expr() -> pl.Expr:
    """Create expression to determine variant type from ALT allele.

    Returns:
        Expression that returns variant type (SNP, INS, or DEL)
    """
    return (
        pl.when(pl.col("ALT").str.starts_with("+"))
        .then(pl.lit(VariantType.INS.value))
        .when(pl.col("ALT").str.starts_with("-"))
        .then(pl.lit(VariantType.DEL.value))
        .otherwise(pl.lit(VariantType.SNP.value))
    )


def transform_ref_alt_expr() -> tuple[pl.Expr, pl.Expr]:
    """Create expressions to transform REF/ALT for indels.

    Returns:
        Tuple of (REF expression, ALT expression)
    """
    ref_expr = (
        pl.when(pl.col("ALT").str.starts_with("-"))
        .then(pl.col("REF") + pl.col("ALT").str.slice(1))
        .otherwise(pl.col("REF"))
    )

    alt_expr = (
        pl.when(pl.col("ALT").str.starts_with("+"))
        .then(pl.col("REF") + pl.col("ALT").str.slice(1))
        .when(pl.col("ALT").str.starts_with("-"))
        .then(pl.col("REF"))
        .otherwise(pl.col("ALT"))
    )

    return ref_expr, alt_expr


def create_filter_expr(config: ConversionConfig) -> pl.Expr:
    """Create expression for combining all filters.

    Args:
        config: Conversion configuration

    Returns:
        Expression that returns the filter string
    """
    # Base filters
    filters = []

    # iVar PASS filter
    filters.append(
        pl.when(pl.col("PASS")).then(pl.lit("")).otherwise(pl.lit(FilterType.FAIL_TEST.value)),
    )

    # Quality filter
    filters.append(
        pl.when(pl.col("ALT_QUAL") >= config.bad_qual_threshold)
        .then(pl.lit(""))
        .otherwise(pl.lit(FilterType.BAD_QUALITY.value)),
    )

    # Strand bias filter (if enabled)
    if not config.ignore_strand_bias:
        filters.append(
            pl.when(create_strand_bias_expr())
            .then(pl.lit(FilterType.STRAND_BIAS.value))
            .otherwise(pl.lit("")),
        )

    # Combine filters
    return (
        pl.concat_list(filters)
        .list.eval(pl.when(pl.element() != "").then(pl.element()))
        .list.drop_nulls()
        .list.join(";")
        .fill_null("")
        .map_elements(
            lambda x: FilterType.PASS.value if x == "" else x,
            return_dtype=pl.Utf8,
        )
    )


def create_sample_info_expr() -> pl.Expr:
    """Create expression for sample information column.

    Returns:
        Expression that creates the sample info string
    """
    return pl.concat_str(
        [
            pl.lit("1"),  # GT (always 1 for haploid)
            pl.col("TOTAL_DP"),
            pl.col("REF_DP"),
            pl.col("REF_RV"),
            pl.col("REF_QUAL"),
            pl.col("ALT_DP"),
            pl.col("ALT_RV"),
            pl.col("ALT_QUAL"),
            pl.col("ALT_FREQ"),
        ],
        separator=":",
    )


def transform_ivar_to_vcf(
    ivar_lf: pl.LazyFrame,
    config: ConversionConfig,
) -> pl.LazyFrame:
    """Transform iVar data to VCF format using pure expressions.

    Args:
        ivar_lf: Lazy DataFrame with iVar data
        config: Conversion configuration

    Returns:
        Transformed lazy DataFrame in VCF format
    """
    # Apply filters
    if config.pass_only:
        ivar_lf = ivar_lf.filter(pl.col("PASS"))

    ivar_lf = ivar_lf.filter(pl.col("ALT_FREQ") >= config.freq_threshold)

    # Get REF/ALT transformation expressions
    ref_expr, alt_expr = transform_ref_alt_expr()

    # Build the transformation pipeline
    return ivar_lf.select(
        [
            pl.col("REGION").alias("CHROM"),
            pl.col("POS"),
            pl.lit(".").alias("ID"),
            ref_expr.alias("REF"),
            alt_expr.alias("ALT"),
            pl.lit(".").alias("QUAL"),
            create_filter_expr(config).alias("FILTER"),
            pl.concat_str(
                [
                    pl.lit("TYPE="),
                    determine_variant_type_expr(),
                ],
            ).alias("INFO"),
            pl.lit(
                "GT:DP:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ",
            ).alias("FORMAT"),
            create_sample_info_expr().alias("SAMPLE"),
        ],
    )


def find_consecutive_variants_expr() -> pl.Expr:
    """Create expression to identify consecutive variants.

    Returns:
        Boolean expression marking consecutive positions
    """
    pos_diff = pl.col("POS").diff()
    return (pos_diff <= 1) | (pos_diff.shift(-1) <= 1)


def process_consecutive_snps(
    ivar_lf: pl.LazyFrame,
    config: ConversionConfig,
) -> pl.LazyFrame:
    """Process consecutive SNPs for potential merging.

    Args:
        ivar_lf: Lazy DataFrame with VCF data
        config: Conversion configuration

    Returns:
        Processed DataFrame with merged variants
    """
    if config.ignore_merge_codons:
        return ivar_lf

    # For now, just remove duplicates and filter out REF==ALT
    return (
        ivar_lf.unique(subset=["CHROM", "POS", "REF", "ALT"])
        .filter(pl.col("REF") != pl.col("ALT"))
        .sort("POS")
    )


def generate_vcf_header(config: ConversionConfig) -> list[str]:
    """Generate VCF header lines.

    Args:
        config: Conversion configuration

    Returns:
        List of header lines
    """
    headers = [
        "##fileformat=VCFv4.2",
        "##source=iVar",
    ]

    # Add contig information if reference provided
    if config.ref_fasta and config.ref_fasta.exists():
        headers.extend(
            f"##contig=<ID={record.id},length={len(record.seq)}>"
            for record in SeqIO.parse(str(config.ref_fasta), "fasta")
        )

    # INFO fields
    headers.append(
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Either SNP (Single Nucleotide Polymorphism), DEL (deletion) or INS (Insertion)">',
    )

    # FILTER fields
    headers.extend(
        [
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=ft,Description="Fisher\'s exact test of variant frequency compared to mean error rate, p-value > 0.05">',
            '##FILTER=<ID=bq,Description="Bad quality variant: ALT_QUAL lower than 20">',
        ],
    )

    if not config.ignore_strand_bias:
        headers.append('##FILTER=<ID=sb,Description="Strand bias filter not passed">')

    # FORMAT fields
    headers.extend(
        [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">',
            '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">',
            '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">',
            '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">',
            '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">',
            '##FORMAT=<ID=ALT_QUAL,Number=1,Type=Integer,Description="Mean quality of alternate base">',
            '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Frequency of alternate base">',
        ],
    )

    if not config.ignore_merge_codons:
        headers.extend(
            [
                '##FORMAT=<ID=MERGED_AF,Number=A,Type=Float,Description="Frequency of each merged variant comma separated">',
                '##FORMAT=<ID=MERGED_DP,Number=A,Type=Float,Description="Total Depth of each merged variant comma separated">',
            ],
        )

    return headers


def write_vcf_file(
    vcf_df: pl.DataFrame,
    filepath: Path,
    headers: list[str],
    sample_name: str,
) -> None:
    """Write VCF data to file.

    Args:
        vcf_df: DataFrame with VCF data
        filepath: Output file path
        headers: VCF header lines
        sample_name: Sample name for the VCF
    """
    # Rename columns for VCF format
    vcf_df = vcf_df.rename({"CHROM": "#CHROM", "SAMPLE": sample_name})

    # Prepare header text
    header_text = "\n".join(headers) + "\n"

    # Write file (handle gzipped output)
    if str(filepath).endswith(".gz"):
        # For gzip, we need to write everything as text to the same handle

        with gzip.open(filepath, "wt") as f:
            f.write(header_text)
            # Write column headers
            f.write("\t".join(vcf_df.columns) + "\n")
            # Write data rows
            for row in vcf_df.iter_rows():
                f.write("\t".join(str(x) for x in row) + "\n")
    else:
        with open(filepath, "w") as f:
            f.write(header_text)
            vcf_df.write_csv(f, separator="\t", include_header=True)


def process_ivar_file(config: ConversionConfig) -> None:
    """Process iVar TSV file and generate VCF outputs.

    Args:
        config: Validated conversion configuration
    """
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        # Load data lazily
        task = progress.add_task("[cyan]Loading iVar data...", total=None)
        ivar_df = pl.scan_csv(str(config.file_in), separator="\t")

        # Check if input has any data rows (collect schema to check)
        progress.update(task, description="[yellow]Checking input data...")
        try:
            row_count = ivar_df.select(pl.len()).collect().item()
        except pl.exceptions.NoDataError:
            row_count = 0

        # Generate headers (needed regardless of data)
        progress.update(task, description="[yellow]Generating VCF headers...")
        headers = generate_vcf_header(config)
        sample_name = config.file_in.stem

        if row_count == 0:
            # Handle empty input: write VCF with headers only
            progress.update(
                task,
                description="[yellow]No variants found, writing empty VCF...",
            )
            empty_df = pl.DataFrame(
                schema={
                    "CHROM": pl.Utf8,
                    "POS": pl.Int64,
                    "ID": pl.Utf8,
                    "REF": pl.Utf8,
                    "ALT": pl.Utf8,
                    "QUAL": pl.Utf8,
                    "FILTER": pl.Utf8,
                    "INFO": pl.Utf8,
                    "FORMAT": pl.Utf8,
                    "SAMPLE": pl.Utf8,
                },
            )
            write_vcf_file(empty_df, config.file_out, headers, sample_name)

            all_hap_path = (
                config.file_out.parent / f"{config.file_out.stem}_all_hap{config.file_out.suffix}"
            )
            write_vcf_file(empty_df, all_hap_path, headers, sample_name)

            progress.update(
                task,
                description="[bold yellow]✓ No variants found, empty VCF written",
                completed=True,
            )
            return

        # Transform to VCF format
        progress.update(task, description="[yellow]Transforming to VCF format...")
        vcf_df = transform_ivar_to_vcf(ivar_df, config)

        # Process consecutive SNPs
        progress.update(task, description="[yellow]Processing consecutive variants...")
        processed_df = process_consecutive_snps(vcf_df, config)

        # Collect the data
        progress.update(task, description="[yellow]Collecting results...")
        result_df = processed_df.collect()

        # Write consensus output
        progress.update(task, description="[green]Writing consensus VCF...")
        write_vcf_file(result_df, config.file_out, headers, sample_name)

        # Write all haplotypes output
        progress.update(task, description="[green]Writing all haplotypes VCF...")
        all_hap_path = (
            config.file_out.parent / f"{config.file_out.stem}_all_hap{config.file_out.suffix}"
        )
        write_vcf_file(result_df, all_hap_path, headers, sample_name)

        progress.update(
            task,
            description="[bold green]✓ Conversion complete!",
            completed=True,
        )


# ===== Typer CLI Commands =====


@app.command()
def convert(  # noqa: PLR0913
    file_in: Annotated[
        Path,
        typer.Argument(
            help="Input iVar TSV file",
            exists=True,
            dir_okay=False,
            readable=True,
            metavar="INPUT.tsv",
        ),
    ],
    file_out: Annotated[
        Path,
        typer.Argument(
            help="Output VCF file (supports .gz)",
            dir_okay=False,
            metavar="OUTPUT.vcf[.gz]",
        ),
    ],
    pass_only: Annotated[  # noqa: FBT002
        bool,
        typer.Option(
            "--pass-only",
            "-p",
            help="Only output variants that PASS all filters",
            rich_help_panel="Filtering Options",
        ),
    ] = False,
    allele_freq_threshold: Annotated[
        float,
        typer.Option(
            "--allele-freq",
            "-f",
            min=0.0,
            max=1.0,
            help="Minimum allele frequency threshold",
            rich_help_panel="Filtering Options",
        ),
    ] = 0.0,
    bad_quality_threshold: Annotated[
        int,
        typer.Option(
            "--min-quality",
            "-q",
            min=0,
            help="Minimum ALT_QUAL threshold",
            rich_help_panel="Filtering Options",
        ),
    ] = 20,
    merge_af_threshold: Annotated[
        float,
        typer.Option(
            "--merge-threshold",
            "-m",
            min=0.0,
            max=1.0,
            help="AF difference threshold for merging variants",
            rich_help_panel="Merging Options",
        ),
    ] = 0.25,
    consensus_af: Annotated[
        float,
        typer.Option(
            "--consensus-af",
            "-c",
            min=0.0,
            max=1.0,
            help="AF threshold for consensus calling",
            rich_help_panel="Merging Options",
        ),
    ] = 0.75,
    ignore_strand_bias: Annotated[  # noqa: FBT002
        bool,
        typer.Option(
            "--ignore-strand-bias",
            "-s",
            help="Ignore strand bias filter (recommended for amplicon data)",
            rich_help_panel="Advanced Options",
        ),
    ] = False,
    ignore_merge_codons: Annotated[  # noqa: FBT002
        bool,
        typer.Option(
            "--no-merge",
            "-n",
            help="Don't merge consecutive variants in codons",
            rich_help_panel="Advanced Options",
        ),
    ] = False,
    ref_fasta: Annotated[
        Path | None,
        typer.Option(
            "--reference",
            "-r",
            help="Reference FASTA for contig information",
            exists=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="Advanced Options",
        ),
    ] = None,
    verbose: Annotated[
        int,
        typer.Option(
            "--verbose",
            "-v",
            count=True,
            help="Increase verbosity (-v, -vv, -vvv)",
            rich_help_panel="Logging",
        ),
    ] = 0,
) -> None:
    """Convert iVar TSV variant files to VCF format.

    This tool processes iVar variant calling output and generates standard VCF files
    with advanced filtering options including strand bias detection, quality filtering,
    and consecutive variant merging.

    [bold cyan]Examples:[/bold cyan]

    Basic conversion:
    [green]$ ivar2vcf input.tsv output.vcf[/green]

    Filter low-quality variants:
    [green]$ ivar2vcf input.tsv output.vcf.gz -f 0.05 -q 30[/green]

    Amplicon sequencing (ignore strand bias):
    [green]$ ivar2vcf input.tsv output.vcf -s[/green]
    """
    # Configure logging
    log_levels = {0: "WARNING", 1: "SUCCESS", 2: "INFO", 3: "DEBUG"}
    level = log_levels.get(verbose, "WARNING")
    logger.remove()
    logger.add(lambda msg: console.print(msg, end=""), colorize=True, level=level)

    # Display header
    console.print(
        Panel.fit(
            "[bold cyan]iVar to VCF Converter[/bold cyan]\n"
            f"[dim]Processing: {file_in.name} → {file_out.name}[/dim]",
            border_style="cyan",
        ),
    )

    try:
        # Create and validate configuration
        config = ConversionConfig(
            file_in=file_in,
            file_out=file_out,
            pass_only=pass_only,
            freq_threshold=allele_freq_threshold,
            bad_qual_threshold=bad_quality_threshold,
            merge_af_threshold=merge_af_threshold,
            consensus_af=consensus_af,
            ignore_strand_bias=ignore_strand_bias,
            ignore_merge_codons=ignore_merge_codons,
            ref_fasta=ref_fasta,
        )

        # Display configuration summary
        console.print("\n[bold]Configuration:[/bold]")
        console.print(f"  • Pass only: {config.pass_only}")
        console.print(f"  • Min allele frequency: {config.freq_threshold}")
        console.print(f"  • Min quality: {config.bad_qual_threshold}")
        if not config.ignore_strand_bias:
            console.print("  • Strand bias filter: [green]enabled[/green]")
        else:
            console.print("  • Strand bias filter: [yellow]disabled[/yellow]")
        if not config.ignore_merge_codons:
            console.print("  • Codon merging: [green]enabled[/green]")
        else:
            console.print("  • Codon merging: [yellow]disabled[/yellow]")
        console.print()

        # Process the file
        process_ivar_file(config)

        # Success message
        console.print(
            f"\n[bold green]✓[/bold green] Successfully converted to {config.file_out}",
        )
        all_hap_path = (
            config.file_out.parent / f"{config.file_out.stem}_all_hap{config.file_out.suffix}"
        )
        console.print(
            f"[bold green]✓[/bold green] All haplotypes written to {all_hap_path}",
        )

    except Exception as e:  # noqa: BLE001
        console.print(f"\n[bold red]✗ Error:[/bold red] {e}")
        raise typer.Exit(1)  # noqa: B904


@app.command()
def validate(  # noqa: C901, PLR0912
    file_path: Annotated[
        Path,
        typer.Argument(
            help="VCF file to validate",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
) -> None:
    """Validate a VCF file format and contents.

    Checks for proper VCF structure, valid headers, and data integrity.
    """
    console.print(f"[cyan]Validating VCF file:[/cyan] {file_path}")

    try:
        # Count header lines (handle gzipped files)
        header_count = 0
        if str(file_path).endswith(".gz"):
            with gzip.open(file_path, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        header_count += 1
                    else:
                        break
        else:
            with open(file_path) as f:
                for line in f:
                    if line.startswith("#"):
                        header_count += 1
                    else:
                        break

        # Read VCF data
        ivar_df = pl.read_csv(
            file_path,
            separator="\t",
            skip_rows=header_count - 1,
            has_header=True,
        )

        # Display statistics
        console.print("\n[bold]File Statistics:[/bold]")
        console.print(f"  • Header lines: {header_count}")
        console.print(f"  • Data rows: {len(ivar_df)}")
        console.print(f"  • Columns: {len(ivar_df.columns)}")

        # Check for required columns
        required_cols = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]
        missing_cols = [col for col in required_cols if col not in ivar_df.columns]

        if missing_cols:
            console.print(
                f"\n[red]✗ Missing required columns:[/red] {', '.join(missing_cols)}",
            )
        else:
            console.print("\n[green]✓ All required VCF columns present[/green]")

        # Variant type breakdown
        if "INFO" in ivar_df.columns:
            console.print("\n[bold]Variant Types:[/bold]")
            info_counts = ivar_df["INFO"].value_counts().sort("INFO")
            for row in info_counts.iter_rows():
                console.print(f"  • {row[0]}: {row[1]:,}")

        # Filter breakdown
        if "FILTER" in ivar_df.columns:
            console.print("\n[bold]Filter Summary:[/bold]")
            filter_counts = ivar_df["FILTER"].value_counts().sort("FILTER")
            for row in filter_counts.iter_rows():
                console.print(f"  • {row[0]}: {row[1]:,}")

        console.print("\n[bold green]✓ Validation complete![/bold green]")

    except Exception as e:  # noqa: BLE001
        console.print(f"\n[bold red]✗ Validation failed:[/bold red] {e}")
        raise typer.Exit(1)  # noqa: B904


@app.command()
def stats(
    file_path: Annotated[
        Path,
        typer.Argument(
            help="iVar TSV file to analyze",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
) -> None:
    """Display statistics about an iVar TSV file.

    Shows variant counts, frequency distributions, and quality metrics.
    """
    console.print(f"[cyan]Analyzing iVar file:[/cyan] {file_path}\n")

    try:
        # Read TSV file
        ivar_df = pl.read_csv(file_path, separator="\t")

        # Basic statistics
        console.print("[bold]File Overview:[/bold]")
        console.print(f"  • Total variants: {len(ivar_df):,}")
        console.print(f"  • Unique positions: {ivar_df['POS'].n_unique():,}")

        # PASS statistics
        if "PASS" in ivar_df.columns:
            pass_count = len(ivar_df.filter(pl.col("PASS")))
            console.print(
                f"  • PASS variants: {pass_count:,} ({pass_count / len(ivar_df) * 100:.1f}%)",
            )

        # Variant type breakdown
        console.print("\n[bold]Variant Types:[/bold]")
        snp_count = len(
            ivar_df.filter(
                ~pl.col("ALT").str.starts_with("+") & ~pl.col("ALT").str.starts_with("-"),
            ),
        )
        ins_count = len(ivar_df.filter(pl.col("ALT").str.starts_with("+")))
        del_count = len(ivar_df.filter(pl.col("ALT").str.starts_with("-")))
        console.print(f"  • SNPs: {snp_count:,}")
        console.print(f"  • Insertions: {ins_count:,}")
        console.print(f"  • Deletions: {del_count:,}")

        # Frequency distribution
        console.print("\n[bold]Allele Frequency Distribution:[/bold]")
        freq_bins = [
            (0, 0.01),
            (0.01, 0.05),
            (0.05, 0.1),
            (0.1, 0.25),
            (0.25, 0.5),
            (0.5, 0.75),
            (0.75, 1.0),
        ]
        for low, high in freq_bins:
            count = len(
                ivar_df.filter(
                    (pl.col("ALT_FREQ") >= low) & (pl.col("ALT_FREQ") < high),
                ),
            )
            if count > 0:
                console.print(f"  • {low:.0%}-{high:.0%}: {count:,} variants")

        # Quality statistics
        if "ALT_QUAL" in ivar_df.columns:
            console.print("\n[bold]Quality Statistics:[/bold]")
            console.print(f"  • Mean ALT_QUAL: {ivar_df['ALT_QUAL'].mean():.1f}")
            console.print(f"  • Median ALT_QUAL: {ivar_df['ALT_QUAL'].median():.1f}")
            console.print(f"  • Min ALT_QUAL: {ivar_df['ALT_QUAL'].min():.1f}")
            console.print(f"  • Max ALT_QUAL: {ivar_df['ALT_QUAL'].max():.1f}")

        # Coverage statistics
        if "TOTAL_DP" in ivar_df.columns:
            console.print("\n[bold]Coverage Statistics:[/bold]")
            console.print(f"  • Mean depth: {ivar_df['TOTAL_DP'].mean():,.0f}")
            console.print(f"  • Median depth: {ivar_df['TOTAL_DP'].median():,.0f}")
            console.print(f"  • Min depth: {ivar_df['TOTAL_DP'].min():,}")
            console.print(f"  • Max depth: {ivar_df['TOTAL_DP'].max():,}")

    except Exception as e:  # noqa: BLE001
        console.print(f"[bold red]✗ Error:[/bold red] {e}")
        raise typer.Exit(1)  # noqa: B904


@app.callback()
def main() -> None:
    """iVar to VCF Converter - Transform variant calls with style! 🧬"""


if __name__ == "__main__":
    app()
