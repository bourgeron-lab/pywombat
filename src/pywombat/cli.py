"""CLI for wombat tool."""

import gzip
import re
import warnings
from pathlib import Path
from typing import Optional

import click
import polars as pl
import yaml


@click.group()
def cli():
    """
    Wombat: A tool for processing bcftools tabulated TSV files.

    \b
    Commands:
        filter   Process and filter variant data
        prepare  Convert TSV to optimized Parquet format
    """
    pass


@cli.command("prepare")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output Parquet file path.",
)
@click.option(
    "--chunk-size",
    type=int,
    default=50000,
    help="Number of rows to process at a time (default: 50000).",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def prepare_cmd(
    input_file: Path,
    output: Path,
    chunk_size: int,
    verbose: bool,
):
    """
    Convert bcftools TSV to optimized Parquet format.

    This command pre-processes a TSV file by:

    \b
    1. Extracting all INFO fields from the '(null)' column into separate columns
    2. Applying memory-efficient data types (Categorical for CHROM, UInt32 for POS)
    3. Writing to Parquet format for efficient columnar access

    The output Parquet file can then be used with 'wombat filter' for much faster
    and more memory-efficient filtering, especially for large files.

    \b
    Examples:
        wombat prepare input.tsv.gz -o prepared.parquet
        wombat prepare input.tsv.gz -o prepared.parquet --chunk-size 100000
    """
    try:
        if verbose:
            click.echo(f"Preparing {input_file} -> {output}", err=True)

        # Ensure output has .parquet extension
        if not str(output).endswith(".parquet"):
            output = Path(f"{output}.parquet")

        # Process the file
        prepare_parquet(input_file, output, chunk_size, verbose)

        if verbose:
            click.echo(f"Successfully created {output}", err=True)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


def prepare_parquet(
    input_file: Path,
    output: Path,
    chunk_size: int = 50000,
    verbose: bool = False,
) -> None:
    """
    Convert a bcftools TSV file to Parquet with pre-expanded INFO fields.

    Processes the file in chunks to handle large files without running out of memory.

    Args:
        input_file: Path to input TSV or TSV.gz file
        output: Path to output Parquet file
        chunk_size: Number of rows to process per chunk
        verbose: Whether to print progress
    """
    from tqdm import tqdm

    # First pass: discover all INFO fields
    if verbose:
        click.echo("Pass 1: Discovering INFO fields...", err=True)

    all_fields = set()
    all_flags = set()
    total_lines = 0

    is_gzipped = str(input_file).endswith(".gz")
    opener = gzip.open if is_gzipped else open

    with opener(input_file, "rt") as f:
        header_line = f.readline().strip()
        header_cols = header_line.split("\t")

        # Find the (null) column index dynamically
        null_col_idx = None
        for i, col in enumerate(header_cols):
            if col == "(null)":
                null_col_idx = i
                break

        if null_col_idx is None:
            if verbose:
                click.echo("Warning: No (null) column found in input", err=True)
        else:
            for line in tqdm(f, desc="Scanning", disable=not verbose):
                total_lines += 1
                parts = line.split("\t")
                if len(parts) > null_col_idx:
                    null_value = parts[null_col_idx]
                    if null_value and null_value != ".":
                        pairs = null_value.split(";")
                        for pair in pairs:
                            if "=" in pair:
                                field_name = pair.split("=", 1)[0]
                                all_fields.add(field_name)
                            elif pair.strip():
                                all_flags.add(pair.strip())

    if verbose:
        click.echo(
            f"Found {len(all_fields)} key-value fields and {len(all_flags)} flags in {total_lines} variants",
            err=True,
        )

    # Second pass: process chunks and write Parquet
    if verbose:
        click.echo("Pass 2: Converting to Parquet...", err=True)

    # Define memory-efficient dtypes
    dtype_overrides = {
        "#CHROM": pl.Categorical,
        "POS": pl.UInt32,
        "FILTER": pl.Categorical,
    }

    # Create a temporary directory for chunk files
    import tempfile
    import shutil

    temp_dir = Path(tempfile.mkdtemp(prefix="wombat_prepare_"))
    part_files = []

    try:
        with opener(input_file, "rt") as f:
            header_line = f.readline().strip()

            # Process in chunks
            chunk_lines = []
            pbar = tqdm(total=total_lines, desc="Converting", disable=not verbose)

            for line in f:
                chunk_lines.append(line)
                if len(chunk_lines) >= chunk_size:
                    df_chunk = _process_chunk(
                        header_line, chunk_lines, all_fields, all_flags, dtype_overrides
                    )
                    part_file = temp_dir / f"part_{len(part_files):06d}.parquet"
                    df_chunk.write_parquet(part_file)
                    part_files.append(part_file)
                    pbar.update(len(chunk_lines))
                    chunk_lines = []

            # Process remaining lines
            if chunk_lines:
                df_chunk = _process_chunk(
                    header_line, chunk_lines, all_fields, all_flags, dtype_overrides
                )
                part_file = temp_dir / f"part_{len(part_files):06d}.parquet"
                df_chunk.write_parquet(part_file)
                part_files.append(part_file)
                pbar.update(len(chunk_lines))

            pbar.close()

        # Combine all parts into final output using lazy scanning
        if verbose:
            click.echo(f"Combining {len(part_files)} parts into final output...", err=True)

        if part_files:
            # Use scan_parquet to lazily read all parts and write combined output
            combined = pl.scan_parquet(part_files).collect()
            combined.write_parquet(output)

        if verbose:
            click.echo(f"Wrote {len(part_files)} chunks to {output}", err=True)

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


def _process_chunk(
    header: str,
    lines: list,
    fields: set,
    flags: set,
    dtype_overrides: dict,
) -> pl.DataFrame:
    """Process a chunk of lines into a DataFrame with expanded INFO fields."""
    import io

    content = header + "\n" + "".join(lines)
    df = pl.read_csv(
        io.StringIO(content),
        separator="\t",
        infer_schema_length=10000,
    )

    # Expand INFO fields from (null) column
    if "(null)" in df.columns:
        # Extract key-value fields
        for field in sorted(fields):
            df = df.with_columns(
                pl.col("(null)").str.extract(f"{field}=([^;]+)").alias(field)
            )

        # Extract boolean flags
        for flag in sorted(flags):
            df = df.with_columns(
                pl.col("(null)").str.contains(f"(^|;){flag}(;|$)").alias(flag)
            )

        # Drop the original (null) column
        df = df.drop("(null)")

    # Drop CSQ column if it exists (redundant after expansion)
    if "CSQ" in df.columns:
        df = df.drop("CSQ")

    # Apply memory-efficient dtypes
    for col, dtype in dtype_overrides.items():
        if col in df.columns:
            try:
                df = df.with_columns(pl.col(col).cast(dtype))
            except Exception:
                pass  # Skip if cast fails

    return df


@cli.command("filter")
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=str,
    help="Output file prefix. If not specified, generates from input filename.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["tsv", "tsv.gz", "parquet"], case_sensitive=False),
    default="tsv",
    help="Output format: tsv (default), tsv.gz (compressed), or parquet.",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
@click.option(
    "-p",
    "--pedigree",
    type=click.Path(exists=True, path_type=Path),
    help="Pedigree file to add father and mother genotype columns.",
)
@click.option(
    "-F",
    "--filter-config",
    type=click.Path(exists=True, path_type=Path),
    help="Filter configuration YAML file to apply quality and impact filters.",
)
@click.option(
    "--debug",
    type=str,
    help="Debug mode: show rows matching chrom:pos (e.g., chr11:70486013).",
)
def filter_cmd(
    input_file: Path,
    output: Optional[str],
    output_format: str,
    verbose: bool,
    pedigree: Optional[Path],
    filter_config: Optional[Path],
    debug: Optional[str],
):
    """
    Process and filter variant data from TSV or Parquet files.

    \b
    Supports two input formats:
    - TSV/TSV.gz: Full processing (INFO expansion + melting)
    - Parquet: Fast processing (melting only, INFO already expanded)

    \b
    For large files, use 'wombat prepare' first to convert to Parquet,
    then use 'wombat filter' on the Parquet file for better performance.

    \b
    This command:
    1. Expands the '(null)' column (TSV only) into separate columns
    2. Melts sample columns into rows with sample names
    3. Splits sample values (GT:DP:GQ:AD format) into separate columns
    4. Applies quality and expression filters (if config provided)

    \b
    Examples:
        wombat filter input.tsv -o output
        wombat filter prepared.parquet -o output -f parquet
        wombat filter input.tsv -p pedigree.tsv -F config.yml
    """
    try:
        if verbose:
            click.echo(f"Reading input file: {input_file}", err=True)

        # Detect input format
        is_parquet = str(input_file).endswith(".parquet")
        is_gzipped = str(input_file).endswith(".gz")

        if verbose:
            if is_parquet:
                click.echo("Detected Parquet input (pre-processed)", err=True)
            elif is_gzipped:
                click.echo("Detected gzipped TSV file", err=True)

        # Read pedigree file if provided
        pedigree_df = None
        if pedigree:
            if verbose:
                click.echo(f"Reading pedigree file: {pedigree}", err=True)
            pedigree_df = read_pedigree(pedigree)

        # Load filter config if provided
        filter_config_data = None
        if filter_config:
            if verbose:
                click.echo(f"Reading filter config: {filter_config}", err=True)
            filter_config_data = load_filter_config(filter_config)

        # Debug mode: show specific variant
        if debug:
            debug_variant(input_file, pedigree_df, filter_config_data, debug, verbose)
            return

        # Determine output prefix
        if output is None:
            # Generate default output prefix from input filename
            input_stem = input_file.name
            # Remove known extensions
            for ext in [".tsv.gz", ".tsv", ".parquet"]:
                if input_stem.endswith(ext):
                    input_stem = input_stem[: -len(ext)]
                    break

            # Add config name if filter is provided
            if filter_config:
                config_name = filter_config.stem  # Get basename without extension
                output = f"{input_stem}.{config_name}"
            else:
                output = input_stem

        # Process using streaming mode
        if verbose:
            click.echo("Processing with streaming mode...", err=True)

        # Build lazy query based on input format
        if is_parquet:
            # Parquet input: INFO fields already expanded by 'wombat prepare'
            lazy_df = pl.scan_parquet(input_file)

            # OPTIMIZATION: Apply expression filter BEFORE melting
            # Expression filters (VEP_IMPACT, etc.) don't depend on sample data
            if filter_config_data and "expression" in filter_config_data:
                expression = filter_config_data["expression"]
                if expression and verbose:
                    click.echo(
                        f"Applying expression filter before melting: {expression}",
                        err=True,
                    )

                # Collect a small sample to get schema for expression parsing
                schema_df = lazy_df.head(1).collect()
                try:
                    filter_expr = parse_impact_filter_expression(expression, schema_df)
                    lazy_df = lazy_df.filter(filter_expr)

                    # Count filtered variants
                    if verbose:
                        filtered_count = lazy_df.select(pl.len()).collect().item()
                        click.echo(
                            f"Variants after expression filter: {filtered_count}",
                            err=True,
                        )
                except ValueError as e:
                    if verbose:
                        click.echo(
                            f"Warning: Could not apply early filter: {e}", err=True
                        )

            # Now collect and melt (on filtered variants only)
            df = lazy_df.collect()
            formatted_df = format_bcftools_tsv_minimal(df, pedigree_df)
            lazy_df = formatted_df.lazy()

            # Remove expression from config so it's not applied again
            if filter_config_data and "expression" in filter_config_data:
                filter_config_data = filter_config_data.copy()
                del filter_config_data["expression"]
        else:
            # TSV input: need full processing (melt + annotation expansion)
            string_columns = [
                "FID",
                "sample_id",
                "father_id",
                "mother_id",
                "FatherBarcode",
                "MotherBarcode",
                "sample",
            ]
            schema_overrides = {col: pl.Utf8 for col in string_columns}
            lazy_df = pl.scan_csv(
                input_file, separator="\t", schema_overrides=schema_overrides
            )

            # Apply formatting transformations (melt + expand annotations)
            lazy_df = format_bcftools_tsv_lazy(lazy_df, pedigree_df)

        # Apply filters if provided
        if filter_config_data:
            lazy_df = apply_filters_lazy(
                lazy_df, filter_config_data, verbose, pedigree_df
            )

        # Write output
        output_path = Path(f"{output}.{output_format}")

        if output_format == "tsv":
            lazy_df.sink_csv(output_path, separator="\t")
        elif output_format == "tsv.gz":
            # For gzip, we need to collect and write
            df = lazy_df.collect()
            csv_content = df.write_csv(separator="\t")
            with gzip.open(output_path, "wt") as f:
                f.write(csv_content)
        elif output_format == "parquet":
            lazy_df.sink_parquet(output_path)

        if verbose:
            click.echo(f"Data written to {output_path}", err=True)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


def debug_variant(
    input_file: Path,
    pedigree_df: Optional[pl.DataFrame],
    filter_config: Optional[dict],
    debug_pos: str,
    verbose: bool,
):
    """Debug mode: display rows matching a specific chrom:pos."""
    # Parse debug position
    if ":" not in debug_pos:
        click.echo(
            "Error: Debug position must be in format 'chrom:pos' (e.g., chr11:70486013)",
            err=True,
        )
        raise click.Abort()

    chrom, pos = debug_pos.split(":", 1)
    try:
        pos = int(pos)
    except ValueError:
        click.echo(f"Error: Position must be an integer, got '{pos}'", err=True)
        raise click.Abort()

    if verbose:
        click.echo(f"Debug mode: searching for {chrom}:{pos}", err=True)

    # Read and format the data
    # Force certain columns to string type
    string_columns = [
        "FID",
        "sample_id",
        "father_id",
        "mother_id",
        "FatherBarcode",
        "MotherBarcode",
        "sample",
    ]
    schema_overrides = {col: pl.Utf8 for col in string_columns}
    df = pl.read_csv(input_file, separator="\t", schema_overrides=schema_overrides)
    formatted_df = format_bcftools_tsv(df, pedigree_df)

    # Filter to matching rows
    matching_rows = formatted_df.filter(
        (pl.col("#CHROM") == chrom) & (pl.col("POS") == pos)
    )

    if matching_rows.shape[0] == 0:
        click.echo(f"No rows found matching {chrom}:{pos}", err=True)
        return

    # Determine which columns to display
    columns_to_show = ["#CHROM", "POS"]

    # Add VEP_SYMBOL if it exists
    if "VEP_SYMBOL" in matching_rows.columns:
        columns_to_show.append("VEP_SYMBOL")

    # Extract column names from expression if filter config provided
    if filter_config and "expression" in filter_config:
        expression = filter_config["expression"]
        # Extract column names from expression using regex
        # Match patterns like "column_name" before operators
        column_pattern = r"\b([A-Za-z_][A-Za-z0-9_]*)\b\s*[=!<>]"
        found_columns = re.findall(column_pattern, expression)

        for col in found_columns:
            if col in matching_rows.columns and col not in columns_to_show:
                columns_to_show.append(col)

    # Select only the columns we want to display
    display_df = matching_rows.select(
        [col for col in columns_to_show if col in matching_rows.columns]
    )

    # Replace null and NaN values with <null> and <NaN> for display
    for col in display_df.columns:
        if display_df[col].dtype in [pl.Float32, pl.Float64]:
            # For numeric columns, handle both NaN and null
            display_df = display_df.with_columns(
                pl.when(pl.col(col).is_null())
                .then(pl.lit("<null>"))
                .when(pl.col(col).is_nan())
                .then(pl.lit("<NaN>"))
                .otherwise(pl.col(col).cast(pl.Utf8))
                .alias(col)
            )
        else:
            # For non-numeric columns, only handle null
            display_df = display_df.with_columns(
                pl.when(pl.col(col).is_null())
                .then(pl.lit("<null>"))
                .otherwise(pl.col(col).cast(pl.Utf8))
                .alias(col)
            )

    # Display the results
    click.echo(f"\nFound {matching_rows.shape[0]} row(s) matching {chrom}:{pos}:\n")
    click.echo(display_df.write_csv(separator="\t"))


def load_filter_config(config_path: Path) -> dict:
    """Load and parse filter configuration from YAML file."""
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def apply_quality_filters(
    df: pl.DataFrame, quality_config: dict, verbose: bool = False
) -> pl.DataFrame:
    """Apply quality filters to the dataframe."""
    if quality_config is None:
        return df

    original_rows = df.shape[0]

    # Check for rare genotypes with '2' and warn
    if "2" in str(df["sample_gt"].to_list()):
        rare_gts = df.filter(pl.col("sample_gt").str.contains("2"))
        if rare_gts.shape[0] > 0:
            warnings.warn(
                f"Found {rare_gts.shape[0]} rows with rare genotypes containing '2'. These will be kept."
            )

    # Filter: sample_gt must contain at least one '1' (default: true)
    filter_no_alt = quality_config.get("filter_no_alt_allele", True)
    if filter_no_alt:
        df = df.filter(
            pl.col("sample_gt").str.contains("1")
            | pl.col("sample_gt").str.contains("2")
        )

    # Apply minimum depth filter
    if "sample_dp_min" in quality_config:
        min_dp = quality_config["sample_dp_min"]
        df = df.filter(pl.col("sample_dp").cast(pl.Float64, strict=False) >= min_dp)

    # Apply minimum GQ filter
    if "sample_gq_min" in quality_config:
        min_gq = quality_config["sample_gq_min"]
        df = df.filter(pl.col("sample_gq").cast(pl.Float64, strict=False) >= min_gq)

    # Determine genotype for VAF filters
    # Het: contains exactly one '1' (0/1 or 1/0)
    # HomAlt: 1/1
    # HomRef: 0/0
    is_het = (pl.col("sample_gt").str.count_matches("1") == 1) & ~pl.col(
        "sample_gt"
    ).str.contains("2")
    is_hom_alt = pl.col("sample_gt") == "1/1"
    is_hom_ref = pl.col("sample_gt") == "0/0"

    # VAF filters for heterozygous
    if "sample_vaf_het_min" in quality_config:
        min_vaf_het = quality_config["sample_vaf_het_min"]
        df = df.filter(~is_het | (pl.col("sample_vaf") >= min_vaf_het))

    if "sample_vaf_het_max" in quality_config:
        max_vaf_het = quality_config["sample_vaf_het_max"]
        df = df.filter(~is_het | (pl.col("sample_vaf") <= max_vaf_het))

    # VAF filters for homozygous alternate
    if "sample_vaf_homalt_min" in quality_config:
        min_vaf_homalt = quality_config["sample_vaf_homalt_min"]
        df = df.filter(~is_hom_alt | (pl.col("sample_vaf") >= min_vaf_homalt))

    # VAF filters for homozygous reference (wild type)
    if "sample_vaf_hom_ref_max" in quality_config:
        max_vaf_hom_ref = quality_config["sample_vaf_hom_ref_max"]
        df = df.filter(~is_hom_ref | (pl.col("sample_vaf") <= max_vaf_hom_ref))

    # Apply filters to parents if they exist and option is enabled
    apply_to_parents = quality_config.get("apply_to_parents", False)
    if apply_to_parents and "father" in df.columns:
        # Apply same filters to father columns
        if "sample_dp_min" in quality_config:
            min_dp = quality_config["sample_dp_min"]
            df = df.filter(
                (pl.col("father_dp").is_null())
                | (pl.col("father_dp").cast(pl.Float64, strict=False) >= min_dp)
            )

        if "sample_gq_min" in quality_config:
            min_gq = quality_config["sample_gq_min"]
            df = df.filter(
                (pl.col("father_gq").is_null())
                | (pl.col("father_gq").cast(pl.Float64, strict=False) >= min_gq)
            )

        # Father genotype checks
        father_is_het = (pl.col("father_gt").str.count_matches("1") == 1) & ~pl.col(
            "father_gt"
        ).str.contains("2")
        father_is_hom_alt = pl.col("father_gt") == "1/1"
        father_is_hom_ref = pl.col("father_gt") == "0/0"

        if "sample_vaf_het_min" in quality_config:
            min_vaf_het = quality_config["sample_vaf_het_min"]
            df = df.filter(
                pl.col("father_vaf").is_null()
                | ~father_is_het
                | (pl.col("father_vaf") >= min_vaf_het)
            )

        if "sample_vaf_het_max" in quality_config:
            max_vaf_het = quality_config["sample_vaf_het_max"]
            df = df.filter(
                pl.col("father_vaf").is_null()
                | ~father_is_het
                | (pl.col("father_vaf") <= max_vaf_het)
            )

        if "sample_vaf_homalt_min" in quality_config:
            min_vaf_homalt = quality_config["sample_vaf_homalt_min"]
            df = df.filter(
                pl.col("father_vaf").is_null()
                | ~father_is_hom_alt
                | (pl.col("father_vaf") >= min_vaf_homalt)
            )

        if "sample_vaf_hom_ref_max" in quality_config:
            max_vaf_hom_ref = quality_config["sample_vaf_hom_ref_max"]
            df = df.filter(
                pl.col("father_vaf").is_null()
                | ~father_is_hom_ref
                | (pl.col("father_vaf") <= max_vaf_hom_ref)
            )

        # Apply same filters to mother columns
        if "sample_dp_min" in quality_config:
            min_dp = quality_config["sample_dp_min"]
            df = df.filter(
                (pl.col("mother_dp").is_null())
                | (pl.col("mother_dp").cast(pl.Float64, strict=False) >= min_dp)
            )

        if "sample_gq_min" in quality_config:
            min_gq = quality_config["sample_gq_min"]
            df = df.filter(
                (pl.col("mother_gq").is_null())
                | (pl.col("mother_gq").cast(pl.Float64, strict=False) >= min_gq)
            )

        # Mother genotype checks
        mother_is_het = (pl.col("mother_gt").str.count_matches("1") == 1) & ~pl.col(
            "mother_gt"
        ).str.contains("2")
        mother_is_hom_alt = pl.col("mother_gt") == "1/1"
        mother_is_hom_ref = pl.col("mother_gt") == "0/0"

        if "sample_vaf_het_min" in quality_config:
            min_vaf_het = quality_config["sample_vaf_het_min"]
            df = df.filter(
                pl.col("mother_vaf").is_null()
                | ~mother_is_het
                | (pl.col("mother_vaf") >= min_vaf_het)
            )

        if "sample_vaf_het_max" in quality_config:
            max_vaf_het = quality_config["sample_vaf_het_max"]
            df = df.filter(
                pl.col("mother_vaf").is_null()
                | ~mother_is_het
                | (pl.col("mother_vaf") <= max_vaf_het)
            )

        if "sample_vaf_homalt_min" in quality_config:
            min_vaf_homalt = quality_config["sample_vaf_homalt_min"]
            df = df.filter(
                pl.col("mother_vaf").is_null()
                | ~mother_is_hom_alt
                | (pl.col("mother_vaf") >= min_vaf_homalt)
            )

        if "sample_vaf_hom_ref_max" in quality_config:
            max_vaf_hom_ref = quality_config["sample_vaf_hom_ref_max"]
            df = df.filter(
                pl.col("mother_vaf").is_null()
                | ~mother_is_hom_ref
                | (pl.col("mother_vaf") <= max_vaf_hom_ref)
            )

    if verbose:
        filtered_rows = df.shape[0]
        click.echo(
            f"Quality filters: {original_rows} -> {filtered_rows} rows ({original_rows - filtered_rows} filtered out)",
            err=True,
        )

    return df


# ------------------ De novo (DNM) filter helpers ------------------


def _chrom_short(chrom: str) -> str:
    """Normalize chromosome name to short form (e.g., 'chrX' -> 'X')."""
    if chrom is None:
        return ""
    chrom = str(chrom)
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _pos_in_par(chrom: str, pos: int, par_regions: dict) -> bool:
    """Return True if (chrom,pos) falls in any PAR region from config.

    Normalizes chromosome names to match both 'X'/'chrX' formats.
    """
    if not par_regions:
        return False

    chrom_short = _chrom_short(chrom)

    for assembly, regions in par_regions.items():
        for region_name, region in regions.items():
            region_chrom = _chrom_short(region.get("chrom", "X"))
            start = int(region.get("start"))
            end = int(region.get("end"))
            # Normalize both to uppercase for comparison
            if region_chrom.upper() == chrom_short.upper() and start <= pos <= end:
                return True
    return False


def apply_de_novo_filter(
    df: pl.DataFrame,
    dnm_config: dict,
    verbose: bool = False,
    pedigree_df: Optional[pl.DataFrame] = None,
) -> pl.DataFrame:
    """Apply de novo detection filters to dataframe using vectorized operations.

    dnm_config expected keys:
      - sample_dp_min, sample_gq_min, sample_vaf_min
      - parent_dp_min, parent_gq_min, parent_vaf_max
      - par_regions: dict with PAR regions keyed by assembly

    This function will read `sex` from `df` when present; otherwise it will use
    the `pedigree_df` (which should contain `sample_id` and `sex`).
    """
    if not dnm_config:
        return df

    # Required thresholds
    s_dp = dnm_config.get("sample_dp_min", 10)
    s_gq = dnm_config.get("sample_gq_min", 18)
    s_vaf = dnm_config.get("sample_vaf_min", 0.15)
    s_vaf_hemizygous = dnm_config.get("sample_vaf_hemizygous_min", 0.85)

    p_dp = dnm_config.get("parent_dp_min", 10)
    p_gq = dnm_config.get("parent_gq_min", 18)
    p_vaf = dnm_config.get("parent_vaf_max", 0.02)

    fafmax_max = dnm_config.get("fafmax_faf95_max_genomes_max", None)
    genomes_filters_pass_only = dnm_config.get("genomes_filters_pass_only", False)

    par_regions = dnm_config.get("par_regions", {})

    original = df.shape[0]

    # Ensure we have parent identifiers (father/mother). Try to add from pedigree if missing.
    if "father_id" not in df.columns or "mother_id" not in df.columns:
        if pedigree_df is not None:
            # Join pedigree to get father/mother/sample sex if needed
            df = df.join(
                pedigree_df, left_on="sample", right_on="sample_id", how="left"
            )
        else:
            raise click.Abort(
                "DNM filtering requires a pedigree (with father/mother IDs and sex)."
            )

    if verbose:
        click.echo(
            f"DNM: Starting with {df.shape[0]} variants after pedigree join", err=True
        )

    # Optimization #6: Early exit for variants with missing parents
    df = df.filter(
        pl.col("father_id").is_not_null() & pl.col("mother_id").is_not_null()
    )

    if verbose:
        click.echo(
            f"DNM: {df.shape[0]} variants after filtering missing parents", err=True
        )

    if df.shape[0] == 0:
        if verbose:
            click.echo(
                f"De novo filter: {original} -> 0 rows (all missing parents)", err=True
            )
        return df

    # Ensure we have sex information
    if "sex" not in df.columns:
        if pedigree_df is not None and "sex" in pedigree_df.columns:
            # Re-join to get sex if not already present
            if "sex" not in df.columns:
                df = df.join(
                    pedigree_df.select(["sample_id", "sex"]),
                    left_on="sample",
                    right_on="sample_id",
                    how="left",
                )
        else:
            raise click.Abort(
                "DNM filtering requires sex information in the pedigree (column 'sex')."
            )

    # Filter out variants with missing/partial genotypes (containing '.')
    # Sample genotype must be fully called (no '.', './0', '1/.', './.' etc.)
    df = df.filter(
        ~pl.col("sample_gt").str.contains(r"\.") & pl.col("sample_gt").is_not_null()
    )

    if verbose:
        click.echo(
            f"DNM: {df.shape[0]} variants after removing samples with missing/partial GT",
            err=True,
        )

    # Ensure proband passes basic sample thresholds (DP, GQ)
    # VAF threshold will be applied differently for hemizygous vs diploid variants
    df = df.filter(
        (pl.col("sample_dp").cast(pl.Float64, strict=False) >= s_dp)
        & (pl.col("sample_gq").cast(pl.Float64, strict=False) >= s_gq)
    )

    if verbose:
        click.echo(
            f"DNM: {df.shape[0]} variants after proband QC (DP>={s_dp}, GQ>={s_gq})",
            err=True,
        )

    if df.shape[0] == 0:
        if verbose:
            click.echo(
                f"De novo filter: {original} -> 0 rows (sample QC failed)", err=True
            )
        return df

    # Optimization #7: Vectorized parent filtering using when/then expressions
    # Build chromosome-specific parent filters

    # Ensure #CHROM is string type for operations (convert from categorical if needed)
    if "#CHROM" in df.columns and df.schema["#CHROM"] == pl.Categorical:
        df = df.with_columns(pl.col("#CHROM").cast(pl.Utf8))

    # Normalize chromosome to short form for comparison
    df = df.with_columns(
        pl.col("#CHROM")
        .str.replace("^chr", "")
        .str.to_uppercase()
        .alias("_chrom_short")
    )

    # Determine if variant is in PAR region (vectorized)
    par_mask = pl.lit(False)
    if par_regions:
        for assembly, regions in par_regions.items():
            for region_name, region in regions.items():
                region_chrom = _chrom_short(region.get("chrom", "X")).upper()
                start = int(region.get("start"))
                end = int(region.get("end"))
                par_mask = par_mask | (
                    (pl.col("_chrom_short") == region_chrom)
                    & (pl.col("POS") >= start)
                    & (pl.col("POS") <= end)
                )

    df = df.with_columns(par_mask.alias("_in_par"))

    # Normalize sex to uppercase string for comparison
    df = df.with_columns(
        pl.col("sex").cast(pl.Utf8).str.to_uppercase().alias("_sex_norm")
    )

    # Determine if variant is hemizygous (X in males outside PAR, or Y in males)
    is_hemizygous = (
        (pl.col("_chrom_short") == "X")
        & (pl.col("_sex_norm").is_in(["1", "M"]))
        & ~pl.col("_in_par")
    ) | (pl.col("_chrom_short") == "Y")

    # Determine if variant is homozygous (1/1, 2/2, etc.)
    is_homozygous = pl.col("sample_gt").is_in(["1/1", "2/2", "3/3"])

    # Apply VAF threshold: hemizygous and homozygous variants require higher VAF (>=0.85)
    sample_vaf_filter = (
        pl.when(is_hemizygous | is_homozygous)
        .then(pl.col("sample_vaf") >= s_vaf_hemizygous)
        .otherwise(pl.col("sample_vaf") > s_vaf)
    )

    df = df.filter(sample_vaf_filter)

    if verbose:
        click.echo(
            f"DNM: {df.shape[0]} variants after VAF filtering (het>{s_vaf}, hom/hemizygous>={s_vaf_hemizygous})",
            err=True,
        )

    # Apply fafmax_faf95_max_genomes filter if specified
    if fafmax_max is not None:
        if "fafmax_faf95_max_genomes" in df.columns:
            df = df.filter(
                (
                    pl.col("fafmax_faf95_max_genomes").cast(pl.Float64, strict=False)
                    <= fafmax_max
                )
                | pl.col("fafmax_faf95_max_genomes").is_null()
            )
            if verbose:
                click.echo(
                    f"DNM: {df.shape[0]} variants after fafmax_faf95_max_genomes filter (<={fafmax_max})",
                    err=True,
                )
        elif verbose:
            click.echo(
                "DNM: Warning - fafmax_faf95_max_genomes column not found, skipping frequency filter",
                err=True,
            )

    # Apply genomes_filters filter if specified
    if genomes_filters_pass_only:
        if "genomes_filters" in df.columns:
            df = df.filter(
                (pl.col("genomes_filters") == ".") | pl.col("genomes_filters").is_null()
            )
            if verbose:
                click.echo(
                    f"DNM: {df.shape[0]} variants after genomes_filters filter (pass only)",
                    err=True,
                )
        elif verbose:
            click.echo(
                "DNM: Warning - genomes_filters column not found, skipping genomes_filters filter",
                err=True,
            )

    # Build parent quality checks (common to all)
    father_qual_ok = (pl.col("father_dp").cast(pl.Float64, strict=False) >= p_dp) & (
        pl.col("father_gq").cast(pl.Float64, strict=False) >= p_gq
    )
    mother_qual_ok = (pl.col("mother_dp").cast(pl.Float64, strict=False) >= p_dp) & (
        pl.col("mother_gq").cast(pl.Float64, strict=False) >= p_gq
    )

    father_vaf_ok = pl.col("father_vaf").is_null() | (pl.col("father_vaf") < p_vaf)
    mother_vaf_ok = pl.col("mother_vaf").is_null() | (pl.col("mother_vaf") < p_vaf)

    # Parent genotype checks: ensure no '.' in genotypes (no ./., 0/., 1/. etc.)
    father_gt_ok = (
        ~pl.col("father_gt").str.contains(r"\.") & pl.col("father_gt").is_not_null()
    )
    mother_gt_ok = (
        ~pl.col("mother_gt").str.contains(r"\.") & pl.col("mother_gt").is_not_null()
    )

    # Build comprehensive parent filter using when/then logic
    parent_filter = (
        pl.when(pl.col("_in_par"))
        # PAR region: both parents must be reference (autosomal-like) with valid GTs
        .then(
            father_qual_ok
            & father_vaf_ok
            & father_gt_ok
            & mother_qual_ok
            & mother_vaf_ok
            & mother_gt_ok
        )
        .when(pl.col("_chrom_short") == "Y")
        # Y chromosome: only check father (mother doesn't have Y), father GT must be valid
        .then(father_qual_ok & father_vaf_ok & father_gt_ok)
        .when((pl.col("_chrom_short") == "X") & (pl.col("_sex_norm").is_in(["1", "M"])))
        # X chromosome, male proband: father is hemizygous, only check mother VAF and GT
        .then(father_qual_ok & mother_qual_ok & mother_vaf_ok & mother_gt_ok)
        # Default (autosomes or X/female): both parents must be reference with valid GTs
        .otherwise(
            father_qual_ok
            & father_vaf_ok
            & father_gt_ok
            & mother_qual_ok
            & mother_vaf_ok
            & mother_gt_ok
        )
    )

    # Apply parent filter
    result = df.filter(parent_filter)

    if verbose:
        click.echo(
            f"DNM: {result.shape[0]} variants after parent genotype filtering (parent DP>={p_dp}, GQ>={p_gq}, VAF<{p_vaf})",
            err=True,
        )

    # Drop temporary columns
    result = result.drop(["_chrom_short", "_in_par", "_sex_norm"])

    if verbose:
        click.echo(f"De novo filter: {original} -> {result.shape[0]} rows", err=True)

    return result


def parse_impact_filter_expression(expression: str, df: pl.DataFrame) -> pl.Expr:
    """Parse a filter expression string into a Polars expression."""
    # Replace operators with Polars equivalents
    # Support: =, !=, <=, >=, <, >, &, |, ()

    expr_str = expression.strip()

    # Split by logical operators while preserving them
    tokens = re.split(r"(\s*[&|]\s*|\(|\))", expr_str)
    tokens = [t.strip() for t in tokens if t.strip()]

    def parse_condition(condition: str) -> pl.Expr:
        """Parse a single condition into a Polars expression."""
        condition = condition.strip()

        # Try different operators in order of specificity
        for op in ["<=", ">=", "!=", "=", "<", ">"]:
            if op in condition:
                parts = condition.split(op, 1)
                if len(parts) == 2:
                    col_name = parts[0].strip()
                    value = parts[1].strip()

                    # Check if column exists
                    if col_name not in df.columns:
                        raise ValueError(f"Column '{col_name}' not found in dataframe")

                    # Check for null value
                    if value.upper() == "NULL":
                        col_expr = pl.col(col_name)
                        if op == "=":
                            return col_expr.is_null()
                        elif op == "!=":
                            return ~col_expr.is_null()
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for null comparison, use = or !="
                            )

                    # Check for NaN value
                    if value.upper() == "NAN":
                        col_expr = pl.col(col_name).cast(pl.Float64, strict=False)
                        if op == "=":
                            return col_expr.is_nan()
                        elif op == "!=":
                            return ~col_expr.is_nan()
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for NaN comparison, use = or !="
                            )

                    # Try to convert value to number, otherwise treat as string
                    try:
                        value_num = float(value)
                        col_expr = pl.col(col_name).cast(pl.Float64, strict=False)

                        if op == "=":
                            return col_expr == value_num
                        elif op == "!=":
                            return col_expr != value_num
                        elif op == "<=":
                            return col_expr <= value_num
                        elif op == ">=":
                            return col_expr >= value_num
                        elif op == "<":
                            return col_expr < value_num
                        elif op == ">":
                            return col_expr > value_num
                    except ValueError:
                        # String comparison (case-insensitive)
                        value = value.strip("'\"")
                        col_expr = pl.col(col_name).str.to_lowercase()
                        value_lower = value.lower()

                        if op == "=":
                            return col_expr == value_lower
                        elif op == "!=":
                            return col_expr != value_lower
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for string comparison"
                            )
                break

        raise ValueError(f"Could not parse condition: {condition}")

    def build_expression(tokens: list, idx: int = 0) -> tuple[pl.Expr, int]:
        """Recursively build expression from tokens."""
        if idx >= len(tokens):
            return None, idx

        result = None
        i = idx

        while i < len(tokens):
            token = tokens[i]

            if token == "(":
                # Parse sub-expression
                sub_expr, new_i = build_expression(tokens, i + 1)
                if result is None:
                    result = sub_expr
                i = new_i
            elif token == ")":
                # End of sub-expression
                return result, i + 1
            elif token == "&":
                # AND operator
                next_expr, new_i = build_expression(tokens, i + 1)
                if next_expr is not None:
                    result = result & next_expr if result is not None else next_expr
                return result, new_i
            elif token == "|":
                # OR operator
                next_expr, new_i = build_expression(tokens, i + 1)
                if next_expr is not None:
                    result = result | next_expr if result is not None else next_expr
                return result, new_i
            else:
                # It's a condition
                cond_expr = parse_condition(token)
                if result is None:
                    result = cond_expr
                else:
                    # If we have a result and encounter another condition without an operator,
                    # assume AND
                    result = result & cond_expr

            i += 1

        return result, i

    expr, _ = build_expression(tokens)
    return expr


def apply_impact_filters(
    df: pl.DataFrame,
    impact_config: list,
    output_prefix: str,
    output_format: str,
    verbose: bool,
):
    """Apply impact filters and create separate output files."""
    if not impact_config:
        return

    # Sort impact filters by priority (lower number = higher priority)
    impact_filters = sorted(impact_config, key=lambda x: x.get("priority", 999))

    # Create a dict to store variants by impact filter
    impact_variants = {}

    # Apply each impact filter
    for impact_filter in impact_filters:
        name = impact_filter["name"]
        priority = impact_filter.get("priority", 999)
        expression = impact_filter["expression"]

        if verbose:
            click.echo(
                f"Applying impact filter '{name}' (priority {priority})...", err=True
            )

        try:
            # Parse and apply the filter expression
            filter_expr = parse_impact_filter_expression(expression, df)
            filtered_df = df.filter(filter_expr)

            if verbose:
                click.echo(
                    f"  Impact filter '{name}': {filtered_df.shape[0]} variants",
                    err=True,
                )

            # Store the filtered dataframe with its priority
            impact_variants[name] = {
                "df": filtered_df,
                "priority": priority,
            }
        except Exception as e:
            click.echo(
                f"Error applying impact filter '{name}': {e}",
                err=True,
            )
            raise

    # Add flag_higher_impact column to each filtered dataframe
    for name, data in impact_variants.items():
        filtered_df = data["df"]
        priority = data["priority"]

        # Find variants that appear in higher priority filters
        higher_priority_filters = []
        for other_name, other_data in impact_variants.items():
            if other_data["priority"] < priority:
                higher_priority_filters.append(other_name)

        if higher_priority_filters:
            # Create a set of variant keys from higher priority filters
            variant_keys_in_higher = set()
            for other_name in higher_priority_filters:
                other_df = impact_variants[other_name]["df"]
                for row in other_df.select(["#CHROM", "POS", "REF", "ALT"]).iter_rows():
                    variant_keys_in_higher.add(tuple(row))

            # Add flag_higher_impact column
            def check_higher_impact(chrom, pos, ref, alt):
                key = (chrom, pos, ref, alt)
                if key in variant_keys_in_higher:
                    # Find which filters it appears in
                    filters_with_variant = []
                    for other_name in higher_priority_filters:
                        other_df = impact_variants[other_name]["df"]
                        match = other_df.filter(
                            (pl.col("#CHROM") == chrom)
                            & (pl.col("POS") == pos)
                            & (pl.col("REF") == ref)
                            & (pl.col("ALT") == alt)
                        )
                        if match.shape[0] > 0:
                            filters_with_variant.append(other_name)
                    return (
                        ", ".join(filters_with_variant) if filters_with_variant else ""
                    )
                return ""

            # Add the flag column
            filtered_df = filtered_df.with_columns(
                [
                    pl.struct(["#CHROM", "POS", "REF", "ALT"])
                    .map_elements(
                        lambda x: check_higher_impact(
                            x["#CHROM"], x["POS"], x["REF"], x["ALT"]
                        ),
                        return_dtype=pl.Utf8,
                    )
                    .alias("flag_higher_impact")
                ]
            )
        else:
            # No higher priority filters
            filtered_df = filtered_df.with_columns(
                [pl.lit("").alias("flag_higher_impact")]
            )

        # Write to file
        output_filename = f"{output_prefix}.{name}.{output_format}"
        output_path = Path(output_filename)

        if output_format == "tsv":
            filtered_df.write_csv(output_path, separator="\t")
        elif output_format == "tsv.gz":
            csv_content = filtered_df.write_csv(separator="\t")
            with gzip.open(output_path, "wt") as f:
                f.write(csv_content)
        elif output_format == "parquet":
            filtered_df.write_parquet(output_path)

        click.echo(
            f"Written {filtered_df.shape[0]} variants to {output_path}", err=True
        )


def apply_filters_and_write(
    df: pl.DataFrame,
    filter_config: dict,
    output_prefix: Optional[str],
    output_format: str,
    verbose: bool,
    pedigree_df: Optional[pl.DataFrame] = None,
):
    """Apply filters and write output files."""
    # If DNM mode is enabled, apply DNM-specific criteria and skip impact/frequency filters
    if filter_config.get("dnm", {}).get("enabled", False):
        dnm_cfg = {}
        # Merge quality & dnm-specific thresholds into a single config for the function
        dnm_cfg.update(filter_config.get("quality", {}))
        dnm_cfg.update(filter_config.get("dnm", {}))

        filtered_df = apply_de_novo_filter(
            df, dnm_cfg, verbose, pedigree_df=pedigree_df
        )

        # Write result (same behavior as non-impact single-output)
        if not output_prefix:
            if output_format != "tsv":
                click.echo(
                    "Error: stdout output only supported for TSV format.",
                    err=True,
                )
                raise click.Abort()
            click.echo(filtered_df.write_csv(separator="\t"), nl=False)
        else:
            output_path = Path(f"{output_prefix}.{output_format}")

            if output_format == "tsv":
                filtered_df.write_csv(output_path, separator="\t")
            elif output_format == "tsv.gz":
                csv_content = filtered_df.write_csv(separator="\t")
                with gzip.open(output_path, "wt") as f:
                    f.write(csv_content)
            elif output_format == "parquet":
                filtered_df.write_parquet(output_path)

            click.echo(f"De novo variants written to {output_path}", err=True)

        return

    # Apply quality filters first
    quality_config = filter_config.get("quality", {})
    filtered_df = apply_quality_filters(df, quality_config, verbose)

    # Get impact filters
    impact_config = filter_config.get("impact", [])

    if not impact_config:
        # No impact filters - write single output file
        if not output_prefix:
            # Write to stdout
            if output_format != "tsv":
                click.echo(
                    "Error: stdout output only supported for TSV format.",
                    err=True,
                )
                raise click.Abort()
            click.echo(filtered_df.write_csv(separator="\t"), nl=False)
        else:
            output_path = Path(f"{output_prefix}.{output_format}")

            if output_format == "tsv":
                filtered_df.write_csv(output_path, separator="\t")
            elif output_format == "tsv.gz":
                csv_content = filtered_df.write_csv(separator="\t")
                with gzip.open(output_path, "wt") as f:
                    f.write(csv_content)
            elif output_format == "parquet":
                filtered_df.write_parquet(output_path)

            click.echo(f"Formatted data written to {output_path}", err=True)
    else:
        # Apply impact filters and create multiple output files
        if not output_prefix:
            click.echo(
                "Error: Output prefix required when using impact filters.",
                err=True,
            )
            raise click.Abort()

        apply_impact_filters(
            filtered_df,
            impact_config,
            output_prefix,
            output_format,
            verbose,
        )


def read_pedigree(pedigree_path: Path) -> pl.DataFrame:
    """
    Read a pedigree file and return a DataFrame with sample relationships.

    Args:
        pedigree_path: Path to the pedigree file

    Returns:
        DataFrame with columns: sample_id, father_id, mother_id
    """
    # Try reading with header first
    # Force certain columns to string type
    string_columns = [
        "FID",
        "sample_id",
        "father_id",
        "mother_id",
        "FatherBarcode",
        "MotherBarcode",
        "sample",
    ]
    schema_overrides = {col: pl.Utf8 for col in string_columns}
    df = pl.read_csv(pedigree_path, separator="\t", schema_overrides=schema_overrides)

    # Check if first row has 'FID' in first column (indicates header)
    if df.columns[0] == "FID" or "sample_id" in df.columns:
        # Has header - use it as-is
        pass
    else:
        # No header - assume standard pedigree format
        # FID, sample_id, father_id, mother_id, sex, phenotype
        df.columns = ["FID", "sample_id", "father_id", "mother_id", "sex", "phenotype"]

    # Ensure we have the required columns (try different possible names)
    if "sample_id" not in df.columns and len(df.columns) >= 4:
        # Try to identify columns by position
        df = df.rename(
            {
                df.columns[1]: "sample_id",
                df.columns[2]: "father_id",
                df.columns[3]: "mother_id",
            }
        )

    # Handle different column names for father/mother
    if "FatherBarcode" in df.columns:
        df = df.rename({"FatherBarcode": "father_id", "MotherBarcode": "mother_id"})

    # Normalize sex column name if present (e.g., 'Sex' or 'sex')
    sex_col = next((c for c in df.columns if c.lower() == "sex"), None)
    if sex_col and sex_col != "sex":
        df = df.rename({sex_col: "sex"})

    # Select only the columns we need (include sex if present)
    select_cols = ["sample_id", "father_id", "mother_id"]
    if "sex" in df.columns:
        select_cols.append("sex")
    pedigree_df = df.select(select_cols)

    # Replace 0 and -9 with null (indicating no parent)
    # Explicit cast to Utf8 ensures type is preserved even when all values become null
    pedigree_df = pedigree_df.with_columns(
        [
            pl.when(pl.col("father_id").cast(pl.Utf8).is_in(["0", "-9"]))
            .then(None)
            .otherwise(pl.col("father_id"))
            .cast(pl.Utf8)
            .alias("father_id"),
            pl.when(pl.col("mother_id").cast(pl.Utf8).is_in(["0", "-9"]))
            .then(None)
            .otherwise(pl.col("mother_id"))
            .cast(pl.Utf8)
            .alias("mother_id"),
        ]
    )

    return pedigree_df


def add_parent_genotypes(df: pl.DataFrame, pedigree_df: pl.DataFrame) -> pl.DataFrame:
    """
    Add father and mother genotype columns to the DataFrame.

    Args:
        df: DataFrame with sample genotype information
        pedigree_df: DataFrame with parent relationships

    Returns:
        DataFrame with added parent genotype columns
    """
    # Join with pedigree to get father and mother IDs for each sample
    df = df.join(pedigree_df, left_on="sample", right_on="sample_id", how="left")

    # Define the core variant-identifying columns for joining parent genotypes
    # We only want to join on genomic position, not on annotation columns
    # This ensures we match parents even if they have different VEP annotations
    core_variant_cols = ["#CHROM", "POS", "REF", "ALT"]
    # Check which columns actually exist in the dataframe
    join_cols = [col for col in core_variant_cols if col in df.columns]

    # Create a self-join friendly version of the data for looking up parent genotypes
    # We select only the join columns + sample genotype information
    parent_lookup = df.select(
        join_cols
        + [
            pl.col("sample"),
            pl.col("sample_gt"),
            pl.col("sample_dp"),
            pl.col("sample_gq"),
            pl.col("sample_ad"),
            pl.col("sample_vaf"),
        ]
    ).unique()

    # Join for father's genotypes
    # Match on genomic position AND father_id == sample
    father_data = parent_lookup.rename(
        {
            "sample": "father_id",
            "sample_gt": "father_gt",
            "sample_dp": "father_dp",
            "sample_gq": "father_gq",
            "sample_ad": "father_ad",
            "sample_vaf": "father_vaf",
        }
    )

    df = df.join(father_data, on=join_cols + ["father_id"], how="left")

    # Join for mother's genotypes
    mother_data = parent_lookup.rename(
        {
            "sample": "mother_id",
            "sample_gt": "mother_gt",
            "sample_dp": "mother_dp",
            "sample_gq": "mother_gq",
            "sample_ad": "mother_ad",
            "sample_vaf": "mother_vaf",
        }
    )

    df = df.join(mother_data, on=join_cols + ["mother_id"], how="left")

    # Rename father_id and mother_id to father and mother for debugging
    df = df.rename({"father_id": "father", "mother_id": "mother"})

    # Replace '.' with '0' for parent DP and GQ columns
    df = df.with_columns(
        [
            pl.when(pl.col("father_dp") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("father_dp"))
            .alias("father_dp"),
            pl.when(pl.col("father_gq") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("father_gq"))
            .alias("father_gq"),
            pl.when(pl.col("mother_dp") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("mother_dp"))
            .alias("mother_dp"),
            pl.when(pl.col("mother_gq") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("mother_gq"))
            .alias("mother_gq"),
        ]
    )

    return df


def format_expand_annotations(df: pl.DataFrame) -> pl.DataFrame:
    """
    Expand the (null) annotation column into separate columns.

    This is a separate step that can be applied after filtering to avoid
    expensive annotation expansion on variants that will be filtered out.

    Handles two types of INFO fields:
    - Key-value pairs (e.g., "DP=30") -> extracted as string values
    - Boolean flags (e.g., "PASS", "DB") -> created as True/False columns

    Args:
        df: DataFrame with (null) column

    Returns:
        DataFrame with expanded annotation columns
    """
    # Find the (null) column
    if "(null)" not in df.columns:
        # Already expanded or missing - return as-is
        return df

    # Extract all unique field names and flags from the (null) column
    null_values = df.select("(null)").to_series()
    all_fields = set()
    all_flags = set()

    for value in null_values:
        if value and not (isinstance(value, float)):  # Skip null/NaN values
            pairs = str(value).split(";")
            for pair in pairs:
                if "=" in pair:
                    field_name = pair.split("=", 1)[0]
                    all_fields.add(field_name)
                elif pair.strip():  # Boolean flag (no '=')
                    all_flags.add(pair.strip())

    # Create expressions to extract each key-value field
    for field in sorted(all_fields):
        # Extract the field value from the (null) column
        # Pattern: extract value after "field=" and before ";" or end of string
        df = df.with_columns(
            pl.col("(null)").str.extract(f"{field}=([^;]+)").alias(field)
        )

    # Create boolean columns for flags
    for flag in sorted(all_flags):
        # Check if flag appears in the (null) column (as whole word)
        # Use regex to match flag as a separate field (not part of another field name)
        df = df.with_columns(
            pl.col("(null)").str.contains(f"(^|;){flag}(;|$)").alias(flag)
        )

    # Drop the original (null) column
    df = df.drop("(null)")

    # Drop CSQ column if it exists (it was extracted from (null) column)
    if "CSQ" in df.columns:
        df = df.drop("CSQ")

    return df


def format_bcftools_tsv_minimal(
    df: pl.DataFrame, pedigree_df: Optional[pl.DataFrame] = None
) -> pl.DataFrame:
    """
    Format a bcftools tabulated TSV with minimal processing.

    This version SKIPS expanding the (null) annotation field and only:
    - Melts sample columns into rows
    - Extracts GT:DP:GQ:AD:VAF from sample values
    - Optionally adds parent genotypes if pedigree provided

    Use this for filtering workflows where you want to apply filters BEFORE
    expensive annotation expansion. Call format_expand_annotations() afterwards
    on filtered results.

    Args:
        df: Input DataFrame from bcftools
        pedigree_df: Optional pedigree DataFrame with parent information

    Returns:
        Formatted DataFrame with melted samples (annotations still in (null) column)
    """
    # Determine which columns are sample columns
    # Sample columns have format "SampleName:GT:SampleName:DP:..." or similar
    # Non-sample columns are standard VCF columns or annotation columns

    # Standard VCF/annotation columns (not samples)
    standard_cols = {
        "#CHROM", "POS", "REF", "ALT", "FILTER", "(null)", "CSQ",
        "QUAL", "ID", "INFO", "FORMAT"
    }

    # Find sample columns by looking for columns with ":" in the name
    # that aren't standard columns
    sample_cols = []
    sample_names = []

    for col in df.columns:
        # Skip standard columns
        if col in standard_cols:
            continue

        # Skip columns that look like VEP annotation fields
        if col.startswith("VEP_") or col.startswith("AF") or col.startswith("AC"):
            continue

        # Sample columns typically have ":" in them (GT:DP:GQ:AD format)
        if ":" in col:
            sample_name = col.split(":", 1)[0]
            sample_cols.append(col)
            sample_names.append(sample_name)
        elif col not in df.columns[:10]:
            # Columns after position 10 that don't match known patterns might be samples
            # This is a heuristic for unusual sample column formats
            pass

    if not sample_cols:
        # No sample columns to melt
        return df

    # Step 2: Melt the sample columns
    # Keep all columns except sample columns as id_vars
    id_vars = [col for col in df.columns if col not in sample_cols]

    # Create a mapping of old column names to sample names
    rename_map = {old: new for old, new in zip(sample_cols, sample_names)}

    # Rename sample columns to just sample names before melting
    df = df.rename(rename_map)

    # Melt the dataframe
    melted_df = df.melt(
        id_vars=id_vars,
        value_vars=sample_names,
        variable_name="sample",
        value_name="sample_value",
    )

    # Step 4: Split sample_value into GT:DP:GQ:AD format
    # Split on ':' to get individual fields
    # Use nullable=True to handle missing fields gracefully
    melted_df = melted_df.with_columns(
        [
            # GT - first field (nullable for missing data)
            pl.col("sample_value")
            .str.split(":")
            .list.get(0, null_on_oob=True)
            .alias("sample_gt"),
            # DP - second field (nullable for missing data)
            pl.col("sample_value")
            .str.split(":")
            .list.get(1, null_on_oob=True)
            .alias("sample_dp"),
            # GQ - third field (nullable for missing data)
            pl.col("sample_value")
            .str.split(":")
            .list.get(2, null_on_oob=True)
            .alias("sample_gq"),
            # AD - fourth field, split on ',' and keep second value (nullable)
            pl.col("sample_value")
            .str.split(":")
            .list.get(3, null_on_oob=True)
            .str.split(",")
            .list.get(1, null_on_oob=True)
            .alias("sample_ad"),
        ]
    )

    # Replace '.' with '0' for DP and GQ columns
    melted_df = melted_df.with_columns(
        [
            pl.when(pl.col("sample_dp") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("sample_dp"))
            .alias("sample_dp"),
            pl.when(pl.col("sample_gq") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("sample_gq"))
            .alias("sample_gq"),
        ]
    )

    # Step 5: Calculate sample_vaf as sample_ad / sample_dp
    # Convert to numeric, calculate ratio, handle division by zero
    melted_df = melted_df.with_columns(
        [
            (
                pl.col("sample_ad").cast(pl.Float64, strict=False)
                / pl.col("sample_dp").cast(pl.Float64, strict=False)
            ).alias("sample_vaf")
        ]
    )

    # Drop the original sample_value column
    melted_df = melted_df.drop("sample_value")

    # Step 6: Add parent genotype information if pedigree is provided
    if pedigree_df is not None:
        melted_df = add_parent_genotypes(melted_df, pedigree_df)

    return melted_df


def format_bcftools_tsv(
    df: pl.DataFrame, pedigree_df: Optional[pl.DataFrame] = None
) -> pl.DataFrame:
    """
    Format a bcftools tabulated TSV DataFrame (full processing).

    This is the complete formatting that:
    1. Melts samples and extracts GT:DP:GQ:AD:VAF
    2. Expands (null) annotation field into separate columns
    3. Adds parent genotypes if pedigree provided

    For DNM filtering workflows, consider using format_bcftools_tsv_minimal()
    + apply_de_novo_filter() + format_expand_annotations() to avoid expanding
    annotations on variants that will be filtered out.

    Args:
        df: Input DataFrame from bcftools
        pedigree_df: Optional pedigree DataFrame with parent information

    Returns:
        Formatted DataFrame with expanded fields and melted samples
    """
    # First do minimal formatting (melt + sample columns)
    melted_df = format_bcftools_tsv_minimal(df, pedigree_df)

    # Then expand annotations
    expanded_df = format_expand_annotations(melted_df)

    return expanded_df


def format_bcftools_tsv_lazy(
    lazy_df: pl.LazyFrame, pedigree_df: Optional[pl.DataFrame] = None
) -> pl.LazyFrame:
    """
    Format a bcftools tabulated TSV using lazy operations for streaming.

    This is a simplified version that collects minimally for complex operations.
    """
    # For complex transformations like melting, we need to collect temporarily
    # but we do this in a streaming fashion
    df = lazy_df.collect(streaming=True)
    formatted_df = format_bcftools_tsv(df, pedigree_df)
    return formatted_df.lazy()


# ------------------ Chunked two-pass processing with progress ------------------
import io
import math


def _open_file_lines(path: Path):
    """Yield header line and an iterator over the remaining lines (text)."""
    if str(path).endswith(".gz"):
        import gzip as _gzip

        f = _gzip.open(path, "rt")
    else:
        f = open(path, "rt")

    try:
        header = f.readline()
        for line in f:
            yield header, line
    finally:
        f.close()


def _line_iterator(path: Path, chunk_size: int = 50000):
    """Yield (header, chunk_lines) tuples where chunk_lines is a list of lines."""
    if str(path).endswith(".gz"):
        import gzip as _gzip

        f = _gzip.open(path, "rt")
    else:
        f = open(path, "rt")

    try:
        header = f.readline()
        while True:
            chunk = []
            for _ in range(chunk_size):
                line = f.readline()
                if not line:
                    break
                chunk.append(line)
            if not chunk:
                break
            yield header, chunk
    finally:
        f.close()


def build_parent_lookup_from_file(
    path: Path,
    pedigree_df: Optional[pl.DataFrame] = None,
    progress_bar=None,
    verbose: bool = False,
    chunk_size: int = 50000,
):
    """First pass: build a minimal parent lookup DataFrame with per-sample genotypes.

    Returns a tuple (lookup_df, total_lines) where total_lines is the approximate number
    of data lines processed (excluding header). If a `progress_bar` is provided it will
    be updated as we process chunks.
    """
    parts = []

    schema_overrides = {
        col: pl.Utf8
        for col in [
            "FID",
            "sample_id",
            "father_id",
            "mother_id",
            "FatherBarcode",
            "MotherBarcode",
            "sample",
        ]
    }

    processed = 0
    chunk_idx = 0
    for header, chunk in _line_iterator(path, chunk_size=chunk_size):
        chunk_idx += 1
        content = header + "".join(chunk)
        try:
            df_chunk = pl.read_csv(
                io.StringIO(content), separator="\t", schema_overrides=schema_overrides
            )
        except Exception:
            # Skip unparsable chunk
            processed += len(chunk)
            if progress_bar is not None:
                progress_bar.update(len(chunk))
            elif verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Building lookup: processed ~{processed} lines...", err=True
                )
            continue

        try:
            # Use minimal format for lookup building (skip annotation expansion)
            formatted = format_bcftools_tsv_minimal(df_chunk, pedigree_df=None)
        except Exception:
            # If chunk cannot be parsed into variants, skip
            processed += len(chunk)
            if progress_bar is not None:
                progress_bar.update(len(chunk))
            elif verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Building lookup: processed ~{processed} lines...", err=True
                )
            continue

        cols = [
            "#CHROM",
            "POS",
            "REF",
            "ALT",
            "sample",
            "sample_gt",
            "sample_dp",
            "sample_gq",
            "sample_ad",
            "sample_vaf",
        ]
        sel = [c for c in cols if c in formatted.columns]
        # Optimization: Only store non-reference genotypes in lookup (skip 0/0)
        part = formatted.select(sel).filter(pl.col("sample_gt") != "0/0").unique()
        parts.append(part)

        # Update progress
        processed += len(chunk)
        if progress_bar is not None:
            progress_bar.update(len(chunk))
        elif verbose and (chunk_idx % 10 == 0):
            click.echo(f"Building lookup: processed ~{processed} lines...", err=True)

    if parts:
        lookup = pl.concat(parts).unique()
    else:
        lookup = pl.DataFrame([])

    # processed currently counts number of data lines seen
    return lookup, processed


def add_parent_genotypes_from_lookup(
    df: pl.DataFrame, parent_lookup: pl.DataFrame
) -> pl.DataFrame:
    """Join father/mother genotype info from the parent_lookup into df.

    Assumes df has columns: #CHROM, POS, REF, ALT, father_id, mother_id
    """
    join_cols = [c for c in ["#CHROM", "POS", "REF", "ALT"] if c in df.columns]

    if parent_lookup.is_empty():
        # Create empty parent columns
        return df

    # Prepare father lookup
    father_lookup = parent_lookup.rename(
        {
            "sample": "father",
            "sample_gt": "father_gt",
            "sample_dp": "father_dp",
            "sample_gq": "father_gq",
            "sample_ad": "father_ad",
            "sample_vaf": "father_vaf",
        }
    )

    # Left join on join_cols + ['father']
    if "father_id" in df.columns:
        df = df.join(
            father_lookup,
            left_on=join_cols + ["father_id"],
            right_on=join_cols + ["father"],
            how="left",
        )
    else:
        # No father id, attempt join on 'father' column
        df = df.join(
            father_lookup,
            on=join_cols + ["father"],
            how="left",
        )

    # Prepare mother lookup
    mother_lookup = parent_lookup.rename(
        {
            "sample": "mother",
            "sample_gt": "mother_gt",
            "sample_dp": "mother_dp",
            "sample_gq": "mother_gq",
            "sample_ad": "mother_ad",
            "sample_vaf": "mother_vaf",
        }
    )

    if "mother_id" in df.columns:
        df = df.join(
            mother_lookup,
            left_on=join_cols + ["mother_id"],
            right_on=join_cols + ["mother"],
            how="left",
        )
    else:
        df = df.join(
            mother_lookup,
            on=join_cols + ["mother"],
            how="left",
        )

    # Normalize '.' to '0' for DP/GQ like previous function
    df = df.with_columns(
        [
            pl.when(pl.col("father_dp") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("father_dp"))
            .alias("father_dp"),
            pl.when(pl.col("father_gq") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("father_gq"))
            .alias("father_gq"),
            pl.when(pl.col("mother_dp") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("mother_dp"))
            .alias("mother_dp"),
            pl.when(pl.col("mother_gq") == ".")
            .then(pl.lit("0"))
            .otherwise(pl.col("mother_gq"))
            .alias("mother_gq"),
        ]
    )

    return df


def process_with_progress(
    input_path: Path,
    output_prefix: str,
    output_format: str,
    pedigree_df: Optional[pl.DataFrame],
    filter_config: Optional[dict],
    verbose: bool,
    chunk_size: int = 50000,
):
    """Process input in two passes and show a progress bar.

    Pass 1: build parent lookup
    Pass 2: process chunks, join parent genotypes, apply filters, and write incrementally
    """
    # tqdm optional
    try:
        from tqdm.auto import tqdm
    except Exception:
        tqdm = None

    # Build parent lookup in a single pass (counts lines while building lookup)
    if verbose:
        click.echo("Pass 1: building parent genotype lookup (single pass)...", err=True)

    pbar_lookup = None
    if tqdm is not None:
        # No known total yet; tqdm will show progress increasing
        pbar_lookup = tqdm(desc="Building parent lookup", unit="lines")

    parent_lookup, total_lines = build_parent_lookup_from_file(
        input_path,
        pedigree_df,
        progress_bar=pbar_lookup,
        verbose=verbose,
        chunk_size=chunk_size,
    )

    if pbar_lookup is not None:
        pbar_lookup.close()

    if verbose:
        click.echo(
            f"Parent lookup contains {parent_lookup.shape[0]} genotype entries (from ~{total_lines} lines)",
            err=True,
        )

    total_chunks = math.ceil(total_lines / chunk_size) if chunk_size > 0 else 1

    # Prepare output paths
    if output_format == "tsv":
        out_path = Path(f"{output_prefix}.tsv")
    elif output_format == "tsv.gz":
        out_path = Path(f"{output_prefix}.tsv.gz")
    else:
        # We'll write parquet parts
        out_path = Path(f"{output_prefix}")

    first_write = True

    # Iterate chunks and process
    iterator = _line_iterator(input_path, chunk_size=chunk_size)
    chunk_idx = 0

    progress_bar = None
    processed_lines = 0
    if tqdm is not None:
        progress_bar = tqdm(total=total_lines, desc="Processing variants")

    for header, chunk in iterator:
        chunk_idx += 1
        chunk_count = len(chunk)

        # Update progress at start of chunk to show we're working
        if progress_bar is not None:
            progress_bar.set_postfix_str(f"chunk {chunk_idx}")

        content = header + "".join(chunk)
        df_chunk = pl.read_csv(io.StringIO(content), separator="\t")

        # Use MINIMAL format (skip annotation expansion for now)
        try:
            melted = format_bcftools_tsv_minimal(df_chunk, pedigree_df=None)
        except Exception:
            # If parse fails for chunk, skip
            processed_lines += chunk_count
            if progress_bar is not None:
                progress_bar.update(chunk_count)
            elif verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Processed {processed_lines}/{total_lines} lines...", err=True
                )
            continue

        # Optimization #1: Early GT filtering for DNM mode - skip reference-only variants
        if filter_config and filter_config.get("dnm", {}).get("enabled", False):
            melted = melted.filter(
                pl.col("sample_gt").str.contains("1")
                | pl.col("sample_gt").str.contains("2")
            )
            if melted.shape[0] == 0:
                # All variants filtered out, skip to next chunk
                processed_lines += chunk_count
                if progress_bar is not None:
                    progress_bar.update(chunk_count)
                continue

        # Attach parent ids from pedigree
        if pedigree_df is not None and "sample" in melted.columns:
            melted = melted.join(
                pedigree_df, left_on="sample", right_on="sample_id", how="left"
            )

        # Add parent genotypes from global lookup
        melted = add_parent_genotypes_from_lookup(melted, parent_lookup)

        # Apply filters BEFORE expanding annotations (key optimization)
        if filter_config and filter_config.get("dnm", {}).get("enabled", False):
            cfg = {}
            cfg.update(filter_config.get("quality", {}))
            cfg.update(filter_config.get("dnm", {}))
            filtered = apply_de_novo_filter(
                melted, cfg, verbose=False, pedigree_df=pedigree_df
            )
        else:
            # Apply standard quality filters
            filtered = melted
            quality_cfg = filter_config.get("quality", {}) if filter_config else {}
            if quality_cfg:
                filtered = apply_quality_filters(filtered, quality_cfg, verbose=False)

            # Apply expression filter if present
            expr = filter_config.get("expression") if filter_config else None
            if expr:
                try:
                    expr_parsed = parse_impact_filter_expression(expr, filtered)
                    filtered = filtered.filter(expr_parsed)
                except Exception:
                    # If expression parsing fails on a chunk, skip applying it
                    pass

        # NOW expand annotations only for variants that passed filters
        if filtered.shape[0] > 0 and "(null)" in filtered.columns:
            if progress_bar is not None:
                progress_bar.set_postfix_str(f"expanding annotations chunk {chunk_idx}")
            filtered = format_expand_annotations(filtered)

        # Update progress after filtering (before write)
        if progress_bar is not None:
            progress_bar.set_postfix_str(f"writing chunk {chunk_idx}")

        # Write filtered chunk to file (skip if empty)
        if filtered.shape[0] > 0:
            if output_format in ("tsv", "tsv.gz"):
                csv_text = filtered.write_csv(separator="\t")
                # First write includes header; subsequent writes skip header
                if first_write:
                    write_text = csv_text
                    first_write = False
                    if output_format == "tsv.gz":
                        with gzip.open(out_path, "wt") as f:
                            f.write(write_text)
                    else:
                        with open(out_path, "wt") as f:
                            f.write(write_text)
                else:
                    # Skip header
                    tail = "\n".join(csv_text.splitlines()[1:])
                    if output_format == "tsv.gz":
                        with gzip.open(out_path, "at") as f:
                            f.write("\n" + tail)
                    else:
                        with open(out_path, "at") as f:
                            f.write("\n" + tail)
            else:
                # Parquet: write part file
                part_path = out_path.with_suffix(f".part{chunk_idx}.parquet")
                filtered.write_parquet(part_path)

        if progress_bar is not None:
            progress_bar.update(chunk_count)
            progress_bar.set_postfix_str("")  # Clear status
        else:
            processed_lines += chunk_count
            if verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Processed {processed_lines}/{total_lines} lines...", err=True
                )

    if progress_bar is not None:
        progress_bar.close()

    if verbose:
        click.echo("Processing complete.", err=True)


def apply_filters_lazy(
    lazy_df: pl.LazyFrame,
    filter_config: dict,
    verbose: bool = False,
    pedigree_df: Optional[pl.DataFrame] = None,
) -> pl.LazyFrame:
    """Apply quality and expression filters using lazy operations."""
    # If DNM mode is enabled, we need to collect and apply DNM logic
    if filter_config.get("dnm", {}).get("enabled", False):
        dnm_cfg = {}
        dnm_cfg.update(filter_config.get("quality", {}))
        dnm_cfg.update(filter_config.get("dnm", {}))

        # Collect minimally and apply DNM filter eagerly, then return lazy frame
        df = lazy_df.collect(streaming=True)
        filtered_df = apply_de_novo_filter(
            df, dnm_cfg, verbose, pedigree_df=pedigree_df
        )
        return filtered_df.lazy()

    quality_config = filter_config.get("quality", {})
    expression = filter_config.get("expression")

    # Apply quality filters
    if quality_config:
        # Filter: sample_gt must contain at least one '1' (default: true)
        filter_no_alt = quality_config.get("filter_no_alt_allele", True)
        if filter_no_alt:
            lazy_df = lazy_df.filter(
                pl.col("sample_gt").str.contains("1")
                | pl.col("sample_gt").str.contains("2")
            )

        # Apply minimum depth filter
        if "sample_dp_min" in quality_config:
            min_dp = quality_config["sample_dp_min"]
            lazy_df = lazy_df.filter(
                pl.col("sample_dp").cast(pl.Float64, strict=False) >= min_dp
            )

        # Apply minimum GQ filter
        if "sample_gq_min" in quality_config:
            min_gq = quality_config["sample_gq_min"]
            lazy_df = lazy_df.filter(
                pl.col("sample_gq").cast(pl.Float64, strict=False) >= min_gq
            )

        # VAF filters for heterozygous (0/1 or 1/0)
        if (
            "sample_vaf_het_min" in quality_config
            or "sample_vaf_het_max" in quality_config
        ):
            # Check if genotype is het (contains one '1' and one '0', no '2')
            is_het = (
                (pl.col("sample_gt").str.count_matches("1") == 1)
                & (pl.col("sample_gt").str.count_matches("0") == 1)
                & (~pl.col("sample_gt").str.contains("2"))
            )

            het_conditions = []
            if "sample_vaf_het_min" in quality_config:
                het_conditions.append(
                    pl.col("sample_vaf") >= quality_config["sample_vaf_het_min"]
                )
            if "sample_vaf_het_max" in quality_config:
                het_conditions.append(
                    pl.col("sample_vaf") <= quality_config["sample_vaf_het_max"]
                )

            if het_conditions:
                het_filter = het_conditions[0]
                for cond in het_conditions[1:]:
                    het_filter = het_filter & cond

                lazy_df = lazy_df.filter(~is_het | het_filter)

        # VAF filter for homozygous alternate (1/1)
        if "sample_vaf_homalt_min" in quality_config:
            is_homalt = pl.col("sample_gt") == "1/1"
            lazy_df = lazy_df.filter(
                ~is_homalt
                | (pl.col("sample_vaf") >= quality_config["sample_vaf_homalt_min"])
            )

        # VAF filter for homozygous reference (0/0)
        if "sample_vaf_hom_ref_max" in quality_config:
            is_hom_ref = pl.col("sample_gt") == "0/0"
            lazy_df = lazy_df.filter(
                ~is_hom_ref
                | (pl.col("sample_vaf") <= quality_config["sample_vaf_hom_ref_max"])
            )

        # Apply same filters to parents if requested
        apply_to_parents = quality_config.get("apply_to_parents", False)
        if apply_to_parents:
            # Father filters
            if "sample_dp_min" in quality_config:
                min_dp = quality_config["sample_dp_min"]
                lazy_df = lazy_df.filter(
                    (pl.col("father_dp").is_null())
                    | (pl.col("father_dp").cast(pl.Float64, strict=False) >= min_dp)
                )

            if "sample_gq_min" in quality_config:
                min_gq = quality_config["sample_gq_min"]
                lazy_df = lazy_df.filter(
                    (pl.col("father_gq").is_null())
                    | (pl.col("father_gq").cast(pl.Float64, strict=False) >= min_gq)
                )

            # Mother filters
            if "sample_dp_min" in quality_config:
                min_dp = quality_config["sample_dp_min"]
                lazy_df = lazy_df.filter(
                    (pl.col("mother_dp").is_null())
                    | (pl.col("mother_dp").cast(pl.Float64, strict=False) >= min_dp)
                )

            if "sample_gq_min" in quality_config:
                min_gq = quality_config["sample_gq_min"]
                lazy_df = lazy_df.filter(
                    (pl.col("mother_gq").is_null())
                    | (pl.col("mother_gq").cast(pl.Float64, strict=False) >= min_gq)
                )

    # Apply expression filter if provided
    if expression:
        if verbose:
            click.echo(f"Applying expression filter: {expression}", err=True)

        # We need to collect temporarily to use parse_impact_filter_expression
        df = lazy_df.collect(streaming=True)
        filter_expr = parse_impact_filter_expression(expression, df)
        lazy_df = df.lazy().filter(filter_expr)

    return lazy_df


if __name__ == "__main__":
    cli()
