"""CLI for wombat tool."""

import gzip
import re
import warnings
from pathlib import Path
from typing import Optional

import click
import polars as pl
import yaml


@click.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=str,
    help="Output file prefix. If not specified, prints to stdout.",
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
    help="Debug mode: show rows matching chrom:pos (e.g., chr11:70486013). Displays #CHROM, POS, VEP_SYMBOL, and columns from filter expression.",
)
def cli(
    input_file: Path,
    output: Optional[str],
    output_format: str,
    verbose: bool,
    pedigree: Optional[Path],
    filter_config: Optional[Path],
    debug: Optional[str],
):
    """
    Wombat: A tool for processing bcftools tabulated TSV files.

    This command:

    \b
    1. Expands the '(null)' column containing NAME=value pairs separated by ';'
    2. Preserves the CSQ (Consequence) column without melting
    3. Melts sample columns into rows with sample names
    4. Splits sample values (GT:DP:GQ:AD format) into separate columns:
       - sample_gt: Genotype
       - sample_dp: Read depth
       - sample_gq: Genotype quality
       - sample_ad: Allele depth (second value from comma-separated list)
       - sample_vaf: Variant allele frequency (sample_ad / sample_dp)

    \b
    Examples:
        wombat input.tsv -o output
        wombat input.tsv -o output -f parquet
        wombat input.tsv > output.tsv
    """
    try:
        if verbose:
            click.echo(f"Reading input file: {input_file}", err=True)

        # Detect if file is gzipped based on extension
        is_gzipped = str(input_file).endswith(".gz")

        if verbose and is_gzipped:
            click.echo("Detected gzipped file", err=True)

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
            # Remove .tsv.gz or .tsv extension
            if input_stem.endswith(".tsv.gz"):
                input_stem = input_stem[:-7]  # Remove .tsv.gz
            elif input_stem.endswith(".tsv"):
                input_stem = input_stem[:-4]  # Remove .tsv

            # Add config name if filter is provided
            if filter_config:
                config_name = filter_config.stem  # Get basename without extension
                output = f"{input_stem}.{config_name}"
            else:
                output = input_stem

        # Use streaming approach with lazy API
        if verbose:
            click.echo("Processing with streaming mode...", err=True)

        # Build lazy query
        lazy_df = pl.scan_csv(input_file, separator="\t")

        # Apply formatting transformations
        lazy_df = format_bcftools_tsv_lazy(lazy_df, pedigree_df)

        # Apply filters if provided
        if filter_config_data:
            lazy_df = apply_filters_lazy(lazy_df, filter_config_data, verbose)

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
    df = pl.read_csv(input_file, separator="\t")
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
):
    """Apply filters and write output files."""
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
    df = pl.read_csv(pedigree_path, separator="\t")

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

    # Select only the columns we need
    pedigree_df = df.select(["sample_id", "father_id", "mother_id"])

    # Replace 0 and -9 with null (indicating no parent)
    pedigree_df = pedigree_df.with_columns(
        [
            pl.when(pl.col("father_id").cast(pl.Utf8).is_in(["0", "-9"]))
            .then(None)
            .otherwise(pl.col("father_id"))
            .alias("father_id"),
            pl.when(pl.col("mother_id").cast(pl.Utf8).is_in(["0", "-9"]))
            .then(None)
            .otherwise(pl.col("mother_id"))
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


def format_bcftools_tsv(
    df: pl.DataFrame, pedigree_df: Optional[pl.DataFrame] = None
) -> pl.DataFrame:
    """
    Format a bcftools tabulated TSV DataFrame.

    Args:
        df: Input DataFrame from bcftools
        pedigree_df: Optional pedigree DataFrame with parent information

    Returns:
        Formatted DataFrame with expanded fields and melted samples
    """
    # Find the (null) column
    if "(null)" not in df.columns:
        raise ValueError("Column '(null)' not found in the input file")

    # Get column index of (null)
    null_col_idx = df.columns.index("(null)")

    # Split columns into: before (null), (null), and after (null)
    cols_after = df.columns[null_col_idx + 1 :]

    # Step 1: Expand the (null) column
    # Split by semicolon and create new columns

    # First, we need to extract all unique field names from the (null) column
    # to know what columns to create
    null_values = df.select("(null)").to_series()
    all_fields = set()

    for value in null_values:
        if value and not (isinstance(value, float)):  # Skip null/NaN values
            pairs = str(value).split(";")
            for pair in pairs:
                if "=" in pair:
                    field_name = pair.split("=", 1)[0]
                    all_fields.add(field_name)

    # Create expressions to extract each field
    for field in sorted(all_fields):
        # Extract the field value from the (null) column
        # Pattern: extract value after "field=" and before ";" or end of string
        df = df.with_columns(
            pl.col("(null)").str.extract(f"{field}=([^;]+)").alias(field)
        )

    # Drop the original (null) column
    df = df.drop("(null)")

    # Drop CSQ column if it exists (it was extracted from (null) column)
    if "CSQ" in df.columns:
        df = df.drop("CSQ")

    # Step 2: Identify sample columns and extract sample names
    # Sample columns have format "sample_name:..." in the header
    # Skip the CSQ column as it should not be melted (handled above)
    sample_cols = []
    sample_names = []

    for col in cols_after:
        # Skip CSQ column
        if col == "CSQ":
            continue

        if ":" in col:
            sample_name = col.split(":", 1)[0]
            sample_cols.append(col)
            sample_names.append(sample_name)
        else:
            # If no colon, treat the whole column name as sample name
            sample_cols.append(col)
            sample_names.append(col)

    if not sample_cols:
        # No sample columns to melt, just return expanded data
        return df

    # Step 3: Melt the sample columns
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


def apply_filters_lazy(
    lazy_df: pl.LazyFrame, filter_config: dict, verbose: bool = False
) -> pl.LazyFrame:
    """Apply quality and expression filters using lazy operations."""
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
