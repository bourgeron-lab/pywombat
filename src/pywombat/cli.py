"""CLI for wombat tool."""

from pathlib import Path
from typing import Optional

import click
import polars as pl


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
    type=click.Choice(["tsv", "parquet"], case_sensitive=False),
    default="tsv",
    help="Output format: tsv (default) or parquet.",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
@click.option(
    "-p",
    "--pedigree",
    type=click.Path(exists=True, path_type=Path),
    help="Pedigree file to add father and mother genotype columns.",
)
def cli(
    input_file: Path,
    output: Optional[str],
    output_format: str,
    verbose: bool,
    pedigree: Optional[Path],
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

        # Read the TSV file
        df = pl.read_csv(input_file, separator="\t")

        if verbose:
            click.echo(
                f"Input shape: {df.shape[0]} rows, {df.shape[1]} columns", err=True
            )

        # Read pedigree file if provided
        pedigree_df = None
        if pedigree:
            if verbose:
                click.echo(f"Reading pedigree file: {pedigree}", err=True)
            pedigree_df = read_pedigree(pedigree)

        # Process the dataframe
        formatted_df = format_bcftools_tsv(df, pedigree_df)

        if verbose:
            click.echo(
                f"Output shape: {formatted_df.shape[0]} rows, {formatted_df.shape[1]} columns",
                err=True,
            )

        # Output the result
        if output:
            # Construct output filename with prefix and format
            output_path = Path(f"{output}.{output_format}")

            if output_format == "tsv":
                formatted_df.write_csv(output_path, separator="\t")
            elif output_format == "parquet":
                formatted_df.write_parquet(output_path)

            click.echo(f"Formatted data written to {output_path}", err=True)
        else:
            # Write to stdout (only for TSV format)
            if output_format != "tsv":
                click.echo(
                    "Error: stdout output only supported for TSV format. Use -o to specify an output prefix for parquet.",
                    err=True,
                )
                raise click.Abort()
            click.echo(formatted_df.write_csv(separator="\t"), nl=False)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


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


if __name__ == "__main__":
    cli()
