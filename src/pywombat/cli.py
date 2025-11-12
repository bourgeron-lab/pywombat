"""CLI for wombat tool."""

from pathlib import Path
from typing import Optional

import click
import polars as pl


@click.group()
def cli():
    """Wombat: A tool for processing bcftools tabulated TSV files."""
    pass


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file path. If not specified, prints to stdout.",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def format(input_file: Path, output: Optional[Path], verbose: bool):
    """
    Format a bcftools tabulated TSV file.

    This command:

    \b
    1. Expands the '(null)' column containing NAME=value pairs separated by ';'
    2. Preserves the CSQ (Consequence) column without melting
    3. Melts sample columns into rows with sample names and values

    \b
    Example:
        wombat format input.tsv -o output.tsv
        wombat format input.tsv > output.tsv
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

        # Process the dataframe
        formatted_df = format_bcftools_tsv(df)

        if verbose:
            click.echo(
                f"Output shape: {formatted_df.shape[0]} rows, {formatted_df.shape[1]} columns",
                err=True,
            )

        # Output the result
        if output:
            formatted_df.write_csv(output, separator="\t")
            click.echo(f"Formatted data written to {output}", err=True)
        else:
            # Write to stdout
            click.echo(formatted_df.write_csv(separator="\t"), nl=False)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


def format_bcftools_tsv(df: pl.DataFrame) -> pl.DataFrame:
    """
    Format a bcftools tabulated TSV DataFrame.

    Args:
        df: Input DataFrame from bcftools

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

    # Step 2: Identify sample columns and extract sample names
    # Sample columns have format "sample_name:..." in the header
    # Skip the CSQ column as it should not be melted
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

    return melted_df


if __name__ == "__main__":
    cli()
