"""TSV/Parquet formatting: melt samples, expand annotations."""

from typing import Optional

import polars as pl

from pywombat.pedigree import add_parent_genotypes


def format_expand_annotations(df: pl.DataFrame) -> pl.DataFrame:
    """
    Expand the (null) annotation column into separate columns.

    Handles two types of INFO fields:
    - Key-value pairs (e.g., "DP=30") -> extracted as string values
    - Boolean flags (e.g., "PASS", "DB") -> created as True/False columns

    Args:
        df: DataFrame with (null) column

    Returns:
        DataFrame with expanded annotation columns
    """
    if "(null)" not in df.columns:
        return df

    null_values = df.select("(null)").to_series()
    all_fields = set()
    all_flags = set()

    for value in null_values:
        if value and not (isinstance(value, float)):
            pairs = str(value).split(";")
            for pair in pairs:
                if "=" in pair:
                    field_name = pair.split("=", 1)[0]
                    all_fields.add(field_name)
                elif pair.strip():
                    all_flags.add(pair.strip())

    for field in sorted(all_fields):
        df = df.with_columns(
            pl.col("(null)").str.extract(f"{field}=([^;]+)").alias(field)
        )

    for flag in sorted(all_flags):
        df = df.with_columns(
            pl.col("(null)").str.contains(f"(^|;){flag}(;|$)").alias(flag)
        )

    df = df.drop("(null)")

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

    Args:
        df: Input DataFrame from bcftools
        pedigree_df: Optional pedigree DataFrame with parent information

    Returns:
        Formatted DataFrame with melted samples (annotations still in (null) column)
    """
    standard_cols = {
        "#CHROM", "POS", "REF", "ALT", "FILTER", "(null)", "CSQ",
        "QUAL", "ID", "INFO", "FORMAT"
    }

    sample_cols = []
    sample_names = []

    for col in df.columns:
        if col in standard_cols:
            continue
        if col.startswith("VEP_") or col.startswith("AF") or col.startswith("AC"):
            continue
        if ":" in col:
            sample_name = col.split(":", 1)[0]
            sample_cols.append(col)
            sample_names.append(sample_name)
        elif col not in df.columns[:10]:
            pass

    if not sample_cols:
        return df

    id_vars = [col for col in df.columns if col not in sample_cols]
    rename_map = {old: new for old, new in zip(sample_cols, sample_names)}
    df = df.rename(rename_map)

    melted_df = df.unpivot(
        index=id_vars,
        on=sample_names,
        variable_name="sample",
        value_name="sample_value",
    )

    melted_df = melted_df.with_columns(
        [
            pl.col("sample_value")
            .str.split(":")
            .list.get(0, null_on_oob=True)
            .alias("sample_gt"),
            pl.col("sample_value")
            .str.split(":")
            .list.get(1, null_on_oob=True)
            .alias("sample_dp"),
            pl.col("sample_value")
            .str.split(":")
            .list.get(2, null_on_oob=True)
            .alias("sample_gq"),
            pl.col("sample_value")
            .str.split(":")
            .list.get(3, null_on_oob=True)
            .str.split(",")
            .list.get(1, null_on_oob=True)
            .alias("sample_ad"),
        ]
    )

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

    melted_df = melted_df.with_columns(
        [
            (
                pl.col("sample_ad").cast(pl.Float64, strict=False)
                / pl.col("sample_dp").cast(pl.Float64, strict=False)
            ).alias("sample_vaf")
        ]
    )

    melted_df = melted_df.drop("sample_value")

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

    Args:
        df: Input DataFrame from bcftools
        pedigree_df: Optional pedigree DataFrame with parent information

    Returns:
        Formatted DataFrame with expanded fields and melted samples
    """
    melted_df = format_bcftools_tsv_minimal(df, pedigree_df)
    expanded_df = format_expand_annotations(melted_df)
    return expanded_df


def format_bcftools_tsv_lazy(
    lazy_df: pl.LazyFrame, pedigree_df: Optional[pl.DataFrame] = None
) -> pl.LazyFrame:
    """
    Format a bcftools tabulated TSV using lazy operations for streaming.

    This is a simplified version that collects minimally for complex operations.
    """
    df = lazy_df.collect(engine="streaming")
    formatted_df = format_bcftools_tsv(df, pedigree_df)
    return formatted_df.lazy()
