"""Pedigree reading and parent genotype handling."""

from pathlib import Path
from typing import Optional

import click
import polars as pl


_PEDIGREE_HEADER_NAMES = {
    "fid",
    "iid",
    "sample_id",
    "father_id",
    "mother_id",
    "fatherbarcode",
    "motherbarcode",
    "sex",
}


def _first_pedigree_line(pedigree_path: Path) -> str:
    """The first non-blank line of a pedigree file, without its line ending."""
    with open(pedigree_path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.strip():
                return line
    return ""


def _pedigree_has_header(pedigree_path: Path) -> bool:
    """
    Whether the pedigree file's first line is a header rather than a member.

    A pedigree may come with or without a header. Reading a headerless one as if it
    had a header loses its first member, whose parents and sex then never reach the
    output (so a child listed first gets no de novo calls).
    """
    fields = _first_pedigree_line(pedigree_path).lstrip("#").split("\t")
    return any(f.strip().lower() in _PEDIGREE_HEADER_NAMES for f in fields)


def read_pedigree(pedigree_path: Path) -> pl.DataFrame:
    """
    Read a pedigree file and return a DataFrame with sample relationships.

    Args:
        pedigree_path: Path to the pedigree file

    Returns:
        DataFrame with columns: sample_id, father_id, mother_id
    """
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

    if _pedigree_has_header(pedigree_path):
        df = pl.read_csv(
            pedigree_path, separator="\t", schema_overrides=schema_overrides
        )
    else:
        n_cols = len(_first_pedigree_line(pedigree_path).split("\t"))
        names = ["FID", "sample_id", "father_id", "mother_id", "sex", "phenotype"]
        names += [f"column_{i + 1}" for i in range(len(names), n_cols)]
        names = names[:n_cols]
        df = pl.read_csv(
            pedigree_path,
            separator="\t",
            has_header=False,
            new_columns=names,
            schema_overrides={c: pl.Utf8 for c in names[:4]},
        )

    if "sample_id" not in df.columns and len(df.columns) >= 4:
        df = df.rename(
            {
                df.columns[1]: "sample_id",
                df.columns[2]: "father_id",
                df.columns[3]: "mother_id",
            }
        )

    if "FatherBarcode" in df.columns:
        df = df.rename({"FatherBarcode": "father_id", "MotherBarcode": "mother_id"})

    sex_col = next((c for c in df.columns if c.lower() == "sex"), None)
    if sex_col and sex_col != "sex":
        df = df.rename({sex_col: "sex"})

    select_cols = ["sample_id", "father_id", "mother_id"]
    if "sex" in df.columns:
        select_cols.append("sex")
    pedigree_df = df.select(select_cols)

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
    df = df.join(pedigree_df, left_on="sample", right_on="sample_id", how="left")

    core_variant_cols = ["#CHROM", "POS", "REF", "ALT"]
    join_cols = [col for col in core_variant_cols if col in df.columns]

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

    df = df.rename({"father_id": "father", "mother_id": "mother"})

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


def add_parent_genotypes_from_lookup(
    df: pl.DataFrame, parent_lookup: pl.DataFrame
) -> pl.DataFrame:
    """Join father/mother genotype info from the parent_lookup into df.

    Assumes df has columns: #CHROM, POS, REF, ALT, father_id, mother_id
    """
    join_cols = [c for c in ["#CHROM", "POS", "REF", "ALT"] if c in df.columns]

    if parent_lookup.is_empty():
        return df

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

    if "father_id" in df.columns:
        df = df.join(
            father_lookup,
            left_on=join_cols + ["father_id"],
            right_on=join_cols + ["father"],
            how="left",
        )
    else:
        df = df.join(
            father_lookup,
            on=join_cols + ["father"],
            how="left",
        )

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
