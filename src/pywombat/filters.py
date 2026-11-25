"""Quality and expression filters, lazy filter pipeline."""

import warnings
from pathlib import Path
from typing import Optional

import click
import polars as pl
import yaml

from pywombat.denovo import apply_de_novo_filter
from pywombat.expression import parse_impact_filter_expression


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

    if "2" in str(df["sample_gt"].to_list()):
        rare_gts = df.filter(pl.col("sample_gt").str.contains("2"))
        if rare_gts.shape[0] > 0:
            warnings.warn(
                f"Found {rare_gts.shape[0]} rows with rare genotypes containing '2'. These will be kept."
            )

    filter_no_alt = quality_config.get("filter_no_alt_allele", True)
    if filter_no_alt:
        df = df.filter(
            pl.col("sample_gt").str.contains("1")
            | pl.col("sample_gt").str.contains("2")
        )

    homalt_only = quality_config.get("homalt_only", False)
    if homalt_only:
        df = df.filter(pl.col("sample_gt") == "1/1")
        if verbose:
            click.echo(
                f"Quality: {df.shape[0]} variants after filtering for homozygous alternative (1/1) only",
                err=True,
            )

    if "sample_dp_min" in quality_config:
        min_dp = quality_config["sample_dp_min"]
        df = df.filter(pl.col("sample_dp").cast(pl.Float64, strict=False) >= min_dp)

    if "sample_gq_min" in quality_config:
        min_gq = quality_config["sample_gq_min"]
        df = df.filter(pl.col("sample_gq").cast(pl.Float64, strict=False) >= min_gq)

    is_het = (pl.col("sample_gt").str.count_matches("1") == 1) & ~pl.col(
        "sample_gt"
    ).str.contains("2")
    is_hom_alt = pl.col("sample_gt") == "1/1"
    is_hom_ref = pl.col("sample_gt") == "0/0"

    if "sample_vaf_het_min" in quality_config:
        min_vaf_het = quality_config["sample_vaf_het_min"]
        df = df.filter(~is_het | (pl.col("sample_vaf") >= min_vaf_het))

    if "sample_vaf_het_max" in quality_config:
        max_vaf_het = quality_config["sample_vaf_het_max"]
        df = df.filter(~is_het | (pl.col("sample_vaf") <= max_vaf_het))

    if "sample_vaf_homalt_min" in quality_config:
        min_vaf_homalt = quality_config["sample_vaf_homalt_min"]
        df = df.filter(~is_hom_alt | (pl.col("sample_vaf") >= min_vaf_homalt))

    if "sample_vaf_hom_ref_max" in quality_config:
        max_vaf_hom_ref = quality_config["sample_vaf_hom_ref_max"]
        df = df.filter(~is_hom_ref | (pl.col("sample_vaf") <= max_vaf_hom_ref))

    apply_to_parents = quality_config.get("apply_to_parents", False)
    if apply_to_parents and "father" in df.columns:
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


def apply_filters_lazy(
    lazy_df: pl.LazyFrame,
    filter_config: dict,
    verbose: bool = False,
    pedigree_df: Optional[pl.DataFrame] = None,
) -> pl.LazyFrame:
    """Apply quality and expression filters using lazy operations."""
    if filter_config.get("dnm", {}).get("enabled", False):
        dnm_cfg = {}
        dnm_cfg.update(filter_config.get("quality", {}))
        dnm_cfg.update(filter_config.get("dnm", {}))

        df = lazy_df.collect(engine="streaming")
        filtered_df = apply_de_novo_filter(
            df, dnm_cfg, verbose, pedigree_df=pedigree_df
        )
        return filtered_df.lazy()

    quality_config = filter_config.get("quality", {})
    expression = filter_config.get("expression")

    if quality_config:
        filter_no_alt = quality_config.get("filter_no_alt_allele", True)
        if filter_no_alt:
            lazy_df = lazy_df.filter(
                pl.col("sample_gt").str.contains("1")
                | pl.col("sample_gt").str.contains("2")
            )

        homalt_only = quality_config.get("homalt_only", False)
        if homalt_only:
            lazy_df = lazy_df.filter(pl.col("sample_gt") == "1/1")

        if "sample_dp_min" in quality_config:
            min_dp = quality_config["sample_dp_min"]
            lazy_df = lazy_df.filter(
                pl.col("sample_dp").cast(pl.Float64, strict=False) >= min_dp
            )

        if "sample_gq_min" in quality_config:
            min_gq = quality_config["sample_gq_min"]
            lazy_df = lazy_df.filter(
                pl.col("sample_gq").cast(pl.Float64, strict=False) >= min_gq
            )

        if (
            "sample_vaf_het_min" in quality_config
            or "sample_vaf_het_max" in quality_config
        ):
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

        if "sample_vaf_homalt_min" in quality_config:
            is_homalt = pl.col("sample_gt") == "1/1"
            lazy_df = lazy_df.filter(
                ~is_homalt
                | (pl.col("sample_vaf") >= quality_config["sample_vaf_homalt_min"])
            )

        if "sample_vaf_hom_ref_max" in quality_config:
            is_hom_ref = pl.col("sample_gt") == "0/0"
            lazy_df = lazy_df.filter(
                ~is_hom_ref
                | (pl.col("sample_vaf") <= quality_config["sample_vaf_hom_ref_max"])
            )

        apply_to_parents = quality_config.get("apply_to_parents", False)
        if apply_to_parents:
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

    if expression:
        if verbose:
            click.echo(f"Applying expression filter: {expression}", err=True)

        df = lazy_df.collect(engine="streaming")
        filter_expr = parse_impact_filter_expression(expression, df)
        lazy_df = df.lazy().filter(filter_expr)

    return lazy_df
