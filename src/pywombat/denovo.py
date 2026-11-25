"""De novo mutation (DNM) detection and filtering."""

import warnings
from pathlib import Path
from typing import Optional

import click
import polars as pl

from pywombat.expression import parse_impact_filter_expression
from pywombat.formatting import format_bcftools_tsv_minimal


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
            if region_chrom.upper() == chrom_short.upper() and start <= pos <= end:
                return True
    return False


def get_unique_chromosomes(parquet_file: Path) -> list[str]:
    """Get list of unique chromosomes from Parquet file, sorted naturally.

    Args:
        parquet_file: Path to Parquet file

    Returns:
        Sorted list of chromosome names (e.g., ['1', '2', ..., '22', 'X', 'Y', 'MT'])
    """
    df = pl.scan_parquet(parquet_file).select("#CHROM").unique().collect()
    chroms = df["#CHROM"].to_list()

    def chrom_sort_key(chrom: str) -> tuple:
        """Sort key for natural chromosome ordering."""
        chrom_norm = chrom.replace("chr", "").replace("Chr", "").replace("CHR", "").upper()

        try:
            return (0, int(chrom_norm), "")
        except ValueError:
            pass

        if chrom_norm in ["X", "Y", "MT", "M"]:
            order = {"X": 23, "Y": 24, "MT": 25, "M": 25}
            return (1, order.get(chrom_norm, 99), chrom_norm)

        return (2, 0, chrom_norm)

    return sorted(chroms, key=chrom_sort_key)


def apply_dnm_prefilters(
    lazy_df: pl.LazyFrame,
    filter_config: dict,
    verbose: bool = False
) -> pl.LazyFrame:
    """Apply variant-level DNM filters before melting.

    These filters don't require sample-level data and can be applied
    on wide-format data to reduce memory usage.

    Applies:
    - Population frequency filters (fafmax_faf95_max_genomes_max)
    - Quality filters (genomes_filters PASS only)

    Args:
        lazy_df: LazyFrame with wide-format data (not melted)
        filter_config: Filter configuration dict
        verbose: Whether to print progress messages

    Returns:
        Filtered LazyFrame
    """
    dnm_config = filter_config.get("dnm", {})

    fafmax_max = dnm_config.get("fafmax_faf95_max_genomes_max")
    if fafmax_max is not None:
        lazy_df = lazy_df.filter(
            (pl.col("fafmax_faf95_max_genomes").cast(pl.Float64, strict=False) <= fafmax_max)
            | pl.col("fafmax_faf95_max_genomes").is_null()
        )
        if verbose:
            click.echo(
                f"DNM prefilter: Applied frequency filter (fafmax <= {fafmax_max})", err=True
            )

    if dnm_config.get("genomes_filters_pass_only", False):
        lazy_df = lazy_df.filter(
            (pl.col("genomes_filters") == ".") | pl.col("genomes_filters").is_null()
        )
        if verbose:
            click.echo(
                "DNM prefilter: Applied genomes_filters PASS filter", err=True
            )

    return lazy_df


def apply_de_novo_filter(
    df: pl.DataFrame,
    dnm_config: dict,
    verbose: bool = False,
    pedigree_df: Optional[pl.DataFrame] = None,
    skip_prefilters: bool = False,
) -> pl.DataFrame:
    """Apply de novo detection filters to dataframe using vectorized operations.

    dnm_config expected keys:
      - sample_dp_min, sample_gq_min, sample_vaf_min
      - parent_dp_min, parent_gq_min, parent_vaf_max
      - par_regions: dict with PAR regions keyed by assembly
      - y_excluded_regions: dict with Y regions to skip (ampliconic/PAR/XTR)

    This function will read ``sex`` from ``df`` when present; otherwise it will use
    the ``pedigree_df`` (which should contain ``sample_id`` and ``sex``).

    Args:
        df: DataFrame with melted samples
        dnm_config: DNM configuration dict
        verbose: Whether to print progress messages
        pedigree_df: Pedigree DataFrame
        skip_prefilters: If True, skips frequency/genomes_filters (assumes already applied)
    """
    if not dnm_config:
        return df

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
    y_excluded_regions = dnm_config.get("y_excluded_regions", {})

    original = df.shape[0]

    if "father_id" not in df.columns or "mother_id" not in df.columns:
        if "father" in df.columns and "mother" in df.columns:
            df = df.rename({"father": "father_id", "mother": "mother_id"})
        elif pedigree_df is not None:
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

    if "sex" not in df.columns:
        if pedigree_df is not None and "sex" in pedigree_df.columns:
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

    df = df.filter(
        ~pl.col("sample_gt").str.contains(r"\.") & pl.col("sample_gt").is_not_null()
    )

    if verbose:
        click.echo(
            f"DNM: {df.shape[0]} variants after removing samples with missing/partial GT",
            err=True,
        )

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

    if "#CHROM" in df.columns and df.schema["#CHROM"] == pl.Categorical:
        df = df.with_columns(pl.col("#CHROM").cast(pl.Utf8))

    df = df.with_columns(
        pl.col("#CHROM")
        .str.replace("^chr", "")
        .str.to_uppercase()
        .alias("_chrom_short")
    )

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

    y_excluded_mask = pl.lit(False)
    if y_excluded_regions:
        for assembly, regions in y_excluded_regions.items():
            for region_name, region in regions.items():
                region_chrom = _chrom_short(region.get("chrom", "Y")).upper()
                start = int(region.get("start"))
                end = int(region.get("end"))
                y_excluded_mask = y_excluded_mask | (
                    (pl.col("_chrom_short") == region_chrom)
                    & (pl.col("POS") >= start)
                    & (pl.col("POS") <= end)
                )
    df = df.with_columns(y_excluded_mask.alias("_y_excluded"))

    y_excluded_count = df.filter(pl.col("_y_excluded")).shape[0]
    df = df.filter(~pl.col("_y_excluded"))

    if verbose and y_excluded_count > 0:
        click.echo(
            f"DNM: Removed {y_excluded_count} variants in Y excluded regions "
            f"(ampliconic/PAR/XTR)",
            err=True,
        )

    df = df.with_columns(
        pl.col("sex").cast(pl.Utf8).str.to_uppercase().alias("_sex_norm")
    )

    is_hemizygous = (
        (pl.col("_chrom_short") == "X")
        & (pl.col("_sex_norm").is_in(["1", "M"]))
        & ~pl.col("_in_par")
    ) | (pl.col("_chrom_short") == "Y")

    is_homozygous = pl.col("sample_gt").is_in(["1/1", "2/2", "3/3"])

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

    if not skip_prefilters:
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

    father_qual_ok = (pl.col("father_dp").cast(pl.Float64, strict=False) >= p_dp) & (
        pl.col("father_gq").cast(pl.Float64, strict=False) >= p_gq
    )
    mother_qual_ok = (pl.col("mother_dp").cast(pl.Float64, strict=False) >= p_dp) & (
        pl.col("mother_gq").cast(pl.Float64, strict=False) >= p_gq
    )

    father_vaf_ok = pl.col("father_vaf").is_null() | (pl.col("father_vaf") < p_vaf)
    mother_vaf_ok = pl.col("mother_vaf").is_null() | (pl.col("mother_vaf") < p_vaf)

    father_gt_ok = (
        ~pl.col("father_gt").str.contains(r"\.") & pl.col("father_gt").is_not_null()
    )
    mother_gt_ok = (
        ~pl.col("mother_gt").str.contains(r"\.") & pl.col("mother_gt").is_not_null()
    )

    parent_filter = (
        pl.when(pl.col("_in_par"))
        .then(
            father_qual_ok
            & father_vaf_ok
            & father_gt_ok
            & mother_qual_ok
            & mother_vaf_ok
            & mother_gt_ok
        )
        .when(pl.col("_chrom_short") == "Y")
        .then(father_qual_ok & father_vaf_ok & father_gt_ok)
        .when((pl.col("_chrom_short") == "X") & (pl.col("_sex_norm").is_in(["1", "M"])))
        .then(father_qual_ok & mother_qual_ok & mother_vaf_ok & mother_gt_ok)
        .otherwise(
            father_qual_ok
            & father_vaf_ok
            & father_gt_ok
            & mother_qual_ok
            & mother_vaf_ok
            & mother_gt_ok
        )
    )

    result = df.filter(parent_filter)

    if verbose:
        click.echo(
            f"DNM: {result.shape[0]} variants after parent genotype filtering (parent DP>={p_dp}, GQ>={p_gq}, VAF<{p_vaf})",
            err=True,
        )

    temp_cols = ["_chrom_short", "_in_par", "_sex_norm", "_y_excluded"]
    result = result.drop([c for c in temp_cols if c in result.columns])
    if "father_id" in result.columns:
        result = result.rename({"father_id": "father", "mother_id": "mother"})

    if verbose:
        click.echo(f"De novo filter: {original} -> {result.shape[0]} rows", err=True)

    return result


def process_dnm_by_chromosome(
    input_file: Path,
    pedigree_df: pl.DataFrame,
    filter_config: dict,
    output_format: str,
    verbose: bool
) -> pl.DataFrame:
    """Process DNM filtering chromosome by chromosome to reduce memory usage.

    Processes each chromosome separately:
    1. Load one chromosome at a time from Parquet
    2. Apply frequency/quality prefilters (before melting)
    3. Melt samples
    4. Apply DNM filters
    5. Combine results from all chromosomes

    This reduces peak memory from (total_variants x samples) to
    (max_chr_variants x samples).

    Args:
        input_file: Path to Parquet file
        pedigree_df: Pedigree DataFrame with sample relationships
        filter_config: Filter configuration dict
        output_format: Output format (tsv, tsv.gz, parquet)
        verbose: Whether to print progress messages

    Returns:
        Combined DataFrame with DNM-filtered variants from all chromosomes
    """
    chromosomes = get_unique_chromosomes(input_file)

    if verbose:
        click.echo(
            f"DNM per-chromosome processing: {len(chromosomes)} chromosomes", err=True
        )

    results = []
    quality_cfg = filter_config.get("quality", {})
    dnm_section = filter_config.get("dnm", {})
    overlap = set(quality_cfg) & set(dnm_section) - {"enabled"}
    if overlap:
        warnings.warn(
            f"Keys {sorted(overlap)} appear in both 'quality' and 'dnm'; "
            f"'dnm' values will be used."
        )
    dnm_cfg = {}
    dnm_cfg.update(quality_cfg)
    dnm_cfg.update(dnm_section)

    expression = filter_config.get("expression")
    expr_filter = None
    if expression:
        schema_df = pl.scan_parquet(input_file).head(1).collect()
        try:
            expr_filter = parse_impact_filter_expression(expression, schema_df)
            if verbose:
                click.echo(
                    f"DNM: applying expression filter before melting: {expression}",
                    err=True,
                )
        except Exception as e:
            if verbose:
                click.echo(f"Warning: Could not apply early filter: {e}", err=True)

    for chrom in chromosomes:
        if verbose:
            click.echo(f"Processing chromosome {chrom}...", err=True)

        lazy_df = pl.scan_parquet(input_file).filter(
            pl.col("#CHROM") == chrom
        )

        if expr_filter is not None:
            lazy_df = lazy_df.filter(expr_filter)

        lazy_df = apply_dnm_prefilters(lazy_df, filter_config, verbose=False)

        if verbose:
            pre_count = lazy_df.select(pl.len()).collect().item()
            click.echo(f"  Chromosome {chrom}: {pre_count} variants after prefilter", err=True)

        df = lazy_df.collect()

        if df.shape[0] == 0:
            if verbose:
                click.echo(f"  Chromosome {chrom}: No variants after prefilter, skipping", err=True)
            continue

        formatted_df = format_bcftools_tsv_minimal(df, pedigree_df)

        if verbose:
            click.echo(
                f"  Chromosome {chrom}: {formatted_df.shape[0]} rows after melting", err=True
            )

        filtered_df = apply_de_novo_filter(
            formatted_df, dnm_cfg, verbose=False, pedigree_df=pedigree_df,
            skip_prefilters=True
        )

        if verbose:
            click.echo(
                f"  Chromosome {chrom}: {filtered_df.shape[0]} variants passed DNM filter", err=True
            )

        if filtered_df.shape[0] > 0:
            results.append(filtered_df)

    if not results:
        if verbose:
            click.echo("No variants passed DNM filters across all chromosomes", err=True)
        return pl.DataFrame()

    final_df = pl.concat(results)

    if verbose:
        click.echo(
            f"DNM filtering complete: {final_df.shape[0]} total variants", err=True
        )

    return final_df
