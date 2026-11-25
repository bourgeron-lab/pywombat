"""CLI for wombat tool."""

import gzip
import re
from pathlib import Path
from typing import Optional

import click
import polars as pl

from pywombat.acmg import annotate_acmg, load_ruleset
from pywombat.denovo import process_dnm_by_chromosome
from pywombat.expression import parse_impact_filter_expression
from pywombat.filters import apply_filters_lazy, load_filter_config
from pywombat.formatting import (
    format_bcftools_tsv,
    format_bcftools_tsv_lazy,
    format_bcftools_tsv_minimal,
)
from pywombat.pedigree import read_pedigree
from pywombat.prepare import prepare_parquet


def _normalize_annotate_mode(value) -> str:
    """Normalize the mnv.annotate config value to a standard string.

    Supports:
        True / "true" -> "vep_api" (backward compatible)
        False / None / "false" -> "false"
        "vep_api" -> "vep_api"
        "external" -> "external"

    Returns:
        One of: "false", "vep_api", "external"

    Raises:
        click.BadParameter: If value is not recognized.
    """
    if value is None or value is False:
        return "false"
    if value is True:
        return "vep_api"
    if isinstance(value, str):
        lower = value.strip().lower()
        if lower in ("false", "no", "off", "0"):
            return "false"
        if lower in ("true", "yes", "on", "1", "vep_api"):
            return "vep_api"
        if lower == "external":
            return "external"
    raise click.BadParameter(
        f"Invalid annotate mode: '{value}'. "
        "Must be false, vep_api, or external."
    )


def write_output(
    data,
    output_prefix: str,
    output_format: str,
    verbose: bool = False,
) -> Path:
    """Write a DataFrame or LazyFrame to the requested output format.

    This is the single exit point for filtered data — the post-filter hook
    where ACMG evidence columns will be attached in a future phase.

    Args:
        data: pl.DataFrame or pl.LazyFrame to write
        output_prefix: Output file path prefix (without extension)
        output_format: One of "tsv", "tsv.gz", "parquet"
        verbose: Whether to print progress messages

    Returns:
        The Path that was written to
    """
    output_path = Path(f"{output_prefix}.{output_format}")

    is_lazy = isinstance(data, pl.LazyFrame)

    if output_format == "tsv":
        if is_lazy:
            data.sink_csv(output_path, separator="\t")
        else:
            data.write_csv(output_path, separator="\t")
    elif output_format == "tsv.gz":
        if is_lazy:
            data = data.collect()
        csv_content = data.write_csv(separator="\t")
        with gzip.open(output_path, "wt") as f:
            f.write(csv_content)
    elif output_format == "parquet":
        if is_lazy:
            data.sink_parquet(output_path)
        else:
            data.write_parquet(output_path)

    if verbose:
        click.echo(f"Data written to {output_path}", err=True)

    return output_path


@click.group()
def cli():
    """
    Wombat: A tool for processing bcftools tabulated TSV files.

    \b
    Commands:
        filter             Process and filter variant data
        prepare            Convert TSV to optimized Parquet format
        vep                Annotate variant using Ensembl VEP REST API
        merge-annotations  Merge external VEP annotations into filtered TSV
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

        if not str(output).endswith(".parquet"):
            output = Path(f"{output}.parquet")

        prepare_parquet(input_file, output, chunk_size, verbose)

        if verbose:
            click.echo(f"Successfully created {output}", err=True)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


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
@click.option(
    "--bcf",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to tabix-indexed BCF file for MNV detection.",
)
@click.option(
    "--fasta",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to indexed FASTA reference for MNV reconstruction.",
)
@click.option(
    "--acmg-ruleset",
    type=str,
    default=None,
    help="ACMG evidence ruleset: built-in name (e.g. acmg-svi-2024.1) or path to YAML.",
)
def filter_cmd(
    input_file: Path,
    output: Optional[str],
    output_format: str,
    verbose: bool,
    pedigree: Optional[Path],
    filter_config: Optional[Path],
    debug: Optional[str],
    bcf: Optional[Path],
    fasta: Optional[Path],
    acmg_ruleset: Optional[str],
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

        is_parquet = str(input_file).endswith(".parquet")
        is_gzipped = str(input_file).endswith(".gz")

        if verbose:
            if is_parquet:
                click.echo("Detected Parquet input (pre-processed)", err=True)
            elif is_gzipped:
                click.echo("Detected gzipped TSV file", err=True)

        pedigree_df = None
        if pedigree:
            if verbose:
                click.echo(f"Reading pedigree file: {pedigree}", err=True)
            pedigree_df = read_pedigree(pedigree)

        filter_config_data = None
        if filter_config:
            if verbose:
                click.echo(f"Reading filter config: {filter_config}", err=True)
            filter_config_data = load_filter_config(filter_config)

        if debug:
            debug_variant(input_file, pedigree_df, filter_config_data, debug, verbose)
            return

        if output is None:
            input_stem = input_file.name
            for ext in [".tsv.gz", ".tsv", ".parquet"]:
                if input_stem.endswith(ext):
                    input_stem = input_stem[: -len(ext)]
                    break

            if filter_config:
                config_name = filter_config.stem
                output = f"{input_stem}.{config_name}"
            else:
                output = input_stem

        if verbose:
            click.echo("Processing with streaming mode...", err=True)

        if is_parquet:
            lazy_df = pl.scan_parquet(input_file)

            # Parquet + DNM: per-chromosome processing for memory efficiency
            if filter_config_data and filter_config_data.get("dnm", {}).get("enabled", False):
                if verbose:
                    click.echo("DNM mode: Using per-chromosome processing for memory efficiency", err=True)

                if pedigree_df is None:
                    click.echo("Error: DNM filtering requires a pedigree file (--pedigree option)", err=True)
                    raise click.Abort()

                formatted_df = process_dnm_by_chromosome(
                    input_file, pedigree_df, filter_config_data, output_format, verbose
                )

                if acmg_ruleset:
                    ruleset = load_ruleset(acmg_ruleset)
                    if verbose:
                        click.echo(f"Annotating ACMG evidence ({ruleset.get('name', acmg_ruleset)})...", err=True)
                    formatted_df = annotate_acmg(formatted_df, ruleset)

                write_output(formatted_df, output, output_format, verbose)
                return

            # Parquet + expression: apply BEFORE melting for 95%+ memory reduction
            if filter_config_data and "expression" in filter_config_data:
                expression = filter_config_data["expression"]
                if expression and verbose:
                    click.echo(
                        f"Applying expression filter before melting: {expression}",
                        err=True,
                    )

                schema_df = lazy_df.head(1).collect()
                try:
                    filter_expr = parse_impact_filter_expression(expression, schema_df)
                    lazy_df = lazy_df.filter(filter_expr)

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

            df = lazy_df.collect()
            formatted_df = format_bcftools_tsv_minimal(df, pedigree_df)
            lazy_df = formatted_df.lazy()

            # Remove expression from config so it's not applied again
            if filter_config_data and "expression" in filter_config_data:
                filter_config_data = filter_config_data.copy()
                del filter_config_data["expression"]
        else:
            # TSV input: need full processing
            string_columns = [
                "FID", "sample_id", "father_id", "mother_id",
                "FatherBarcode", "MotherBarcode", "sample",
            ]
            schema_overrides = {col: pl.Utf8 for col in string_columns}
            lazy_df = pl.scan_csv(
                input_file, separator="\t", schema_overrides=schema_overrides
            )

            lazy_df = format_bcftools_tsv_lazy(lazy_df, pedigree_df)

        if filter_config_data:
            lazy_df = apply_filters_lazy(
                lazy_df, filter_config_data, verbose, pedigree_df
            )

        # MNV detection
        mnv_config = filter_config_data.get("mnv", {}) if filter_config_data else {}
        if mnv_config.get("candidate", False):
            if bcf is None:
                raise click.UsageError(
                    "MNV detection requires --bcf option. "
                    "Provide path to tabix-indexed BCF file."
                )
            if fasta is None:
                raise click.UsageError(
                    "MNV detection requires --fasta option. "
                    "Provide path to indexed FASTA reference."
                )

            if verbose:
                click.echo("Detecting MNV candidates...", err=True)

            from pywombat.mnv import detect_mnv_probas

            df = lazy_df.collect()

            # Exclude MHC region if configured
            exclude_mhc = mnv_config.get("exclude_mhc", {})
            if exclude_mhc.get("enabled", False):
                mhc_chrom = exclude_mhc.get("chrom", "6")
                mhc_start = exclude_mhc.get("start", 28510120)
                mhc_end = exclude_mhc.get("end", 33480577)
                before_count = df.shape[0]

                chrom_col = pl.col("#CHROM").cast(pl.Utf8)
                is_mhc = (
                    (chrom_col == mhc_chrom)
                    | (chrom_col == f"chr{mhc_chrom}")
                ) & (pl.col("POS") >= mhc_start) & (pl.col("POS") <= mhc_end)
                df = df.filter(~is_mhc)

                if verbose:
                    removed = before_count - df.shape[0]
                    click.echo(
                        f"  Excluded MHC region (chr{mhc_chrom}:{mhc_start}-{mhc_end}): "
                        f"removed {removed} variants",
                        err=True,
                    )

            df = detect_mnv_probas(df, str(bcf), str(fasta), mnv_config)

            if verbose:
                n_candidates = df.filter(pl.col("mnv_proba").is_not_null()).shape[0]
                click.echo(f"  Found {n_candidates} MNV candidates", err=True)

            # VEP annotation of MNV candidates
            annotate_mode = _normalize_annotate_mode(
                mnv_config.get("annotate", False)
            )

            if annotate_mode == "vep_api":
                if verbose:
                    click.echo("Annotating MNV candidates via VEP API...", err=True)

                from pywombat.vep import annotate_mnv_variants

                df = annotate_mnv_variants(
                    df,
                    verbose=verbose,
                    batch_size=mnv_config.get("annotate_batch_size", 200),
                    timeout=mnv_config.get("annotate_timeout", 30),
                )
            elif annotate_mode == "external":
                if verbose:
                    click.echo(
                        "Creating MNV placeholders and VCF for external annotation...",
                        err=True,
                    )

                from pywombat.vep import create_mnv_placeholders_and_vcf

                df, vcf_path = create_mnv_placeholders_and_vcf(
                    df, output_prefix=output, verbose=verbose
                )

                if verbose and vcf_path.stat().st_size > 0:
                    click.echo(
                        f"  VCF for external annotation: {vcf_path}", err=True
                    )

            # Filter to MNV-only if configured
            if mnv_config.get("only", False):
                before_count = df.shape[0]
                df = df.filter(pl.col("mnv_proba").is_not_null())
                if verbose:
                    click.echo(
                        f"  MNV-only filter: {before_count} → {df.shape[0]} variants",
                        err=True,
                    )

            lazy_df = df.lazy()

        if acmg_ruleset:
            ruleset = load_ruleset(acmg_ruleset)
            if verbose:
                click.echo(f"Annotating ACMG evidence ({ruleset.get('name', acmg_ruleset)})...", err=True)
            df = lazy_df.collect() if isinstance(lazy_df, pl.LazyFrame) else lazy_df
            df = annotate_acmg(df, ruleset)
            lazy_df = df.lazy()

        write_output(lazy_df, output, output_format, verbose)

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

    string_columns = [
        "FID", "sample_id", "father_id", "mother_id",
        "FatherBarcode", "MotherBarcode", "sample",
    ]
    schema_overrides = {col: pl.Utf8 for col in string_columns}
    df = pl.read_csv(input_file, separator="\t", schema_overrides=schema_overrides)
    formatted_df = format_bcftools_tsv(df, pedigree_df)

    matching_rows = formatted_df.filter(
        (pl.col("#CHROM") == chrom) & (pl.col("POS") == pos)
    )

    if matching_rows.shape[0] == 0:
        click.echo(f"No rows found matching {chrom}:{pos}", err=True)
        return

    columns_to_show = ["#CHROM", "POS"]

    if "VEP_SYMBOL" in matching_rows.columns:
        columns_to_show.append("VEP_SYMBOL")

    if filter_config and "expression" in filter_config:
        expression = filter_config["expression"]
        column_pattern = r"\b([A-Za-z_][A-Za-z0-9_]*)\b\s*[=!<>]"
        found_columns = re.findall(column_pattern, expression)

        for col in found_columns:
            if col in matching_rows.columns and col not in columns_to_show:
                columns_to_show.append(col)

    display_df = matching_rows.select(
        [col for col in columns_to_show if col in matching_rows.columns]
    )

    for col in display_df.columns:
        if display_df[col].dtype in [pl.Float32, pl.Float64]:
            display_df = display_df.with_columns(
                pl.when(pl.col(col).is_null())
                .then(pl.lit("<null>"))
                .when(pl.col(col).is_nan())
                .then(pl.lit("<NaN>"))
                .otherwise(pl.col(col).cast(pl.Utf8))
                .alias(col)
            )
        else:
            display_df = display_df.with_columns(
                pl.when(pl.col(col).is_null())
                .then(pl.lit("<null>"))
                .otherwise(pl.col(col).cast(pl.Utf8))
                .alias(col)
            )

    click.echo(f"\nFound {matching_rows.shape[0]} row(s) matching {chrom}:{pos}:\n")
    click.echo(display_df.write_csv(separator="\t"))


@cli.command("vep")
@click.argument("variant", type=str)
@click.option(
    "--all",
    "show_all",
    is_flag=True,
    default=False,
    help="Show all transcripts instead of only canonical.",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def vep_cmd(variant: str, show_all: bool, verbose: bool):
    """
    Annotate a variant using Ensembl VEP REST API (GRCh38).

    \b
    Accepts variant in format: chrom:pos:ref:alt
    Examples:
        wombat vep chr1:151776102:GC:AA
        wombat vep 17:41245466:G:A --all
    """
    from pywombat.vep import (
        extract_annotations,
        format_annotations,
        parse_variant,
        query_vep,
    )

    try:
        parsed = parse_variant(variant)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()

    if verbose:
        click.echo(f"Querying Ensembl VEP API for {variant}...", err=True)

    try:
        responses = query_vep([parsed])
    except (ConnectionError, RuntimeError) as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()

    if not responses:
        click.echo("No results returned from VEP API.", err=True)
        raise click.Abort()

    if verbose:
        click.echo("Extracting annotations...", err=True)

    annotations = extract_annotations(responses[0], canonical_only=not show_all)
    output = format_annotations(annotations, variant)
    click.echo(output)


@cli.command("merge-annotations")
@click.option(
    "--tsv",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="TSV file with placeholder MNV rows (from annotate: external).",
)
@click.option(
    "--annotations",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="VEP tab-delimited output file with annotations.",
)
@click.option(
    "-o",
    "--output",
    type=str,
    required=True,
    help="Output file prefix.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    default="tsv",
    type=click.Choice(["tsv", "tsv.gz", "parquet"]),
    help="Output format (default: tsv).",
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def merge_annotations_cmd(
    tsv: Path,
    annotations: Path,
    output: str,
    output_format: str,
    verbose: bool,
):
    """
    Merge external VEP annotations back into a filtered TSV.

    \b
    After running 'wombat filter' with annotate: external, you get:
      1. A filtered TSV with placeholder MNV rows (VEP columns empty)
      2. A .to_annotate.vcf file with MNV variants

    Run VEP externally on the VCF, then use this command to merge
    the VEP annotations back into the TSV.

    \b
    Example workflow:
        # Step 1: Filter with external annotation
        wombat filter prepared.parquet -F config.yml --bcf c.bcf --fasta ref.fa -o out -v

        # Step 2: Run VEP externally
        vep --input_file out.to_annotate.vcf --tab --output_file out.annotated.tsv ...

        # Step 3: Merge annotations
        wombat merge-annotations --tsv out.tsv --annotations out.annotated.tsv -o merged -v
    """
    try:
        if verbose:
            click.echo(
                f"Merging VEP annotations from {annotations} into {tsv}...",
                err=True,
            )

        from pywombat.vep import merge_vep_annotations

        df = merge_vep_annotations(tsv, annotations, verbose=verbose)

        write_output(df, output, output_format, verbose)

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    cli()
