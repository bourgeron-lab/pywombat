"""Prepare command: TSV → Parquet conversion with chunked processing."""

import gzip
import io
import math
from pathlib import Path
from typing import Optional

import click
import polars as pl

from pywombat.denovo import apply_de_novo_filter
from pywombat.expression import parse_impact_filter_expression
from pywombat.filters import apply_quality_filters
from pywombat.formatting import format_bcftools_tsv_minimal, format_expand_annotations
from pywombat.pedigree import add_parent_genotypes_from_lookup


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

    if verbose:
        click.echo("Pass 2: Converting to Parquet...", err=True)

    dtype_overrides = {
        "#CHROM": pl.Categorical,
        "POS": pl.UInt32,
        "FILTER": pl.Categorical,
    }

    import tempfile
    import shutil

    temp_dir = Path(tempfile.mkdtemp(prefix="wombat_prepare_"))
    part_files = []

    try:
        with opener(input_file, "rt") as f:
            header_line = f.readline().strip()

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

            if chunk_lines:
                df_chunk = _process_chunk(
                    header_line, chunk_lines, all_fields, all_flags, dtype_overrides
                )
                part_file = temp_dir / f"part_{len(part_files):06d}.parquet"
                df_chunk.write_parquet(part_file)
                part_files.append(part_file)
                pbar.update(len(chunk_lines))

            pbar.close()

        if verbose:
            click.echo(f"Combining {len(part_files)} parts into final output...", err=True)

        if part_files:
            combined = pl.scan_parquet(part_files).collect()
            combined.write_parquet(output)

        if verbose:
            click.echo(f"Wrote {len(part_files)} chunks to {output}", err=True)

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def _process_chunk(
    header: str,
    lines: list,
    fields: set,
    flags: set,
    dtype_overrides: dict,
) -> pl.DataFrame:
    """Process a chunk of lines into a DataFrame with expanded INFO fields."""
    content = header + "\n" + "".join(lines)
    df = pl.read_csv(
        io.StringIO(content),
        separator="\t",
        infer_schema_length=10000,
    )

    if "(null)" in df.columns:
        for field in sorted(fields):
            df = df.with_columns(
                pl.col("(null)").str.extract(f"{field}=([^;]+)").alias(field)
            )

        for flag in sorted(flags):
            df = df.with_columns(
                pl.col("(null)").str.contains(f"(^|;){flag}(;|$)").alias(flag)
            )

        df = df.drop("(null)")

    if "CSQ" in df.columns:
        df = df.drop("CSQ")

    for col, dtype in dtype_overrides.items():
        if col in df.columns:
            try:
                df = df.with_columns(pl.col(col).cast(dtype))
            except Exception:
                pass

    return df


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
    of data lines processed (excluding header). If a ``progress_bar`` is provided it will
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
            processed += len(chunk)
            if progress_bar is not None:
                progress_bar.update(len(chunk))
            elif verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Building lookup: processed ~{processed} lines...", err=True
                )
            continue

        try:
            formatted = format_bcftools_tsv_minimal(df_chunk, pedigree_df=None)
        except Exception:
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
        part = formatted.select(sel).filter(pl.col("sample_gt") != "0/0").unique()
        parts.append(part)

        processed += len(chunk)
        if progress_bar is not None:
            progress_bar.update(len(chunk))
        elif verbose and (chunk_idx % 10 == 0):
            click.echo(f"Building lookup: processed ~{processed} lines...", err=True)

    if parts:
        lookup = pl.concat(parts).unique()
    else:
        lookup = pl.DataFrame([])

    return lookup, processed


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
    try:
        from tqdm.auto import tqdm
    except Exception:
        tqdm = None

    if verbose:
        click.echo("Pass 1: building parent genotype lookup (single pass)...", err=True)

    pbar_lookup = None
    if tqdm is not None:
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

    if output_format == "tsv":
        out_path = Path(f"{output_prefix}.tsv")
    elif output_format == "tsv.gz":
        out_path = Path(f"{output_prefix}.tsv.gz")
    else:
        out_path = Path(f"{output_prefix}")

    first_write = True

    iterator = _line_iterator(input_path, chunk_size=chunk_size)
    chunk_idx = 0

    progress_bar = None
    processed_lines = 0
    if tqdm is not None:
        progress_bar = tqdm(total=total_lines, desc="Processing variants")

    for header, chunk in iterator:
        chunk_idx += 1
        chunk_count = len(chunk)

        if progress_bar is not None:
            progress_bar.set_postfix_str(f"chunk {chunk_idx}")

        content = header + "".join(chunk)
        df_chunk = pl.read_csv(io.StringIO(content), separator="\t")

        try:
            melted = format_bcftools_tsv_minimal(df_chunk, pedigree_df=None)
        except Exception:
            processed_lines += chunk_count
            if progress_bar is not None:
                progress_bar.update(chunk_count)
            elif verbose and (chunk_idx % 10 == 0):
                click.echo(
                    f"Processed {processed_lines}/{total_lines} lines...", err=True
                )
            continue

        if filter_config and filter_config.get("dnm", {}).get("enabled", False):
            melted = melted.filter(
                pl.col("sample_gt").str.contains("1")
                | pl.col("sample_gt").str.contains("2")
            )
            if melted.shape[0] == 0:
                processed_lines += chunk_count
                if progress_bar is not None:
                    progress_bar.update(chunk_count)
                continue

        if pedigree_df is not None and "sample" in melted.columns:
            melted = melted.join(
                pedigree_df, left_on="sample", right_on="sample_id", how="left"
            )

        melted = add_parent_genotypes_from_lookup(melted, parent_lookup)

        if filter_config and filter_config.get("dnm", {}).get("enabled", False):
            cfg = {}
            cfg.update(filter_config.get("quality", {}))
            cfg.update(filter_config.get("dnm", {}))
            filtered = apply_de_novo_filter(
                melted, cfg, verbose=False, pedigree_df=pedigree_df
            )
        else:
            filtered = melted
            quality_cfg = filter_config.get("quality", {}) if filter_config else {}
            if quality_cfg:
                filtered = apply_quality_filters(filtered, quality_cfg, verbose=False)

            expr = filter_config.get("expression") if filter_config else None
            if expr:
                try:
                    expr_parsed = parse_impact_filter_expression(expr, filtered)
                    filtered = filtered.filter(expr_parsed)
                except Exception:
                    pass

        if filtered.shape[0] > 0 and "(null)" in filtered.columns:
            if progress_bar is not None:
                progress_bar.set_postfix_str(f"expanding annotations chunk {chunk_idx}")
            filtered = format_expand_annotations(filtered)

        if progress_bar is not None:
            progress_bar.set_postfix_str(f"writing chunk {chunk_idx}")

        if filtered.shape[0] > 0:
            if output_format in ("tsv", "tsv.gz"):
                csv_text = filtered.write_csv(separator="\t")
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
                    tail = "\n".join(csv_text.splitlines()[1:])
                    if output_format == "tsv.gz":
                        with gzip.open(out_path, "at") as f:
                            f.write("\n" + tail)
                    else:
                        with open(out_path, "at") as f:
                            f.write("\n" + tail)
            else:
                part_path = out_path.with_suffix(f".part{chunk_idx}.parquet")
                filtered.write_parquet(part_path)

        if progress_bar is not None:
            progress_bar.update(chunk_count)
            progress_bar.set_postfix_str("")
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
