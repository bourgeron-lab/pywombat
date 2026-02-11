"""VEP (Variant Effect Predictor) annotation via Ensembl REST API."""

import json
import re
import sys
import urllib.error
import urllib.request

import click

VEP_API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"

DEFAULT_PARAMS = {
    "canonical": 1,
    "LoF": 1,
    "hgvs": 1,
    "protein": 1,
    "mane": 1,
    "numbers": 1,
}

_VALID_BASES = re.compile(r"^[ACGTacgt]+$")

# Maps VEP annotation dict keys (from extract_annotations) to DataFrame VEP_* column names.
# Only columns that exist in the DataFrame will be populated; others are ignored.
VEP_COLUMN_MAPPING = {
    "gene_symbol": "VEP_SYMBOL",
    "gene_id": "VEP_Gene",
    "transcript_id": "VEP_Feature",
    "consequence": "VEP_Consequence",
    "impact": "VEP_IMPACT",
    "biotype": "VEP_BIOTYPE",
    "hgvsc": "VEP_HGVSc",
    "hgvsp": "VEP_HGVSp",
    "protein_position": "VEP_Protein_position",
    "amino_acids": "VEP_Amino_acids",
    "codons": "VEP_Codons",
    "exon": "VEP_EXON",
    "intron": "VEP_INTRON",
    "mane_select": "VEP_MANE_SELECT",
    "mane_plus_clinical": "VEP_MANE_PLUS_CLINICAL",
    "lof": "VEP_LoF",
    "lof_filter": "VEP_LoF_filter",
    "lof_flags": "VEP_LoF_flags",
    "lof_info": "VEP_LoF_info",
    "canonical": "VEP_CANONICAL",
    "sift_prediction": "VEP_SIFT",
    "polyphen_prediction": "VEP_PolyPhen",
}

# Sample columns to copy from the original row to the new MNV row.
_SAMPLE_COLUMNS = [
    "sample", "sample_gt", "sample_dp", "sample_gq", "sample_ad", "sample_vaf",
]


def parse_variant(variant_str: str) -> dict:
    """Parse variant string into components.

    Args:
        variant_str: Variant in format 'chr1:151776102:GC:AA' or '1:151776102:GC:AA'

    Returns:
        Dict with keys: chrom, pos, ref, alt, chrom_raw

    Raises:
        ValueError: If format is invalid.
    """
    parts = variant_str.strip().split(":")
    if len(parts) != 4:
        raise ValueError(
            f"Invalid variant format: '{variant_str}'. "
            "Expected chrom:pos:ref:alt (e.g. chr1:151776102:GC:AA)"
        )

    chrom_raw, pos_str, ref, alt = parts

    # Strip 'chr' prefix for VEP API
    chrom = chrom_raw.removeprefix("chr").removeprefix("Chr")
    if not chrom:
        raise ValueError(f"Invalid chromosome: '{chrom_raw}'")

    try:
        pos = int(pos_str)
    except ValueError:
        raise ValueError(f"Invalid position: '{pos_str}' (must be integer)")

    if pos < 1:
        raise ValueError(f"Invalid position: {pos} (must be positive)")

    ref = ref.upper()
    alt = alt.upper()

    if not _VALID_BASES.match(ref):
        raise ValueError(f"Invalid REF allele: '{ref}' (must be A/C/G/T)")
    if not _VALID_BASES.match(alt):
        raise ValueError(f"Invalid ALT allele: '{alt}' (must be A/C/G/T)")

    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "chrom_raw": chrom_raw,
    }


def _format_variant_for_vep(variant: dict) -> str:
    """Convert parsed variant to VEP API input format.

    Returns VCF-like string: '1 151776102 . GC AA . . .'
    """
    return f"{variant['chrom']} {variant['pos']} . {variant['ref']} {variant['alt']} . . ."


def query_vep(
    variants: list[dict],
    params: dict | None = None,
    timeout: int = 30,
) -> list[dict]:
    """Query Ensembl VEP REST API for variant annotations.

    Args:
        variants: List of parsed variant dicts (from parse_variant).
        params: Additional API parameters. Merged with DEFAULT_PARAMS.
        timeout: Request timeout in seconds.

    Returns:
        List of VEP response dicts (one per input variant).

    Raises:
        ConnectionError: If the API is unreachable.
        RuntimeError: If the API returns an error response.
    """
    merged_params = {**DEFAULT_PARAMS, **(params or {})}
    variant_strings = [_format_variant_for_vep(v) for v in variants]

    body = json.dumps({"variants": variant_strings, **merged_params}).encode("utf-8")

    req = urllib.request.Request(
        VEP_API_URL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        error_body = ""
        try:
            error_body = e.read().decode("utf-8")
        except Exception:
            pass
        raise RuntimeError(
            f"VEP API error {e.code}: {error_body or e.reason}"
        ) from e
    except urllib.error.URLError as e:
        raise ConnectionError(
            f"Cannot reach VEP API: {e.reason}"
        ) from e


def extract_annotations(
    vep_response: dict,
    canonical_only: bool = True,
) -> list[dict]:
    """Extract transcript annotations from a VEP response.

    Args:
        vep_response: Single variant response dict from VEP API.
        canonical_only: If True, return only the canonical transcript(s).

    Returns:
        List of annotation dicts with standardized keys.
    """
    transcripts = vep_response.get("transcript_consequences", [])

    if canonical_only:
        transcripts = [t for t in transcripts if t.get("canonical") == 1]

    annotations = []
    for tc in transcripts:
        consequence_terms = tc.get("consequence_terms", [])
        annotations.append(
            {
                "gene_symbol": tc.get("gene_symbol"),
                "gene_id": tc.get("gene_id"),
                "transcript_id": tc.get("transcript_id"),
                "consequence": "&".join(consequence_terms),
                "impact": tc.get("impact"),
                "biotype": tc.get("biotype"),
                "canonical": tc.get("canonical") == 1,
                "hgvsc": tc.get("hgvsc"),
                "hgvsp": tc.get("hgvsp"),
                "protein_position": tc.get("protein_position"),
                "amino_acids": tc.get("amino_acids"),
                "codons": tc.get("codons"),
                "exon": tc.get("exon"),
                "intron": tc.get("intron"),
                "mane_select": tc.get("mane_select"),
                "mane_plus_clinical": tc.get("mane_plus_clinical"),
                "lof": tc.get("lof"),
                "lof_filter": tc.get("lof_filter"),
                "lof_flags": tc.get("lof_flags"),
                "lof_info": tc.get("lof_info"),
                "sift_prediction": tc.get("sift_prediction"),
                "sift_score": tc.get("sift_score"),
                "polyphen_prediction": tc.get("polyphen_prediction"),
                "polyphen_score": tc.get("polyphen_score"),
            }
        )

    return annotations


def _transcript_base_id(transcript_id: str) -> str:
    """Extract base transcript ID without version (ENST00000357654.8 -> ENST00000357654)."""
    return transcript_id.split(".")[0] if transcript_id else ""


def _find_transcript_annotation(
    vep_response: dict,
    target_transcript_id: str,
) -> dict | None:
    """Find annotation for a specific transcript in a VEP response.

    Searches transcript_consequences for one matching the target ENST ID.
    Handles version mismatch (e.g. ENST00000357654.8 vs ENST00000357654)
    by comparing base IDs.

    Args:
        vep_response: Single variant VEP response dict.
        target_transcript_id: Ensembl transcript ID to match.

    Returns:
        Annotation dict from extract_annotations() or None if not found.
    """
    all_annotations = extract_annotations(vep_response, canonical_only=False)
    target_base = _transcript_base_id(target_transcript_id)

    for ann in all_annotations:
        ann_base = _transcript_base_id(ann.get("transcript_id", ""))
        if ann_base == target_base:
            return ann

    return None


def format_annotations(annotations: list[dict], variant_str: str) -> str:
    """Format annotations for CLI display.

    Args:
        annotations: List of annotation dicts from extract_annotations().
        variant_str: Original variant input string for header.

    Returns:
        Formatted string for terminal output.
    """
    lines = [f"Variant: {variant_str}", ""]

    if not annotations:
        lines.append("No transcript annotations found.")
        return "\n".join(lines)

    for ann in annotations:
        # Transcript header
        tid = ann.get("transcript_id", "unknown")
        label = f"--- {tid}"
        if ann.get("canonical"):
            label += " (canonical)"
        label += " ---"
        lines.append(label)

        # Key-value fields to display
        fields = [
            ("Gene", ann.get("gene_symbol")),
            ("Consequence", ann.get("consequence")),
            ("Impact", ann.get("impact")),
            ("HGVS.c", ann.get("hgvsc")),
            ("HGVS.p", ann.get("hgvsp")),
            ("Protein pos", ann.get("protein_position")),
            ("Amino acids", ann.get("amino_acids")),
            ("Codons", ann.get("codons")),
            ("Exon", ann.get("exon")),
            ("Intron", ann.get("intron")),
            ("Biotype", ann.get("biotype")),
            ("MANE Select", ann.get("mane_select")),
            ("MANE+ Clinical", ann.get("mane_plus_clinical")),
        ]

        # SIFT/PolyPhen
        sift_pred = ann.get("sift_prediction")
        sift_score = ann.get("sift_score")
        if sift_pred is not None:
            sift_val = sift_pred
            if sift_score is not None:
                sift_val += f" ({sift_score})"
            fields.append(("SIFT", sift_val))

        polyphen_pred = ann.get("polyphen_prediction")
        polyphen_score = ann.get("polyphen_score")
        if polyphen_pred is not None:
            pp_val = polyphen_pred
            if polyphen_score is not None:
                pp_val += f" ({polyphen_score})"
            fields.append(("PolyPhen", pp_val))

        # LOFTEE fields
        lof = ann.get("lof")
        if lof is not None:
            fields.append(("LoF", lof))
            fields.append(("LoF filter", ann.get("lof_filter") or "-"))
            fields.append(("LoF flags", ann.get("lof_flags") or "-"))
            lof_info = ann.get("lof_info")
            if lof_info:
                fields.append(("LoF info", lof_info))

        # Render fields, skip None values
        max_label_len = max(len(label) for label, val in fields if val is not None)
        for label, val in fields:
            if val is not None:
                lines.append(f"  {label:<{max_label_len}}  {val}")

        lines.append("")

    return "\n".join(lines)


def annotate_mnv_variants(
    df: "pl.DataFrame",
    verbose: bool = False,
    batch_size: int = 200,
    timeout: int = 30,
) -> "pl.DataFrame":
    """Annotate MNV candidates via VEP and append new rows to the DataFrame.

    For each row where mnv_proba is not null and mnv_variant is non-empty,
    queries the Ensembl VEP REST API for the MNV variant, finds the transcript
    matching the original row's VEP_Feature, and builds a new row.

    New rows have:
      - #CHROM/POS/REF/ALT from the mnv_variant (the reconstructed MNV)
      - sample columns copied from the original row
      - parent columns set to null
      - mnv_proba = 1.0 - original value
      - VEP_* columns populated from the VEP response
      - all other columns null

    Args:
        df: DataFrame with MNV columns and VEP_Feature column.
        verbose: Print progress to stderr.
        batch_size: Max variants per VEP API batch (max 200).
        timeout: HTTP timeout per request in seconds.

    Returns:
        DataFrame with original rows plus new MNV-annotated rows appended.
    """
    import polars as pl

    all_columns = df.columns

    # Check prerequisites
    if "VEP_Feature" not in all_columns:
        if verbose:
            click.echo(
                "Warning: VEP_Feature column not found, skipping MNV annotation.",
                err=True,
            )
        return df

    if "mnv_proba" not in all_columns or "mnv_variant" not in all_columns:
        return df

    # Step 1: Find rows eligible for annotation
    mnv_rows = df.filter(
        pl.col("mnv_proba").is_not_null()
        & pl.col("mnv_variant").is_not_null()
        & (pl.col("mnv_variant") != "")
        & (pl.col("mnv_variant") != ":::")
        & pl.col("VEP_Feature").is_not_null()
        & (pl.col("VEP_Feature") != "")
        & (pl.col("VEP_Feature") != ".")
    )

    if mnv_rows.is_empty():
        if verbose:
            click.echo("  No MNV candidates eligible for VEP annotation.", err=True)
        return df

    # Step 2: Deduplicate variants for efficient querying
    unique_variant_strs = (
        mnv_rows.select("mnv_variant").unique().to_series().to_list()
    )

    # Step 3: Parse variants
    parsed_variants = {}  # variant_str -> parsed dict
    for v_str in unique_variant_strs:
        try:
            parsed_variants[v_str] = parse_variant(v_str)
        except ValueError as e:
            if verbose:
                click.echo(
                    f"  Warning: Skipping unparseable MNV variant '{v_str}': {e}",
                    err=True,
                )

    if not parsed_variants:
        return df

    # Step 4: Query VEP in batches
    items = list(parsed_variants.items())  # [(variant_str, parsed_dict), ...]
    vep_results = {}  # variant_str -> VEP response dict

    for batch_start in range(0, len(items), batch_size):
        batch = items[batch_start : batch_start + batch_size]
        batch_variants = [pv for _, pv in batch]

        if verbose:
            click.echo(
                f"  Querying VEP API for {len(batch)} MNV variant(s)...", err=True
            )

        try:
            responses = query_vep(batch_variants, timeout=timeout)
            # Map responses back by matching the input field
            for resp in responses:
                input_str = resp.get("input", "")
                # Find which variant string matches this response
                for v_str, pv in batch:
                    vep_input = _format_variant_for_vep(pv)
                    if vep_input == input_str:
                        vep_results[v_str] = resp
                        break
        except (ConnectionError, RuntimeError) as e:
            if verbose:
                click.echo(f"  Warning: VEP API batch failed: {e}", err=True)

    if not vep_results:
        if verbose:
            click.echo(
                "  No VEP results obtained, skipping MNV annotation.", err=True
            )
        return df

    # Step 5: Determine which columns to map
    vep_columns_in_df = {
        ann_key: col_name
        for ann_key, col_name in VEP_COLUMN_MAPPING.items()
        if col_name in all_columns
    }

    parent_cols = [
        c for c in all_columns if c.startswith("father_") or c.startswith("mother_")
    ]

    sample_cols_present = [c for c in _SAMPLE_COLUMNS if c in all_columns]

    # Step 6: Build new rows
    new_rows = []
    skipped = 0

    for row in mnv_rows.iter_rows(named=True):
        v_str = row["mnv_variant"]
        vep_resp = vep_results.get(v_str)
        if vep_resp is None:
            skipped += 1
            continue

        target_transcript = row["VEP_Feature"]
        annotation = _find_transcript_annotation(vep_resp, target_transcript)
        if annotation is None:
            skipped += 1
            if verbose:
                click.echo(
                    f"  Warning: Transcript {target_transcript} not found "
                    f"in VEP response for MNV {v_str}",
                    err=True,
                )
            continue

        parsed = parsed_variants.get(v_str)
        if parsed is None:
            skipped += 1
            continue

        # Initialize all columns to None
        new_row = {col: None for col in all_columns}

        # Set genomic coordinates from MNV variant
        new_row["#CHROM"] = parsed["chrom_raw"]
        new_row["POS"] = parsed["pos"]
        new_row["REF"] = parsed["ref"]
        new_row["ALT"] = parsed["alt"]

        # Copy sample columns from original row
        for sc in sample_cols_present:
            new_row[sc] = row[sc]

        # Parent columns stay None (already set)

        # Set MNV-specific columns
        new_row["mnv_proba"] = 1.0 - row["mnv_proba"]
        new_row["mnv_variant"] = None
        new_row["mnv_inframe"] = None

        # Populate VEP columns from annotation
        for ann_key, col_name in vep_columns_in_df.items():
            value = annotation.get(ann_key)
            if ann_key == "canonical":
                new_row[col_name] = "YES" if value else None
            elif ann_key == "sift_prediction" and value is not None:
                score = annotation.get("sift_score")
                new_row[col_name] = (
                    f"{value}({score})" if score is not None else value
                )
            elif ann_key == "polyphen_prediction" and value is not None:
                score = annotation.get("polyphen_score")
                new_row[col_name] = (
                    f"{value}({score})" if score is not None else value
                )
            elif value is not None:
                new_row[col_name] = value

        new_rows.append(new_row)

    if not new_rows:
        if verbose:
            click.echo(
                f"  No MNV variants could be annotated (skipped: {skipped})",
                err=True,
            )
        return df

    # Step 7: Concatenate new rows with the original DataFrame
    new_df = pl.DataFrame(new_rows, schema=df.schema)
    result = pl.concat([df, new_df], how="vertical_relaxed")

    if verbose:
        click.echo(f"  Added {len(new_rows)} MNV-annotated row(s)", err=True)
        if skipped > 0:
            click.echo(
                f"  Skipped {skipped} MNV variant(s) (no VEP match)", err=True
            )

    return result
