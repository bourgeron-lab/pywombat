"""MNV (Multi-Nucleotide Variant) detection module.

This module provides functionality to detect MNV candidates by:
1. Finding nearby variants in the same sample (SNVs within 2bp, indels within configurable window)
2. Calculating phasing probability based on allelic depth (AD) values
3. Reconstructing MNV notation for SNV pairs
4. Determining inframe status for indel clusters
"""

import numpy as np
import polars as pl
from scipy.stats import norm

try:
    import cyvcf2
except ImportError:
    cyvcf2 = None

try:
    import pysam
except ImportError:
    pysam = None


def _check_dependencies():
    """Check if required dependencies are installed."""
    if cyvcf2 is None:
        raise ImportError(
            "MNV detection requires cyvcf2. Install with: pip install pywombat[mnv]"
        )
    if pysam is None:
        raise ImportError(
            "MNV detection requires pysam. Install with: pip install pywombat[mnv]"
        )


def calculate_phasing_probability(
    var1_ad: tuple[int, int],
    var2_ad: tuple[int, int]
) -> float:
    """Calculate probability that two variants are in-phase (cis).

    Uses a Gaussian model to compare allelic depth patterns:
    - In-phase: AD values should be similar (var1_ref≈var2_ref, var1_alt≈var2_alt)
    - Out-of-phase: AD values swapped (var1_ref≈var2_alt, var1_alt≈var2_ref)

    Based on Snakefile-MNVs phasing logic (lines 52-67).

    Args:
        var1_ad: (ref_depth, alt_depth) for variant 1
        var2_ad: (ref_depth, alt_depth) for variant 2

    Returns:
        Probability between 0.0 and 1.0 (1.0 for uncertain/low coverage)
    """
    a_0, a_1 = var1_ad
    b_0, b_1 = var2_ad

    # Handle missing or very low coverage
    if np.isnan(b_0) or b_0 <= 1:
        return 1.0

    scale = 1  # Gaussian scale parameter

    # Calculate differences for in-phase (IP) and out-of-phase (OP) models
    d_0 = abs(a_0 - b_0)  # Ref depth difference (in-phase)
    d_1 = abs(a_1 - b_1)  # Alt depth difference (in-phase)
    c_0 = abs(a_0 - b_1)  # Ref vs alt difference (out-of-phase)
    c_1 = abs(a_1 - b_0)  # Alt vs ref difference (out-of-phase)

    # Use log probabilities for numerical stability with large differences
    if d_0 < 25 and d_1 < 25 and c_0 < 25 and c_1 < 25:
        # Small differences - can use regular probability
        p_ip = norm.sf(d_0, scale=scale) * norm.sf(d_1, scale=scale)
        p_op = norm.sf(c_0, scale=scale) * norm.sf(c_1, scale=scale)
        return p_ip / (p_ip + p_op)
    else:
        # Large differences - use log probability to avoid underflow
        p_ip = norm.logsf(d_0, scale=scale) + norm.logsf(d_1, scale=scale)
        p_op = norm.logsf(c_0, scale=scale) + norm.logsf(c_1, scale=scale)
        return np.exp(p_ip - max(p_ip, p_op))


def calculate_indel_inframe(indels: list[dict]) -> bool:
    """Determine if indel cluster maintains reading frame.

    Args:
        indels: List of indel dicts with 'ref' and 'alt' keys

    Returns:
        True if combined indels are inframe (sum is multiple of 3)
    """
    total_net_length = 0
    for indel in indels:
        net_length = len(indel['alt']) - len(indel['ref'])
        total_net_length += net_length

    return total_net_length % 3 == 0


def get_window_variants(
    vcf: cyvcf2.VCF,
    chrom: str,
    start: int,
    end: int,
    sample: str
) -> list[dict]:
    """Query BCF for variants in window carried by sample.

    Args:
        vcf: cyvcf2.VCF object (already opened)
        chrom: Chromosome name
        start: Window start position (0-based)
        end: Window end position (0-based, exclusive)
        sample: Sample name to check genotype

    Returns:
        List of variant dicts with pos, ref, alt, gt, ad fields
    """
    # Get sample index
    try:
        sample_idx = vcf.samples.index(sample)
    except ValueError:
        # Sample not found in VCF
        return []

    variants = []
    region = f"{chrom}:{start}-{end}"

    try:
        for variant in vcf(region):
            # Skip non-ACTG variants (e.g., *)
            if variant.REF == '*' or any(alt == '*' for alt in variant.ALT):
                continue

            # Check if sample carries the variant (GT contains "1")
            gt = variant.genotypes[sample_idx]
            gt_string = '/'.join(str(g) for g in gt[:2])  # GT is [allele1, allele2, phased]

            if '1' not in gt_string:
                continue

            # Extract AD field
            ad_format_idx = variant.FORMAT.index('AD') if 'AD' in variant.FORMAT else None
            if ad_format_idx is not None:
                ad_values = variant.format('AD')[sample_idx]
                ad = tuple(int(ad_values[i]) for i in range(min(2, len(ad_values))))
            else:
                ad = (0, 0)  # Missing AD

            variants.append({
                'pos': variant.POS,
                'ref': variant.REF,
                'alt': variant.ALT[0] if variant.ALT else '',
                'gt': gt_string,
                'ad': ad
            })

    except Exception:
        # Region query failed (e.g., chromosome not in VCF)
        return []

    return variants


def reconstruct_snv_variant(
    chrom: str,
    var1: dict,
    var2: dict,
    fasta: pysam.FastaFile
) -> str:
    """Reconstruct MNV notation by filling gaps with reference.

    Based on Snakefile-MNVs rebuild_variant logic (lines 70-82).

    Args:
        chrom: Chromosome name
        var1: Variant dict {'pos': int, 'ref': str, 'alt': str}
        var2: Variant dict {'pos': int, 'ref': str, 'alt': str}
        fasta: pysam.FastaFile object

    Returns:
        Variant string "chr:pos:ref:alt" or ":::" if positions identical
    """
    p1 = var1['pos']
    p2 = var2['pos']
    r1 = var1['ref']
    a1 = var1['alt']
    r2 = var2['ref']
    a2 = var2['alt']

    # Positions must be different
    if p1 == p2:
        return ":::"

    # Order variants by position
    if p1 < p2:
        # var1 is leftmost
        # Fetch gap sequence (0-based coordinates for pysam)
        # Gap is from end of var1 to start of var2
        gap = fasta.fetch(chrom, p1 + len(r1) - 1, p2 - 1)
        mnv_id = f"{chrom}:{p1}:{r1}{gap}{r2}:{a1}{gap}{a2}"
    else:
        # var2 is leftmost
        gap = fasta.fetch(chrom, p2 + len(r2) - 1, p1 - 1)
        mnv_id = f"{chrom}:{p2}:{r2}{gap}{r1}:{a2}{gap}{a1}"

    return mnv_id


def _is_snv(ref: str, alt: str) -> bool:
    """Check if variant is a SNV (single nucleotide variant)."""
    valid_bases = {'A', 'C', 'G', 'T'}
    return (len(ref) == 1 and len(alt) == 1 and
            ref in valid_bases and alt in valid_bases)


def _is_indel(ref: str, alt: str) -> bool:
    """Check if variant is an indel."""
    return not _is_snv(ref, alt)


def process_variant_for_mnv(
    row: dict,
    vcf: cyvcf2.VCF,
    fasta: pysam.FastaFile,
    indel_window: int
) -> dict:
    """Process single variant-sample row to detect MNV.

    Args:
        row: Dict with keys: #CHROM, POS, REF, ALT, sample, sample_ad
        vcf: cyvcf2.VCF object
        fasta: pysam.FastaFile object
        indel_window: Window size for indel clustering

    Returns:
        Dict with mnv_proba, mnv_inframe, mnv_variant
    """
    chrom = row['#CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    sample = row['sample']

    # Parse sample AD
    try:
        # AD might be stored as string like "10,5" or as integer (just alt count)
        if isinstance(row.get('sample_ad'), str) and ',' in row['sample_ad']:
            ad_parts = row['sample_ad'].split(',')
            var_ad = (int(ad_parts[0]), int(ad_parts[1]))
        else:
            # If only alt count available, need ref count from DP
            ad_alt = int(row.get('sample_ad', 0))
            dp = int(row.get('sample_dp', 0))
            ad_ref = dp - ad_alt
            var_ad = (ad_ref, ad_alt)
    except (ValueError, TypeError):
        var_ad = (0, 0)

    # Determine variant type and search window
    if _is_snv(ref, alt):
        # SNV: search within 2bp window
        window_start = pos - 2
        window_end = pos + 2
        is_snv = True
    else:
        # Indel: search within configurable window
        indel_len = abs(len(ref) - len(alt))
        window_start = pos - indel_window
        window_end = pos + indel_len + indel_window
        is_snv = False

    # Query window for other variants
    window_vars = get_window_variants(vcf, chrom, window_start, window_end, sample)

    # Filter out the current variant and same type
    if is_snv:
        other_vars = [v for v in window_vars
                      if v['pos'] != pos and _is_snv(v['ref'], v['alt'])]
    else:
        other_vars = [v for v in window_vars
                      if v['pos'] != pos and _is_indel(v['ref'], v['alt'])]

    # No MNV candidate found
    if not other_vars:
        return {
            'mnv_proba': None,
            'mnv_inframe': None,
            'mnv_variant': ''
        }

    # Calculate phasing probabilities
    probabilities = [
        calculate_phasing_probability(var_ad, other['ad'])
        for other in other_vars
    ]

    # Take the best match (highest probability)
    best_idx = np.argmax(probabilities)
    best_prob = probabilities[best_idx]
    best_other = other_vars[best_idx]

    if is_snv:
        # Reconstruct SNV variant
        try:
            mnv_var = reconstruct_snv_variant(
                chrom,
                {'pos': pos, 'ref': ref, 'alt': alt},
                best_other,
                fasta
            )
        except Exception:
            # FASTA fetch failed
            mnv_var = ''

        return {
            'mnv_proba': best_prob,
            'mnv_inframe': None,
            'mnv_variant': mnv_var
        }
    else:
        # Calculate inframe status for indels
        all_indels = [
            {'ref': ref, 'alt': alt},
            {'ref': best_other['ref'], 'alt': best_other['alt']}
        ]
        inframe = calculate_indel_inframe(all_indels)

        return {
            'mnv_proba': best_prob,
            'mnv_inframe': inframe,
            'mnv_variant': ''
        }


def detect_mnv_probas(
    df: pl.DataFrame,
    bcf_path: str,
    fasta_path: str,
    config: dict
) -> pl.DataFrame:
    """Main entry point for MNV detection.

    Adds 3 columns to DataFrame:
    - mnv_proba: float or None (phasing probability)
    - mnv_inframe: bool or None (for indels only)
    - mnv_variant: str (reconstructed variant for SNVs)

    Args:
        df: Filtered variants DataFrame (one row per variant-sample)
        bcf_path: Path to tabix-indexed BCF file
        fasta_path: Path to indexed FASTA file
        config: MNV config dict with optional 'indel_window' key (default 10)

    Returns:
        DataFrame with 3 additional columns

    Raises:
        ImportError: If cyvcf2 or pysam not installed
        FileNotFoundError: If BCF or FASTA files not found
    """
    _check_dependencies()

    # Get config parameters
    indel_window = config.get('indel_window', 10)

    # Open BCF and FASTA files
    try:
        vcf = cyvcf2.VCF(bcf_path)
    except Exception as e:
        raise FileNotFoundError(f"Cannot open BCF file {bcf_path}: {e}")

    try:
        fasta = pysam.FastaFile(fasta_path)
    except Exception as e:
        raise FileNotFoundError(f"Cannot open FASTA file {fasta_path}: {e}")

    # Process each row
    results = []
    for row in df.iter_rows(named=True):
        mnv_result = process_variant_for_mnv(row, vcf, fasta, indel_window)
        results.append(mnv_result)

    # Add results as new columns
    df = df.with_columns([
        pl.Series('mnv_proba', [r['mnv_proba'] for r in results], dtype=pl.Float64),
        pl.Series('mnv_inframe', [r['mnv_inframe'] for r in results], dtype=pl.Boolean),
        pl.Series('mnv_variant', [r['mnv_variant'] for r in results], dtype=pl.Utf8),
    ])

    # Close files
    fasta.close()

    return df
