"""Pytest fixtures for pywombat tests."""

import gzip
from pathlib import Path

import pytest


@pytest.fixture
def tmp_test_dir(tmp_path):
    """Create a temporary directory for test files."""
    return tmp_path


@pytest.fixture
def small_tsv_content():
    """Return content for a small test TSV file."""
    # Minimal bcftools-style TSV with (null) INFO column and sample columns
    header = (
        "#CHROM\tPOS\tREF\tALT\tFILTER\t(null)\t"
        "Sample1:GT:Sample1:DP:Sample1:GQ:Sample1:AD\t"
        "Sample2:GT:Sample2:DP:Sample2:GQ:Sample2:AD\t"
        "Sample3:GT:Sample3:DP:Sample3:GQ:Sample3:AD"
    )

    # Row 1: Variant with VEP annotations
    row1 = (
        "chr1\t100\tA\tT\t.\t"
        "AF=0.001;VEP_SYMBOL=GENE1;VEP_IMPACT=HIGH;VEP_CANONICAL=YES\t"
        "0/1:20:99:15,5\t0/0:25:99:25,0\t0/0:30:99:30,0"
    )

    # Row 2: Another variant
    row2 = (
        "chr1\t200\tG\tC\t.\t"
        "AF=0.01;VEP_SYMBOL=GENE2;VEP_IMPACT=MODERATE;VEP_CANONICAL=YES\t"
        "0/0:22:99:22,0\t0/1:28:99:20,8\t0/1:35:99:25,10"
    )

    # Row 3: Rare variant
    row3 = (
        "chr2\t500\tC\tG\t.\t"
        "AF=0.0001;VEP_SYMBOL=GENE3;VEP_IMPACT=LOW;VEP_CANONICAL=YES\t"
        "0/1:18:99:14,4\t0/0:20:99:20,0\t0/1:25:99:20,5"
    )

    return "\n".join([header, row1, row2, row3]) + "\n"


@pytest.fixture
def small_tsv(tmp_test_dir, small_tsv_content):
    """Create a small test TSV file."""
    tsv_path = tmp_test_dir / "test_input.tsv"
    tsv_path.write_text(small_tsv_content)
    return tsv_path


@pytest.fixture
def small_tsv_gz(tmp_test_dir, small_tsv_content):
    """Create a small gzipped test TSV file."""
    gz_path = tmp_test_dir / "test_input.tsv.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write(small_tsv_content)
    return gz_path


@pytest.fixture
def small_pedigree_content():
    """Return content for a small test pedigree file."""
    header = "FID\tsample_id\tFatherBarcode\tMotherBarcode\tSex\tPheno"
    # Simple trio: father (Sample1), mother (Sample2), child (Sample3)
    father = "FAM1\tSample1\t0\t0\t1\t1"
    mother = "FAM1\tSample2\t0\t0\t2\t1"
    child = "FAM1\tSample3\tSample1\tSample2\t1\t2"
    return "\n".join([header, father, mother, child]) + "\n"


@pytest.fixture
def small_pedigree(tmp_test_dir, small_pedigree_content):
    """Create a small test pedigree file."""
    ped_path = tmp_test_dir / "test_pedigree.tsv"
    ped_path.write_text(small_pedigree_content)
    return ped_path


@pytest.fixture
def filter_config_content():
    """Return content for a simple filter config."""
    return """
quality:
  min_dp: 10
  min_gq: 20
  filter_no_alt_allele: true
"""


@pytest.fixture
def filter_config(tmp_test_dir, filter_config_content):
    """Create a filter config file."""
    config_path = tmp_test_dir / "filter_config.yml"
    config_path.write_text(filter_config_content)
    return config_path


@pytest.fixture
def real_test_data():
    """Return paths to real test data if available."""
    tests_dir = Path(__file__).parent
    tsv_gz = tests_dir / "C0733-011-068.rare.VEP_GRCh38_112_standard.annotated.tsv.gz"
    pedigree = tests_dir / "C0733-011-068.pedigree.tsv"

    if tsv_gz.exists() and pedigree.exists():
        return {"tsv_gz": tsv_gz, "pedigree": pedigree}
    return None


@pytest.fixture
def fam_93_data():
    """Return paths to fam_93 test data if available."""
    tests_dir = Path(__file__).parent / "fam_93"
    tsv_gz = tests_dir / "C0733-011-093.rare.VEP_GRCh38_112_standard.annotated.tsv.gz"
    pedigree = tests_dir / "C0733-011-093.pedigree.tsv"
    config = tests_dir / "rare_variants_high_impact.yml"

    if tsv_gz.exists() and pedigree.exists() and config.exists():
        return {"tsv_gz": tsv_gz, "pedigree": pedigree, "config": config}
    return None
