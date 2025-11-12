#!/usr/bin/env python3
"""Test script for wombat CLI."""

import sys

sys.path.insert(0, "/Users/fcliquet/Workspace/pywombat/src")

import polars as pl

from pywombat.cli import format_bcftools_tsv

# Create test data
test_data = {
    "CHROM": ["chr1", "chr1", "chr2"],
    "POS": [100, 200, 150],
    "REF": ["A", "G", "C"],
    "ALT": ["T", "C", "G"],
    "(null)": ["DP=30;AF=0.5;AC=2", "DP=45;AF=0.25;AC=1", "DP=60;AF=0.75;AC=3"],
    "Sample1:GT": ["0/1", "0/0", "1/1"],
    "Sample2:GT": ["1/1", "0/1", "0/1"],
    "Sample3:GT": ["0/0", "0/0", "1/1"],
}

df = pl.DataFrame(test_data)
print("Input DataFrame:")
print(df)
print("\n" + "=" * 80 + "\n")

result = format_bcftools_tsv(df)
print("Output DataFrame:")
print(result)
