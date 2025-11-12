#!/usr/bin/env python3
"""
Example script demonstrating the wombat format command.
This shows the transformation that happens to bcftools TSV files.
"""

import polars as pl

# Create example input data
print("=" * 80)
print("EXAMPLE: BCFTools TSV Formatting")
print("=" * 80)

print("\n1. ORIGINAL DATA (bcftools output):")
print("-" * 80)

original_data = {
    "CHROM": ["chr1", "chr1", "chr2"],
    "POS": [100, 200, 150],
    "REF": ["A", "G", "C"],
    "ALT": ["T", "C", "G"],
    "(null)": ["DP=30;AF=0.5;AC=2", "DP=45;AF=0.25;AC=1", "DP=60;AF=0.75;AC=3"],
    "Sample1:GT": ["0/1", "0/0", "1/1"],
    "Sample2:GT": ["1/1", "0/1", "0/1"],
    "Sample3:GT": ["0/0", "0/0", "1/1"],
}

df_original = pl.DataFrame(original_data)
print(df_original)

print("\n2. AFTER EXPANDING (null) COLUMN:")
print("-" * 80)
print("The (null) column 'DP=30;AF=0.5;AC=2' becomes three columns:")
print("  - DP: 30")
print("  - AF: 0.5")
print("  - AC: 2")

print("\n3. AFTER MELTING SAMPLE COLUMNS:")
print("-" * 80)
print("Each row is duplicated for each sample:")
print("  - Sample1:GT, Sample2:GT, Sample3:GT")
print("Becomes:")
print("  - 3 rows with 'sample' column containing: Sample1, Sample2, Sample3")
print("  - 'sample_value' column containing the values: 0/1, 1/1, 0/0")

print("\n4. FINAL RESULT:")
print("-" * 80)
print("Each variant position now has one row per sample,")
print("with all fields from (null) expanded into separate columns.")

print("\n" + "=" * 80)
print("To run this transformation on your data:")
print("  wombat format input.tsv -o output.tsv")
print("  wombat format input.tsv -o output.tsv --verbose")
print("=" * 80)
