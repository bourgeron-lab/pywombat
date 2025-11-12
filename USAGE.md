"""
Wombat CLI - BCFtools TSV Formatter
====================================

This module provides a command-line tool for reformatting bcftools tabulated TSV files.

Key Features:
-------------
1. Expands the (null) column containing semicolon-separated NAME=value pairs
2. Melts sample columns from wide to long format
3. Uses Polars for fast, efficient data processing

Usage Examples:
--------------

Basic usage:
    $ wombat format input.tsv -o output.tsv

With verbose output:
    $ wombat format input.tsv -o output.tsv --verbose

Output to stdout:
    $ wombat format input.tsv > output.tsv
    $ wombat format input.tsv | head -20

Pipeline usage:
    $ wombat format input.tsv | wombat format - -o final.tsv


Input Format:
------------
The input TSV file should be a bcftools tabulated output with:
- Variant information columns (CHROM, POS, REF, ALT, etc.)
- A (null) column with format: "FIELD1=value1;FIELD2=value2;..."
- Sample columns with format: "SampleName:FORMAT" (e.g., "Sample1:GT")

Example Input:
CHROM  POS  REF  ALT  (null)              Sample1:GT  Sample2:GT  Sample3:GT
chr1   100  A    T    DP=30;AF=0.5;AC=2   0/1         1/1         0/0
chr1   200  G    C    DP=45;AF=0.25;AC=1  0/0         0/1         0/0


Output Format:
-------------
The output TSV file will have:
- Original variant columns preserved
- (null) column expanded into individual field columns
- One row per sample per variant
- 'sample' column with sample names
- 'sample_value' column with the genotype/value for that sample

Example Output:
CHROM  POS  REF  ALT  AC    AF    DP  sample   sample_value
chr1   100  A    T    2     0.5   30  Sample1  0/1
chr1   100  A    T    2     0.5   30  Sample2  1/1
chr1   100  A    T    2     0.5   30  Sample3  0/0
chr1   200  G    C    1     0.25  45  Sample1  0/0
chr1   200  G    C    1     0.25  45  Sample2  0/1
chr1   200  G    C    1     0.25  45  Sample3  0/0


Implementation Details:
----------------------
The transformation happens in three steps:

Step 1 - Parse (null) column:
  - Extract all unique field names (DP, AF, AC, etc.)
  - Use regex to extract values: "FIELD=([^;]+)"
  - Create new columns for each field

Step 2 - Identify sample columns:
  - Find columns after (null)
  - Extract sample names (part before ':')
  - Handle columns without ':' as-is

Step 3 - Melt to long format:
  - Keep variant info as id_vars
  - Melt sample columns
  - Create 'sample' and 'sample_value' columns


Performance Notes:
-----------------
- Uses Polars for fast CSV reading and processing
- Handles large files efficiently
- Memory usage scales with number of rows × samples
- For very large files (>1GB), consider processing in chunks

Error Handling:
--------------
The tool will fail with helpful messages if:
- Input file doesn't exist
- No (null) column is found
- File is not valid TSV format
- Memory is insufficient


Dependencies:
------------
- polars >= 0.19.0: Fast DataFrame library
- click >= 8.1.0: CLI framework
"""

# This file serves as documentation only
# See cli.py for implementation
