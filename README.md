# PyWombat

A CLI tool for processing bcftools tabulated TSV files.

## Installation

This is a UV-managed Python package. To install:

```bash
uv sync
```

## Usage

The main command is `wombat` with a `format` subcommand:

```bash
# Format a bcftools TSV file and print to stdout
wombat format input.tsv

# Format and save to output file
wombat format input.tsv -o output.tsv
wombat format input.tsv --output output.tsv
```

### What does `format` do?

The `format` command processes bcftools tabulated TSV files by:

1. **Expanding the `(null)` column**: This column contains multiple fields in the format `NAME=value` separated by semicolons (e.g., `DP=30;AF=0.5;AC=2`). Each field is extracted into its own column.

2. **Melting sample columns**: After the `(null)` column, there are typically sample columns with headers like `Sample1:GT`, `Sample2:GT`, etc. The tool:
   - Extracts the sample name (the part before the `:` character)
   - Transforms the wide format into long format
   - Creates a `sample` column with the sample names
   - Creates a `sample_value` column with the corresponding values

### Example

**Input:**
```tsv
CHROM	POS	REF	ALT	(null)	Sample1:GT	Sample2:GT
chr1	100	A	T	DP=30;AF=0.5	0/1	1/1
```

**Output:**
```tsv
CHROM	POS	REF	ALT	AC	AF	DP	sample	sample_value
chr1	100	A	T	null	0.5	30	Sample1	0/1
chr1	100	A	T	null	0.5	30	Sample2	1/1
```

## Development

This project uses:
- **UV** for package management
- **Polars** for fast data processing
- **Click** for CLI interface

## Testing

Test files are available in the `tests/` directory:
- `test.tabulated.tsv` - Real bcftools output
- `test_small.tsv` - Small example for quick testing
