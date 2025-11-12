# Wombat Quick Start

## Installation

```bash
cd /Users/fcliquet/Workspace/pywombat
uv sync
```

## Basic Usage

```bash
# Format a bcftools TSV file
uv run wombat format tests/test.tabulated.tsv -o output.tsv

# With verbose output to see progress
uv run wombat format tests/test.tabulated.tsv -o output.tsv --verbose

# Output to stdout
uv run wombat format tests/test.tabulated.tsv

# View help
uv run wombat --help
uv run wombat format --help
```

## What It Does

### Input (Wide Format)

```
CHROM  POS  REF  ALT  (null)              Sample1:GT  Sample2:GT
chr1   100  A    T    DP=30;AF=0.5;AC=2   0/1         1/1
```

### Output (Long Format)

```
CHROM  POS  REF  ALT  AC  AF   DP  sample   sample_value
chr1   100  A    T    2   0.5  30  Sample1  0/1
chr1   100  A    T    2   0.5  30  Sample2  1/1
```

## Features

1. **Expands (null) column**: Splits `DP=30;AF=0.5;AC=2` into separate columns
2. **Melts samples**: Converts wide sample columns to long format
3. **Fast processing**: Uses Polars for efficient data handling

## Testing

```bash
# Run the test suite
uv run python tests/test_format.py

# Test with small example
uv run wombat format tests/test_small.tsv -o tests/output_small.tsv

# Test with real data
uv run wombat format tests/test.tabulated.tsv -o tests/output.tsv
```

## Project Structure

```
pywombat/
├── pyproject.toml          # Project configuration with dependencies
├── README.md               # Main documentation
├── USAGE.md               # Detailed usage guide
├── src/
│   └── pywombat/
│       ├── __init__.py    # Package init
│       └── cli.py         # Main CLI implementation
└── tests/
    ├── test_format.py     # Unit tests
    ├── test_small.tsv     # Small test file
    └── test.tabulated.tsv # Real bcftools output
```

## Dependencies

- **polars** (>=0.19.0): Fast DataFrame library for data processing
- **click** (>=8.1.0): Command-line interface framework

## Troubleshooting

If you encounter issues:

1. Make sure dependencies are installed: `uv sync`
2. Check that input file is valid TSV with a `(null)` column
3. Use `--verbose` flag to see processing details
4. For large files, ensure sufficient memory is available
