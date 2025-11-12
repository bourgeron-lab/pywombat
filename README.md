# PyWombat

A CLI tool for processing bcftools tabulated TSV files.

## Installation

This is a UV-managed Python package. To install:

```bash
uv sync
```

## Usage

The `wombat` command processes bcftools tabulated TSV files:

```bash
# Format a bcftools TSV file and print to stdout
wombat input.tsv

# Format and save to output file (creates output.tsv by default)
wombat input.tsv -o output

# Format and save as parquet
wombat input.tsv -o output -f parquet
wombat input.tsv -o output --format parquet

# Format with pedigree information to add parent genotypes
wombat input.tsv --pedigree pedigree.tsv -o output
```

### What does `wombat` do?

The `wombat` command processes bcftools tabulated TSV files by:

1. **Expanding the `(null)` column**: This column contains multiple fields in the format `NAME=value` separated by semicolons (e.g., `DP=30;AF=0.5;AC=2`). Each field is extracted into its own column.

2. **Preserving the `CSQ` column**: The CSQ (Consequence) column is preserved as-is and not melted, allowing VEP annotations to remain intact.

3. **Melting and splitting sample columns**: After the `(null)` column, there are typically sample columns with values in `GT:DP:GQ:AD` format. The tool:
   - Extracts the sample name (the part before the first `:` character)
   - Transforms the wide format into long format
   - Creates a `sample` column with the sample names
   - Splits the sample values into separate columns:
     - `sample_gt`: Genotype (e.g., 0/1, 1/1)
     - `sample_dp`: Read depth
     - `sample_gq`: Genotype quality
     - `sample_ad`: Allele depth (takes the second value from comma-separated list)
     - `sample_vaf`: Variant allele frequency (calculated as sample_ad / sample_dp)

### Example

**Input:**

```tsv
CHROM POS REF ALT (null) Sample1:GT:Sample1:DP:Sample1:GQ:Sample1:AD Sample2:GT:Sample2:DP:Sample2:GQ:Sample2:AD
chr1 100 A T DP=30;AF=0.5;AC=2 0/1:15:99:5,10 1/1:18:99:0,18
```

**Output:**

```tsv
CHROM POS REF ALT AC AF DP sample sample_gt sample_dp sample_gq sample_ad sample_vaf
chr1 100 A T 2 0.5 30 Sample1 0/1 15 99 10 0.6667
chr1 100 A T 2 0.5 30 Sample2 1/1 18 99 18 1.0
```

Notes:

- The `sample_ad` column contains the second value from the AD field (e.g., from `5,10` it extracts `10`)
- The `sample_vaf` column is the variant allele frequency calculated as `sample_ad / sample_dp`
- By default, output is in TSV format. Use `-f parquet` to output as Parquet files
- The `-o` option specifies an output prefix (e.g., `-o output` creates `output.tsv` or `output.parquet`)

### Pedigree Support

You can provide a pedigree file with the `--pedigree` option to add parent genotype information to the output. This enables trio analysis by including the father's and mother's genotypes for each sample.

**Pedigree File Format:**

The pedigree file should be a tab-separated file with the following columns:

- `FID`: Family ID
- `sample_id`: Sample identifier (matches the sample names in the VCF)
- `FatherBarcode`: Father's sample identifier (use `0` or `-9` if unknown)
- `MotherBarcode`: Mother's sample identifier (use `0` or `-9` if unknown)
- `Sex`: Sex of the sample (optional)
- `Pheno`: Phenotype information (optional)

Example pedigree file:

```tsv
FID sample_id FatherBarcode MotherBarcode Sex Pheno
FAM1 Child1 Father1 Mother1 1 2
FAM1 Father1 0 0 1 1
FAM1 Mother1 0 0 2 1
```

**Output with Pedigree:**

When using `--pedigree`, the output will include additional columns for each parent:

- `father_gt`, `father_dp`, `father_gq`, `father_ad`, `father_vaf`: Father's genotype information
- `mother_gt`, `mother_dp`, `mother_gq`, `mother_ad`, `mother_vaf`: Mother's genotype information

These columns will contain the parent's genotype data for the same variant, allowing you to analyze inheritance patterns.

## Development

This project uses:

- **UV** for package management
- **Polars** for fast data processing
- **Click** for CLI interface

## Testing

Test files are available in the `tests/` directory:

- `test.tabulated.tsv` - Real bcftools output
- `test_small.tsv` - Small example for quick testing
