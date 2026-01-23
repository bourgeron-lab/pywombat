# PyWombat Development Guide

## Project Overview

PyWombat is a bioinformatics CLI tool that transforms bcftools tabulated TSV files into analysis-ready long-format data. It specializes in:
- Expanding VCF INFO fields from `(null)` column into individual columns
- Converting wide-format sample data to long-format (melting)
- Extracting and calculating genotype metrics (GT, DP, GQ, AD, VAF)
- De novo mutation (DNM) detection with sex-chromosome awareness
- Trio/family analysis with pedigree-based parent genotype joining

## Architecture

### Single-File Monolith
All functionality lives in [src/pywombat/cli.py](src/pywombat/cli.py) (~2100 lines). This is intentional for deployment simplicity with `uvx` (one-command installation-free execution).

### Core Data Flow
1. **Input**: bcftools tabulated TSV (wide format, one row per variant)
2. **Transform Pipeline**:
   - `format_bcftools_tsv_lazy()` → Expand INFO, melt samples, parse genotypes
   - `apply_filters_lazy()` → Apply quality/expression filters OR DNM detection
3. **Output**: Long-format TSV/Parquet (one row per variant-sample)

### Key Functions
- **`format_bcftools_tsv()`**: Eager transformation (collects data)
- **`format_bcftools_tsv_lazy()`**: Streaming wrapper using Polars LazyFrame
- **`apply_filters_lazy()`**: Conditional filter dispatcher (quality/expression OR DNM)
- **`apply_de_novo_filter()`**: Complex vectorized DNM logic with PAR region handling
- **`parse_impact_filter_expression()`**: Expression parser for flexible YAML-based filtering

## Development Workflow

### Setup & Running
```bash
# Install dependencies (uses uv package manager)
uv sync

# Run CLI locally
uv run wombat input.tsv -o output

# Production usage (no installation)
uvx pywombat input.tsv -o output
```

### Testing
- Test data: `tests/C0733-011-068.*` files (real variant data)
- Compare script: [tests/compare_dnm_results.py](tests/compare_dnm_results.py) validates DNM detection
- No formal test framework yet—manual validation against reference outputs in `tests/check/`

## Critical Domain Logic

### Genotype Parsing
- Input format: `GT:DP:GQ:AD` (e.g., `0/1:25:99:10,15`)
- Extracted fields:
  - `sample_gt`: Genotype (0/1, 1/1, etc.)
  - `sample_dp`: Total depth (25)
  - `sample_gq`: Quality (99)
  - `sample_ad`: **Second value** from AD (15 = alternate allele depth)
  - `sample_vaf`: Calculated as `sample_ad / sample_dp` (0.6)

### De Novo Mutation Detection
Sex-chromosome-aware logic in `apply_de_novo_filter()`:
- **Autosomes**: Both parents must have VAF < 2% (reference genotype)
- **X chromosome (male proband)**: Mother must be reference; father is hemizygous (ignored)
- **Y chromosome**: Only father checked (mother doesn't have Y)
- **PAR regions**: Treated as autosomal (both parents checked)
- **Hemizygous variants** (X/Y in males): Require VAF ≥ 85% (not ~50% like het)

PAR regions defined in [examples/de_novo_mutations.yml](examples/de_novo_mutations.yml):
```yaml
par_regions:
  GRCh38:
    PAR1: { chrom: "X", start: 10001, end: 2781479 }
    PAR2: { chrom: "X", start: 155701383, end: 156030895 }
```

## Conventions & Patterns

### Polars Usage
- **Lazy evaluation**: Use `scan_csv()` → `LazyFrame` → `sink_csv()` for streaming
- **Streaming mode**: `collect(streaming=True)` for memory-efficient processing
- **Schema overrides**: Force string type for sample IDs to prevent numeric inference:
  ```python
  schema_overrides = {col: pl.Utf8 for col in ["sample", "FatherBarcode", ...]}
  ```

### YAML Configuration Files
Filter configs in [examples/](examples/) define:
- `quality`: Thresholds (DP, GQ, VAF ranges for het/hom)
- `expression`: Polars-compatible filter expressions (see `parse_impact_filter_expression()`)
- `dnm`: De novo detection parameters (parent thresholds, PAR regions)

Example expression syntax:
```yaml
expression: "VEP_IMPACT = HIGH & fafmax_faf95_max_genomes <= 0.001"
```
Supports: `=`, `!=`, `<=`, `>=`, `<`, `>`, `&`, `|`, `()`, NULL, NaN

### Error Handling
- Use `click.echo(..., err=True)` for verbose messages (stderr)
- Raise `click.Abort()` instead of generic exceptions for CLI errors
- Validate pedigree/config requirements early (e.g., DNM mode requires pedigree)

## Common Pitfalls

1. **VAF Calculation**: Always use `sample_ad / sample_dp`, not first AD value
2. **Genotype Filtering**: Use `str.contains("1")` not regex—handles 0/1, 1/0, 1/1, 1/2
3. **Parent Joining**: Pedigree must have `sample_id`, `FatherBarcode`, `MotherBarcode` columns
4. **Sex Normalization**: DNM filter accepts "1"/"2" or "M"/"F" for sex; normalizes to uppercase
5. **Categorical Types**: Convert `#CHROM` from Categorical to Utf8 before string operations

## File References

- Main CLI: [src/pywombat/cli.py](src/pywombat/cli.py)
- Config examples: [examples/de_novo_mutations.yml](examples/de_novo_mutations.yml), [examples/rare_variants_high_impact.yml](examples/rare_variants_high_impact.yml)
- Test data: [tests/C0733-011-068.pedigree.tsv](tests/C0733-011-068.pedigree.tsv)
- Documentation: [README.md](README.md), [QUICKSTART.md](QUICKSTART.md)

## Adding New Features

When adding filters:
1. Update YAML schema in config examples
2. Extend `apply_filters_lazy()` or create specialized function (like `apply_de_novo_filter()`)
3. Use lazy operations when possible; collect only if vectorization requires it
4. Add verbose logging with `--verbose` flag support

When modifying transforms:
- Keep `format_bcftools_tsv()` and `format_bcftools_tsv_lazy()` in sync
- Test with real data from `tests/` directory
- Verify output with `compare_dnm_results.py` for DNM changes
