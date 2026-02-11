# CLAUDE.md - PyWombat Developer Guide for AI Assistants

This document provides context and guidelines for AI assistants working on the pywombat codebase.

## Project Overview

**PyWombat** is a high-performance CLI tool for processing and filtering bcftools tabulated TSV files with advanced filtering capabilities, pedigree support, and de novo mutation (DNM) detection.

- **Version**: 1.5.1 (current)
- **Python**: 3.12+
- **License**: MIT
- **Repository**: https://github.com/bourgeron-lab/pywombat

### Key Features

- **Two-step workflow**: prepare (TSV→Parquet) + filter for optimal performance
- **Memory optimization**: Up to 95% reduction with expression filtering, 88% with per-chromosome DNM processing
- **Flexible filtering**: Quality filters, expression-based filters, DNM detection, MNV detection
- **Pedigree support**: Parent genotype integration and inheritance pattern analysis
- **Multiple output formats**: TSV, TSV.gz, Parquet

### Development Environment

**IMPORTANT**: This project uses [UV](https://github.com/astral-sh/uv) for package management and Python execution.

**Key UV Commands**:

- **Add dependencies**: `uv add <package>` (NOT `pip install`)
- **Install project**: `uv pip install -e .` or `uv pip install -e ".[dev,mnv]"`
- **Run commands**: `uv run wombat <args>` or `uv run pytest`
- **Run Python**: `uv run python script.py`

**Why UV?**

- Fast dependency resolution and installation
- Reproducible builds with `uv.lock`
- Consistent development environment across contributors

**For AI Assistants**: When modifying this project:

- ✅ Use `uv add <package>` to add new dependencies
- ✅ Use `uv run <command>` to execute scripts/tests
- ❌ DO NOT use `pip install` directly (breaks lock file consistency)
- ❌ DO NOT run `python` or `pytest` without `uv run` prefix

**Note for end users**: When using the published package, users can run via:

- `uvx pywombat <args>` (recommended, no installation needed)
- `pipx install pywombat && wombat <args>` (isolated install)
- But for **development/modification**, always use `uv run wombat <args>`

## Project Structure

```
pywombat/
├── src/pywombat/
│   ├── __init__.py           # Package initialization
│   ├── cli.py                # Main CLI implementation
│   ├── mnv.py                # MNV detection module (optional feature)
│   └── vep.py                # VEP REST API annotation module
├── tests/
│   ├── conftest.py           # Pytest fixtures
│   ├── test_cli.py           # CLI structure tests
│   ├── test_filter.py        # Filter command tests (includes DNM, homalt)
│   ├── test_prepare.py       # Prepare command tests
│   ├── test_mnv.py           # MNV detection tests (requires optional deps)
│   ├── test_vep.py           # VEP annotation tests
│   └── test_expression_parser.py  # Expression parser tests
├── examples/
│   ├── README.md             # Configuration examples guide
│   ├── rare_high_impact.yml  # High-impact variants
│   ├── rare_homalt.yml       # Homozygous alternative filtering
│   ├── rare_coding.yml       # All coding variants
│   ├── rare_missense.yml     # Missense with pathogenicity scores
│   ├── rare_spliceai.yml     # Splicing variants
│   ├── de_novo_mutations.yml # DNM detection config
│   └── rare_homalt_mnv.yml   # Homozygous + MNV detection
├── pyproject.toml            # Project metadata and dependencies
├── README.md                 # User-facing documentation
├── CHANGELOG.md              # Version history
└── .github/workflows/
    └── publish.yml          # PyPI publishing workflow
```

## Dependencies

### Core Dependencies
- **polars** >=0.19.0 - Fast DataFrame operations (primary data processing)
- **pyarrow** >=14.0.0 - Parquet file I/O
- **click** >=8.1.0 - CLI framework
- **pyyaml** >=6.0 - Configuration file parsing
- **tqdm** >=4.67.1 - Progress bars

### Optional Dependencies (extras: `mnv`)
- **cyvcf2** >=0.31.0 - VCF/BCF parsing for MNV detection
- **pysam** >=0.22.0 - SAM/BAM/FASTA file handling
- **scipy** >=1.11.0 - Statistical functions (phasing probability)

### Development Dependencies (extras: `dev`)
- **pytest** >=7.0.0 - Testing framework
- **pytest-cov** >=4.0.0 - Code coverage

## Core Functionality

### Command 1: `wombat prepare`

**Purpose**: Converts TSV/TSV.gz to optimized Parquet format with pre-expanded INFO fields.

**Location**: cli.py, line 27

**Process**:
1. First pass: Discovers all INFO fields (key-value pairs and boolean flags)
2. Second pass: Processes file in chunks (default 50k rows)
3. Expands INFO fields from `(null)` column into separate columns
4. Applies memory-efficient dtypes (Categorical, UInt32)
5. Writes to Parquet format for columnar access

**Key functions**:
- `prepare_cmd()` - Command entry point
- `prepare_parquet()` (line 86) - Main preprocessing function
- `_process_chunk()` (line 218) - Processes chunks of TSV data

**Usage**:
```bash
uv run wombat prepare input.tsv -o prepared.parquet --chunk-size 50000 -v
```

### Command 2: `wombat filter`

**Purpose**: Transforms and filters variant data (works with TSV or Parquet input).

**Location**: cli.py, line 372

**Command signature**:
```python
@cli.command("filter")
def filter_cmd(
    input_file: Path,
    output: Optional[str],
    output_format: str,  # tsv, tsv.gz, or parquet
    verbose: bool,
    pedigree: Optional[Path],
    filter_config: Optional[Path],
    debug: Optional[str],  # chr:pos format
    bcf: Optional[Path],    # For MNV detection
    fasta: Optional[Path],  # For MNV detection
)
```

**Key optimizations**:
- **For Parquet + DNM mode**: Per-chromosome processing (88% memory reduction)
- **For Parquet + expression filters**: Applies filters BEFORE melting (95%+ memory reduction)
- **Auto-detection**: Automatically detects input format (TSV, TSV.gz, or Parquet)

**Data Flow**:
```
Input (TSV/Parquet)
  ↓
[TSV only] Expand INFO fields → Melt samples
[Parquet] Apply expression filters → Melt samples
  ↓
Split sample values (GT:DP:GQ:AD)
  ↓
Calculate VAF (variant allele frequency)
  ↓
Add parent genotypes (if pedigree provided)
  ↓
Apply quality filters
  ↓
Apply DNM filters (if enabled)
  ↓
Apply MNV detection (if enabled)
  ↓
Annotate MNV via VEP API (if mnv.annotate enabled)
  ↓
Write output
```

### Command 3: `wombat vep`

**Purpose**: Annotate a single variant using Ensembl VEP REST API (GRCh38).

**Location**: cli.py (registered as `@cli.command("vep")`)

**Module**: `vep.py` — contains all VEP-related functions (no external dependencies beyond stdlib).

**Usage**:
```bash
uv run wombat vep chr17:43094433:G:A          # Canonical transcript only
uv run wombat vep chr17:43094433:G:A --all    # All transcripts
uv run wombat vep chr17:43094433:G:A -v       # Verbose
```

**Key functions** (in `vep.py`):
- `parse_variant()` — Parses `chr:pos:ref:alt` format
- `query_vep()` — POST to Ensembl REST API (batch support, up to 200 variants)
- `extract_annotations()` — Extracts per-transcript annotations from VEP response
- `format_annotations()` — Formats annotations for CLI display
- `annotate_mnv_variants()` — Annotates MNV candidates in a DataFrame (used by `filter` command when `mnv.annotate: true`)
- `_find_transcript_annotation()` — Finds a specific transcript in a VEP response by ENST ID

**VEP API parameters**: `canonical=1, LoF=1, hgvs=1, protein=1, mane=1, numbers=1`

**LOFTEE support**: Enabled via `LoF=1`. Returns `LoF` (HC/LC/OS), `LoF_filter`, `LoF_flags`, `LoF_info` for loss-of-function variants.

**MNV annotation** (via `mnv.annotate: true` in YAML config):
- After MNV detection, queries VEP for reconstructed SNV MNV variants
- Matches response to the original row's `VEP_Feature` transcript
- Appends new rows with MNV VEP annotations to the output DataFrame
- `mnv_proba` column for new rows = `1.0 - original_value`

## Filter Types

### 1. Quality Filters

**Function**: `apply_quality_filters()` (cli.py, line 772)

**Applies depth, quality, and VAF thresholds based on genotype type.**

**Available filters**:
- `filter_no_alt_allele` (bool) - Exclude 0/0 genotypes
- `homalt_only` (bool) - Keep only 1/1 genotypes (recessive analysis)
- `sample_dp_min` (int) - Minimum read depth
- `sample_gq_min` (int) - Minimum genotype quality
- `sample_vaf_het_min` (float) - Min VAF for 0/1 (typically 0.25)
- `sample_vaf_het_max` (float) - Max VAF for 0/1 (typically 0.75)
- `sample_vaf_homalt_min` (float) - Min VAF for 1/1 (typically 0.85-0.98)
- `sample_vaf_hom_ref_max` (float) - Max VAF for 0/0
- `apply_to_parents` (bool) - Apply same filters to parents

**Genotype detection pattern**:
```python
is_het = (pl.col("sample_gt").str.count_matches("1") == 1) & ~pl.col("sample_gt").str.contains("2")
is_hom_alt = pl.col("sample_gt") == "1/1"
is_hom_ref = pl.col("sample_gt") == "0/0"
```

### 2. Expression-Based Filters

**Function**: `parse_impact_filter_expression()` (cli.py, line 1327)

**Parses YAML filter expressions into Polars expressions.**

**Supported operators**:
- **Comparison**: `=`, `!=`, `<`, `>`, `<=`, `>=`
- **Logical**: `&` (AND), `|` (OR)
- **Grouping**: `(`, `)`
- **Special predicates**:
  - `is_snv` - REF and ALT are single nucleotides (A/C/G/T)
  - `is_indel` - NOT a SNV
  - `contains` - Case-insensitive substring matching
  - `is_empty` - Null, empty string, or "."
  - `not_empty` - NOT (null/empty/".")
  - `= null` / `!= null` - Null checks

**Column name support**: Supports dots in column names (e.g., `cadd_v1.7`, `faf95.joint`)

**Example expression**:
```yaml
expression: >
  VEP_CANONICAL = YES &
  (VEP_IMPACT = HIGH | VEP_IMPACT = MODERATE) &
  (fafmax_faf95_max_genomes = null | fafmax_faf95_max_genomes <= 0.001) &
  ((is_snv & AQ > 15) | (is_indel & AQ > 22)) &
  genomes_filters is_empty
```

### 3. De Novo Mutation (DNM) Filter

**Function**: `apply_de_novo_filter()` (cli.py, line 1034)

**Detects de novo mutations in trios/families with sex-chromosome awareness.**

**DNM detection logic**:
```
Autosomes & PAR regions:
  - Both parents must be 0/0 (homozygous reference)
  - Parent VAF < threshold (typically 0.02 = 2%)
  - Proband must have alternate allele

X chromosome (males, non-PAR):
  - Mother must be 0/0
  - Father not informative

Y chromosome (males):
  - Father must be 0/0
  - Mother N/A

Hemizygous variants:
  - VAF must be ≥ 85% (configurable)
```

**Configuration example**:
```yaml
dnm:
  enabled: true
  parent_dp_min: 10
  parent_gq_min: 18
  parent_vaf_max: 0.02
  sample_vaf_hemizygous_min: 0.85
  fafmax_faf95_max_genomes_max: 0.001
  genomes_filters_pass_only: true
  par_regions:
    grch38:
      PAR1: {chrom: X, start: 10000, end: 2781479}
      PAR2: {chrom: X, start: 155701383, end: 156030895}
```

### 4. MNV (Multi-Nucleotide Variant) Detection

**Module**: mnv.py

**Purpose**: Identifies cases where multiple nearby variants may combine to create different functional effects.

**Key functions**:
- `calculate_phasing_probability()` - Calculates P(in-phase) using allelic depth Gaussian model
- `calculate_indel_inframe()` - Determines if indel cluster maintains reading frame (sum % 3 == 0)
- `reconstruct_snv_variant()` - Rebuilds MNV notation using reference genome
- `get_window_variants()` - Queries BCF for variants in genomic windows

**Detection windows**:
- **SNVs**: pos-2 to pos+2 (5bp total window)
- **INDELs**: pos-window to pos+length+window (configurable, default 10bp)

**Output columns**:
- `mnv_proba` - Phasing probability (0.0-1.0) or null
- `mnv_inframe` - True/False (indels only) or null
- `mnv_variant` - Reconstructed variant notation (SNVs only)

**Dependencies**: Requires `cyvcf2`, `pysam`, `scipy` and command-line options `--bcf` and `--fasta`

## Code Patterns

### Polars Lazy vs Eager Evaluation

**Best practice**: Use lazy evaluation to push filters down before `.collect()`

```python
# Good: Lazy evaluation
lazy_df = pl.scan_parquet(input_file)
lazy_df = lazy_df.filter(condition)  # Applied lazily
result = lazy_df.collect()  # Only then materialize

# Avoid: Eager then filter
df = pl.scan_parquet(input_file).collect()  # Loads everything
df = df.filter(condition)
```

**Early filtering optimization** (cli.py, lines 509-564):
```python
# For Parquet input, apply expression filters BEFORE melting
if filter_config_data and "expression" in filter_config_data:
    lazy_df = lazy_df.filter(filter_expr)  # 95%+ memory reduction

# Then melt samples
lazy_df = format_bcftools_tsv_minimal(lazy_df, pedigree_df)
```

### Column Detection Patterns

**Sample columns** are detected by pattern matching:
```python
# Sample columns contain ":" in the name and match GT:DP:GQ:AD format
# Examples: "Sample1:GT:Sample1:DP:Sample1:GQ:Sample1:AD"
```

**INFO field extraction** from `(null)` column:
```python
# Key-value pairs: "DP=30;AF=0.5;AC=2"
for pair in null_value.split(";"):
    if "=" in pair:
        field_name, field_value = pair.split("=", 1)

# Boolean flags: "PASS;DB;SOMATIC"
if flag.strip() and "=" not in flag:
    all_flags.add(flag.strip())
```

### Genotype Parsing

**Format**: `GT:DP:GQ:AD` (e.g., "0/1:15:99:5,10")

```python
# Extraction pattern
sample_gt = "0/1"        # Genotype
sample_dp = "15"         # Depth
sample_gq = "99"         # Quality
sample_ad = "10"         # Alt allele depth (from "5,10")
sample_vaf = 10/15       # Calculated: alt_depth / total_depth
```

### Pedigree Handling

**File format** (tab-separated):
```
FID  sample_id  FatherBarcode  MotherBarcode  Sex  Pheno
FAM1 Child1     Father1        Mother1        1    2
```

**Function**: `add_parent_genotypes()` (cli.py, line 1799)

**Process**:
1. Join pedigree to variant dataframe on `sample_id`
2. Extract parent sample IDs from `FatherBarcode` and `MotherBarcode`
3. For each variant-sample row, find matching father/mother rows
4. Add columns: `father_gt`, `father_dp`, `father_gq`, `father_ad`, `father_vaf`
5. Add columns: `mother_gt`, `mother_dp`, `mother_gq`, `mother_ad`, `mother_vaf`

## Configuration

### YAML Structure

```yaml
quality:
  # Quality thresholds
  sample_dp_min: 10
  sample_gq_min: 19
  sample_vaf_het_min: 0.25
  sample_vaf_het_max: 0.75
  sample_vaf_homalt_min: 0.85
  homalt_only: false
  filter_no_alt_allele: true
  apply_to_parents: false

expression: |
  # Expression-based filter
  VEP_IMPACT = HIGH &
  VEP_CANONICAL = YES &
  gnomad_AF < 0.001

dnm:
  # De novo mutation detection (optional)
  enabled: false
  parent_dp_min: 10
  parent_gq_min: 18
  # ... more options

mnv:
  # Multi-nucleotide variant detection (optional)
  candidate: false
  indel_window: 10
  only: false          # Keep only MNV candidates in output
  annotate: false      # Re-annotate MNV candidates via VEP API
```

**Configuration loading**: `load_filter_config()` (cli.py, line 720)

Returns dict with keys: `quality`, `expression`, `dnm`, `mnv`

## Input/Output Formats

### Input: bcftools Tabulated TSV

```
#CHROM  POS  REF  ALT  FILTER  (null)              Sample1:GT:Sample1:DP:...
chr1    100  A    T    .       DP=30;AF=0.5;PASS   0/1:15:99:5,10
```

**Required format**:
- Use `bcftools +split-vep` plugin to expand VEP annotations
- Sample columns must have format: `SampleName:GT:DP:GQ:AD`

### Input: Parquet (from `wombat prepare`)

- Same columns as TSV but in columnar format
- Pre-expanded INFO fields as separate columns
- Optimized dtypes (Categorical, UInt32)
- Enables lazy evaluation and push-down predicates

### Output Formats

**TSV** (default):
```bash
uv run wombat filter input.tsv -o output  # Creates: output.tsv
```

**Gzipped TSV**:
```bash
uv run wombat filter input.tsv -o output -f tsv.gz  # Creates: output.tsv.gz
```

**Parquet**:
```bash
uv run wombat filter input.tsv -o output -f parquet  # Creates: output.parquet
```

### Output Columns

**Standard columns**:
- `#CHROM`, `POS`, `REF`, `ALT`, `FILTER` - Variant info
- All extracted INFO fields (e.g., `DP`, `AF`, `VEP_SYMBOL`, etc.)
- Boolean flags from INFO (e.g., `PASS`, `DB`)

**Sample columns**:
- `sample`, `sample_gt`, `sample_dp`, `sample_gq`, `sample_ad`, `sample_vaf`

**Parent columns** (if pedigree provided):
- `father_gt`, `father_dp`, `father_gq`, `father_ad`, `father_vaf`
- `mother_gt`, `mother_dp`, `mother_gq`, `mother_ad`, `mother_vaf`

**MNV columns** (if MNV detection enabled):
- `mnv_proba`, `mnv_inframe`, `mnv_variant`

## Testing

### Test Structure

**Configuration** (pyproject.toml):
```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = ["test_*.py"]
python_functions = ["test_*"]
```

**Test fixtures** (conftest.py):
- `tmp_test_dir` - Temporary directory
- `small_tsv` - Synthetic test TSV file
- `small_pedigree` - Test pedigree file

**Test coverage**:
- **test_cli.py** - CLI structure (5 tests)
- **test_filter.py** - Filter functionality (29 tests)
- **test_prepare.py** - Prepare command (7 tests)
- **test_mnv.py** - MNV detection (24 tests, requires optional deps)
- **test_expression_parser.py** - Expression parsing (9 tests)

### Running Tests

```bash
# All tests
uv run pytest

# Specific file
uv run pytest tests/test_filter.py -v

# With coverage
uv run pytest --cov=pywombat --cov-report=html

# MNV tests (requires optional deps)
uv pip install -e ".[dev,mnv]"
uv run pytest tests/test_mnv.py
```

**Coverage expectations**: >80% for new code

## Release Process

### GitHub Workflow

**.github/workflows/publish.yml**:
- Triggered on GitHub release publication or manual trigger
- Builds package using `hatchling`
- Publishes to PyPI using trusted publishing (OIDC token)
- No credentials stored

### Version Management

- **Location**: pyproject.toml (current: 1.4.0)
- **Changelog**: CHANGELOG.md (Keep a Changelog format)
- **Versioning**: Semantic versioning (MAJOR.MINOR.PATCH)

**Recent releases**:
- v1.4.0: MNV (multi-nucleotide variant) detection
- v1.3.4: homalt_only filter bug fix
- v1.3.3: Expression parser support for dotted column names
- v1.2.0: Per-chromosome DNM processing

## Special Considerations

### Memory Optimization

**For expression filtering** (Parquet input):
- Filters applied BEFORE melting
- Reduces memory by 95%+ (200GB → 1.2GB)
- Example: 38 samples, 4.2M variants → < 1 second

**For DNM filtering** (Parquet input):
- Per-chromosome processing
- Reduces memory by 88% (200GB → ~24GB)
- Example: 38 samples, 4.2M variants → 20 seconds

### Lazy Evaluation

**Key principle**: Push filters down to Polars before `.collect()`

```python
# Good pattern
lazy_df = pl.scan_parquet(file)
lazy_df = lazy_df.filter(condition)  # Lazy
result = lazy_df.collect()           # Materialize once

# Avoid
df = pl.scan_parquet(file).collect()  # Load all
df = df.filter(condition)             # Filter after load
```

### Error Handling

**Pattern used throughout**:
```python
try:
    # Main operation
except Exception as e:
    click.echo(f"Error: {e}", err=True)
    raise click.Abort()
```

**Verbose mode**: Print progress to stderr
```python
if verbose:
    click.echo(f"Message", err=True)
```

## Common Modification Patterns

### Adding a New Quality Filter

1. **Add parameter to YAML config**:
```yaml
quality:
  new_filter_threshold: 50
```

2. **Add filter logic** in `apply_quality_filters()` (cli.py, ~line 772):
```python
if "new_filter_threshold" in quality_config:
    threshold = quality_config["new_filter_threshold"]
    df = df.filter(pl.col("column_name") >= threshold)
```

3. **Add test case** in `tests/test_filter.py`:
```python
def test_new_filter():
    config = {"quality": {"new_filter_threshold": 50}}
    # Test implementation
```

4. **Update documentation** in README.md

### Adding a New Expression Operator

1. **Add parsing** in `parse_impact_filter_expression()` (cli.py, ~line 1327):
```python
# Add to parse_condition() nested function
operator_match = re.match(r"^([\w.]+)\s+new_operator\s+(.+)$", condition)
if operator_match:
    col_name = operator_match.group(1)
    value = operator_match.group(2)
    return pl.col(col_name).custom_operation(value)
```

2. **Add test cases** in `tests/test_expression_parser.py`:
```python
def test_new_operator():
    df = pl.DataFrame({"col": [1, 2, 3]})
    expr = parse_impact_filter_expression("col new_operator value", df)
    # Verify behavior
```

3. **Update README.md** with operator documentation and examples

### Debugging Variants

Use `--debug chr:pos` flag to inspect specific variants:

```bash
wombat filter input.tsv -F config.yml --debug chr11:70486013
```

**Shows**:
- All rows matching the position
- VEP_SYMBOL if available
- All columns mentioned in filter expression
- Useful for understanding why variants pass/fail

## Troubleshooting

### Common Errors

#### 1. Import Errors for MNV Dependencies

**Error**:
```
ImportError: No module named 'cyvcf2'
```

**Solution**: Install optional MNV dependencies
```bash
# For end users (published package)
uv pip install pywombat[mnv]

# For development (from source)
uv pip install -e ".[mnv]"
```

#### 2. MNV Detection Errors

**Error**:
```
FileNotFoundError: Cannot open BCF file
UsageError: MNV detection requires --bcf option
```

**Solution**: Ensure BCF is tabix-indexed and FASTA has .fai index
```bash
bcftools index cohort.bcf
samtools faidx genome.fa
```

#### 3. Memory Errors

**Error**:
```
MemoryError: Unable to allocate array
```

**Solutions**:
- Use `wombat prepare` first to create Parquet
- Enable DNM mode for per-chromosome processing
- Apply expression filters (applied before melting)
- Reduce chunk_size in prepare command

#### 4. Expression Parsing Errors

**Error**:
```
ValueError: Column 'field' not found in dataframe
```

**Solution**:
- Check column names in input file
- Ensure VEP annotations are present
- Column names are case-sensitive
- Use `--debug` to inspect available columns

#### 5. Test Failures After Code Changes

**Error**:
```
AssertionError: Expected X rows, got Y
```

**Common causes**:
- Changed filtering logic affects row counts
- Deprecated Polars methods (check warnings)
- Test expectations need updating

### Performance Issues

#### Slow Filtering on Large TSV Files

**Diagnosis**: Direct TSV filtering can be slow

**Solution**: Use two-step workflow
```bash
# Step 1: Prepare once
uv run wombat prepare input.tsv -o prepared.parquet

# Step 2: Filter many times (fast)
uv run wombat filter prepared.parquet -F config.yml -o output
```

#### High Memory Usage

**Diagnosis**: Loading all data into memory

**Solutions**:
1. Use Parquet input (enables lazy evaluation)
2. Enable per-chromosome DNM processing (88% memory reduction)
3. Apply expression filters early (95%+ reduction)

#### Filter Expression Not Reducing Dataset

**Diagnosis**: Filter not applied or syntax error

**Debug approach**:
1. Test expression on small subset
2. Use `--debug chr:pos` to inspect specific variants
3. Check verbose output (`-v`) for filter statistics
4. Verify column names match (case-sensitive)

### Data Format Issues

#### VCF → TSV Conversion Problems

**Issue**: bcftools query format incorrect

**Required format**: Use bcftools +split-vep plugin

**Example**:
```bash
bcftools +split-vep input.vcf.gz \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO\t[%GT:%DP:%GQ:%AD\t]\n' \
  -d -A tab -o output.tsv
```

#### Missing Columns in Output

**Cause**: INFO fields not in input

**Check**: Verify VEP annotation ran with correct flags

**Fix**: Re-annotate VCF with required fields

#### Genotype Parsing Failures

**Issue**: Sample columns not detected

**Requirement**: Sample columns must have format `SampleName:GT:DP:GQ:AD`

**Check**: Column names in input file header (first line)

### Dependency Issues

#### cyvcf2 Installation Fails

**Cause**: Requires htslib C library

**Solution**: Install system dependencies first

```bash
# macOS
brew install htslib

# Ubuntu/Debian
sudo apt-get install libhtslib-dev

# Then install cyvcf2 via UV
uv pip install cyvcf2
```

#### Polars Version Incompatibility

**Symptom**: DeprecationWarning or AttributeError

**Check**: `polars>=0.19.0` required

**Update**:
```bash
uv pip install -U polars
```

#### PyArrow Errors with Parquet

**Error**: `pyarrow.lib.ArrowInvalid`

**Solution**: Ensure `pyarrow>=14.0.0` installed explicitly

```bash
uv pip install pyarrow>=14.0.0
```

### Testing Issues

#### Tests Skipped (MNV Tests)

**Cause**: Optional dependencies not installed

**Fix**: Install full test suite dependencies
```bash
uv pip install -e ".[dev,mnv]"
```

#### Fixture Errors

**Cause**: Test data files missing or wrong format

**Check**: `tests/conftest.py` for fixture definitions

**Regenerate**: Delete tmp files and re-run tests

#### Coverage Too Low

**Target**: >80% for new code

**Check coverage**:
```bash
uv run pytest --cov=pywombat --cov-report=html
```

**View report**: Open `htmlcov/index.html` in browser

## Performance Tips

### For Large Files (>1GB or >50 samples)

1. **Use two-step workflow**:
   ```bash
   uv run wombat prepare large.tsv -o prepared.parquet  # Once
   uv run wombat filter prepared.parquet -F config.yml -o output  # Many times
   ```

2. **For DNM workflows**: Parquet input enables per-chromosome processing automatically (88% memory reduction)

3. **For expression filtering**: Parquet input applies filters before melting (95%+ memory reduction)

4. **Column selection**: Polars can select columns lazily (avoid loading unused columns)

### Memory Optimization Strategies

1. **Reduce chunk size** in prepare command if running out of memory:
   ```bash
   uv run wombat prepare large.tsv -o output.parquet --chunk-size 25000
   ```

2. **Enable streaming** for very large files (already used by default in lazy evaluation)

3. **Process by chromosome** for DNM analysis (automatic with Parquet + DNM config)

---

## Important File Locations

| File | Purpose | Key Functions/Content |
|------|---------|----------------------|
| `src/pywombat/cli.py` | Main CLI | `filter_cmd`, `prepare_cmd`, `vep_cmd`, `apply_quality_filters`, `apply_de_novo_filter`, `parse_impact_filter_expression` |
| `src/pywombat/mnv.py` | MNV detection | `calculate_phasing_probability`, `calculate_indel_inframe`, `detect_mnv_candidates` |
| `src/pywombat/vep.py` | VEP annotation | `parse_variant`, `query_vep`, `extract_annotations`, `annotate_mnv_variants` |
| `pyproject.toml` | Project config | Version, dependencies, build config |
| `examples/rare_homalt.yml` | Recessive analysis | homalt_only filter |
| `examples/de_novo_mutations.yml` | DNM config | DNM settings with PAR regions |
| `examples/rare_homalt_mnv.yml` | Homozygous + MNV | Combined filtering + MNV VEP annotation |

---

This guide should provide AI assistants with the context needed to understand, modify, and extend the pywombat codebase effectively.
