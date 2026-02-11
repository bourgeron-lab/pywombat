# Changelog

All notable changes to PyWombat will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.1] - 2026-02-11

### Added

- **MNV-only output filter** (`mnv.only: true`): New YAML option to keep only MNV candidates in the output
  - Filters to variants with a non-null `mnv_proba` value (both SNV and indel MNV candidates)
  - Applied after MNV detection and optional VEP annotation
  - Useful for focused MNV analysis without non-MNV variants

### Testing

- Added 4 new tests for `mnv.only` filter in `tests/test_filter.py`

## [1.5.0] - 2026-02-11

### Added

- **VEP Annotation CLI** (`wombat vep`): New command to query the Ensembl VEP REST API for individual variants
  - Accepts variant format `chr:pos:ref:alt` (e.g., `chr17:43094464:G:A`)
  - Annotates on GRCh38 with LOFTEE support (LoF, LoF_filter, LoF_flags, LoF_info)
  - Shows canonical transcript by default; use `--all` for all transcripts
  - Displays gene symbol, consequence, impact, SIFT, PolyPhen, and LOFTEE annotations

- **MNV VEP Annotation** (`mnv.annotate: true`): Automatic re-annotation of MNV candidates via VEP
  - When enabled in YAML config, queries VEP for SNV-based MNV candidates (`mnv_variant` column)
  - Matches VEP response to the same transcript (`VEP_Feature`) as the original variant
  - Appends new annotated rows to the output with MNV variant coordinates
  - Batches API calls (up to 200 per request) with configurable timeout
  - Non-fatal error handling: API failures produce warnings, not aborts

- **New VEP Module** (`src/pywombat/vep.py`):
  - `parse_variant()` - Parses `chr:pos:ref:alt` variant strings
  - `query_vep()` - Queries VEP REST API with batching support
  - `extract_annotations()` - Extracts transcript annotations from VEP responses
  - `format_annotations()` - Formats annotations for CLI display
  - `annotate_mnv_variants()` - Annotates MNV candidates in DataFrames
  - `VEP_COLUMN_MAPPING` - Maps VEP API keys to DataFrame column names

### Changed

- **Renamed `mnv_candidate` to `mnv_proba`**: The MNV phasing probability column has been renamed for clarity
  - Updated across all source files, tests, examples, and documentation
  - The column now consistently represents the phasing probability (0.0-1.0)

### Testing

- Added 61 new tests in `tests/test_vep.py`:
  - `TestParseVariant` (12 tests) - Variant string parsing
  - `TestFormatVariantForVep` (3 tests) - VEP API format conversion
  - `TestExtractAnnotations` (8 tests) - VEP response extraction
  - `TestFormatAnnotations` (6 tests) - CLI display formatting
  - `TestQueryVep` (5 tests) - API querying with mocked responses
  - `TestVepCLICommand` (6 tests) - CLI integration
  - `TestVepColumnMapping` (2 tests) - Column mapping validation
  - `TestFindTranscriptAnnotation` (4 tests) - Transcript matching
  - `TestAnnotateMnvVariants` (15 tests) - MNV annotation pipeline

## [1.2.1] - 2026-02-05

### Fixed

- **Missing Dependency**: Added `pyarrow>=14.0.0` as an explicit dependency
  - Required for Parquet file operations (`scan_parquet`, `write_parquet`)
  - Previously was an implicit dependency through Polars
  - Ensures proper installation on all systems

## [1.2.0] - 2026-02-05

### Added

- **Per-Chromosome DNM Processing**: Dramatically reduced memory usage for de novo mutation (DNM) filtering
  - Processes one chromosome at a time instead of loading all variants into memory
  - Reduces peak memory from (total_variants × samples) to (max_chr_variants × samples)
  - Example: 38 samples, 4.2M variants
    - Before: 200GB+ (OOM failure)
    - After: ~24GB (completes successfully in 20 seconds)
  - **88% memory reduction** for DNM workflows

- **Early Frequency Filtering for DNM**: Applies population frequency filters BEFORE melting
  - Frequency filters (fafmax_faf95_max_genomes) applied on wide-format data
  - Quality filters (genomes_filters PASS) applied before melting
  - Reduces data expansion by filtering variants early in the pipeline

- **New Helper Functions**:
  - `get_unique_chromosomes()`: Discovers and naturally sorts chromosomes from Parquet files
  - `apply_dnm_prefilters()`: Applies variant-level filters before melting
  - `process_dnm_by_chromosome()`: Orchestrates per-chromosome DNM filtering

### Changed

- **DNM Filter Architecture**: Refactored `apply_de_novo_filter()` to support `skip_prefilters` parameter
  - Allows separation of variant-level filters (applied before melting) from sample-level filters
  - Prevents double-filtering when prefilters already applied

- **Filter Command Routing**: Automatically detects DNM mode and routes to per-chromosome processing
  - Transparent to users - no command syntax changes required
  - Optimized memory usage is automatic when using DNM config with Parquet input

### Performance

- **DNM Memory Usage**: 88% reduction in peak memory (200GB+ → ~24GB)
- **DNM Processing Time**: 20 seconds for 38-sample, 4.2M variant dataset (previously failed with OOM)
- **Throughput**: Successfully processes 6,788 DNM variants from 4.2M input variants

### Testing

- Added 3 new test cases for DNM optimization:
  - `test_get_unique_chromosomes()`: Verifies chromosome discovery and natural sorting
  - `test_apply_dnm_prefilters()`: Validates frequency prefiltering logic
  - `test_dnm_skip_prefilters()`: Ensures skip_prefilters parameter works correctly
- Total test suite: 25 tests (all passing)

## [1.1.0] - 2026-02-05

### Added

- **Memory-Optimized Two-Step Workflow**: New `wombat prepare` command for preprocessing large files
  - Converts TSV/TSV.gz to Parquet format with pre-expanded INFO fields
  - Processes files in chunks (50k rows default) to handle files of any size
  - Applies memory-efficient data types (Categorical, UInt32, etc.)
  - Reduces file size by ~30% compared to gzipped TSV

- **Parquet Input Support**: `wombat filter` now accepts both TSV and Parquet input
  - Auto-detects input format (TSV, TSV.gz, or Parquet)
  - Pre-filtering optimization: Applies expression filters BEFORE melting samples
  - Reduces memory usage by 95%+ for large files (e.g., 200GB → 1.2GB for 38-sample, 4.2M variant dataset)
  - Processing time improved from minutes/OOM to <1 second for filtered datasets

- **Subcommand Architecture**: Converted CLI to use Click groups
  - `wombat prepare`: Preprocess TSV to Parquet
  - `wombat filter`: Process and filter data (replaces old direct command)
  - **Breaking Change**: Old syntax `wombat input.tsv` no longer works, use `wombat filter input.tsv`

- **Test Suite**: Added comprehensive pytest test suite
  - 22 tests covering CLI structure, prepare command, and filter command
  - Test fixtures for creating synthetic test data
  - Integration tests with real data validation
  - Added pytest and pytest-cov to dev dependencies

### Changed

- **CLI Architecture**: Restructured from single command to group-based subcommands
- **Filter Command**: Now applies expression filters before melting when using Parquet input
- **Sample Column Detection**: Improved heuristics to work with both TSV and Parquet formats
- **Documentation**: Updated README with two-step workflow examples and memory comparison table

### Fixed

- **INFO Field Extraction**: Fixed column index detection in `prepare` command (was using hardcoded index)
- **Type Casting**: Added explicit `.cast(pl.Utf8)` to preserve string types when all values are NULL
- **Parquet Processing**: Fixed `format_bcftools_tsv_minimal` to work without `(null)` column

### Performance

- **Memory Usage**: 95%+ reduction for large files with expression filters
  - Example: 38 samples, 4.2M variants
  - Before: 200GB+ (OOM failure)
  - After: ~1.2GB peak memory
- **Processing Speed**: <1 second for filtered datasets (vs minutes or failure before)
- **Pre-filtering**: Expression filters applied before melting reduces data expansion

### Documentation

- Added memory optimization workflow section to README
- Added performance comparison table showing memory/time improvements
- Updated all examples to use new `wombat filter` syntax
- Added section explaining when to use `prepare` command
- Documented two-step workflow benefits and use cases

## [1.0.1] - 2026-01-24

### Added

- **Boolean Flag Support in INFO Fields**: INFO field entries without `=` signs (e.g., `PASS`, `DB`, `SOMATIC`) are now extracted as boolean columns with `True`/`False` values instead of being ignored. This enables filtering on VCF flag fields.

### Fixed

- INFO field parsing now handles all field types correctly, including standalone boolean flags commonly used in VCF files.

## [1.0.0] - 2026-01-23

First stable release of PyWombat! 🎉

### Added

#### Core Features

- **Fast TSV Processing**: Efficient processing of bcftools tabulated TSV files using Polars
- **Flexible Output Formats**: Support for TSV, compressed TSV (`.gz`), and Parquet formats
- **Streaming Mode**: Memory-efficient processing for large files
- **Pedigree Support**: Trio and family analysis with automatic parent genotype joining
- **Multiple Sample Formats**: Handles various genotype formats (GT:DP:GQ:AD and variants)

#### Filtering Capabilities

- **Quality Filters**: Configurable thresholds for depth (DP), genotype quality (GQ), and variant allele frequency (VAF)
- **Genotype-Specific VAF Filters**: Separate thresholds for heterozygous, homozygous alternate, and homozygous reference calls
- **Expression-Based Filtering**: Complex logical expressions with comparison operators (`=`, `!=`, `<`, `>`, `<=`, `>=`) and logical operators (`&`, `|`)
- **Parent Quality Filtering**: Optional quality filter application to parent genotypes

#### De Novo Mutation Detection

- **Sex-Chromosome Aware Logic**: Proper handling of X and Y chromosomes in males
- **PAR Region Support**: Configurable pseudo-autosomal region (PAR) coordinates for GRCh37 and GRCh38
- **Hemizygous Variant Detection**: Specialized VAF thresholds for X chromosome in males (non-PAR) and Y chromosome
- **Homozygous VAF Thresholds**: Higher VAF requirements (≥85%) for homozygous variants
- **Parent Genotype Validation**: Ensures parents are homozygous reference with low VAF (<2%)
- **Missing Genotype Filtering**: Removes variants with partial/missing genotypes (`./.`, `0/.`, etc.)
- **Population Frequency Filtering**: Maximum allele frequency thresholds (gnomAD fafmax_faf95_max_genomes)
- **Quality Filter Support**: gnomAD genomes_filters PASS-only option

#### User Experience

- **Debug Mode**: Inspect specific variants by chromosome:position for troubleshooting
- **Verbose Mode**: Detailed filtering step information with variant counts
- **Automatic Output Naming**: Intelligent output file naming based on input and filter config
- **Configuration Examples**: Two comprehensive example configurations with extensive documentation
  - `rare_variants_high_impact.yml`: Ultra-rare, high-impact variant filtering
  - `de_novo_mutations.yml`: De novo mutation detection with full documentation

#### Documentation

- **Comprehensive README**: Complete usage guide with examples for all features
- **Example Workflows**: Real-world usage scenarios (rare disease, autism trios, etc.)
- **Input Requirements**: Detailed bcftools command examples for generating input files
- **VEP Annotation Guide**: Complete workflow from VEP annotation to PyWombat processing
- **Examples Directory**: Dedicated directory with configuration files and detailed README
- **Troubleshooting Section**: Common issues and solutions

#### Installation Methods

- **uvx Support**: One-line execution without installation (`uvx pywombat`)
- **uv Development Mode**: Local installation for repeated use (`uv sync`, `uv run wombat`)

### Changed

- Improved performance with streaming lazy operations
- Optimized parent genotype lookup (excludes 0/0 genotypes from storage)
- Enhanced error messages for better user experience
- Normalized chromosome names for PAR region matching (handles both 'X' and 'chrX')

### Fixed

- Sex column reading from pedigree file
- Parent genotype column naming consistency (father_id/mother_id)
- Genotype filtering to catch all partial genotypes (`./.`, `0/.`, `1/.`)
- PAR region matching for different chromosome naming conventions
- Empty chunk handling in output to avoid blank lines

### Performance Optimizations

- Delayed annotation expansion (filter before expanding `(null)` field)
- Vectorized filtering operations (no Python loops)
- Early genotype filtering (skip 0/0 before parent lookup)
- Optimized parent lookup (stores only non-reference genotypes)
- Streaming mode by default for memory efficiency

### Removed

- **Progress bar options**: Removed `--progress`/`--no-progress` and `--chunk-size` options for simplicity
- **Chunked processing mode**: Simplified to use only efficient streaming mode

## [0.5.0] - 2026-01-20

### Added

- Initial de novo mutation detection implementation
- Pedigree file support
- Basic quality filtering
- Expression-based filtering

### Known Issues

- Progress bar had reliability issues (removed in 1.0.0)
- Chunked processing was complex (simplified in 1.0.0)

---

## Release Notes

### v1.0.0 - Production Ready

This release marks PyWombat as production-ready for:

- Rare disease gene discovery
- De novo mutation detection in autism and developmental disorders
- Trio and family-based variant analysis
- High-throughput variant filtering workflows

**Recommended for**: Research groups working with rare variants, de novo mutations, and family-based genomic studies.

**Breaking Changes**: None from 0.5.0, but removed progress bar options for cleaner interface.

---

[1.0.0]: https://github.com/bourgeron-lab/pywombat/releases/tag/v1.0.0
[0.5.0]: https://github.com/bourgeron-lab/pywombat/releases/tag/v0.5.0
