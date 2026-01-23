# Changelog

All notable changes to PyWombat will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
