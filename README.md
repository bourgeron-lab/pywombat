# PyWombat 🦘

A high-performance CLI tool for processing and filtering bcftools tabulated TSV files with advanced filtering capabilities, pedigree support, and de novo mutation detection.

[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- ✨ **Fast Processing**: Uses Polars for efficient data handling
- 🔬 **Quality Filtering**: Configurable depth, quality, and VAF thresholds
- 👨‍👩‍👧 **Pedigree Support**: Trio and family analysis with parent genotypes
- 🧬 **De Novo Detection**: Sex-chromosome-aware DNM identification
- 📊 **Flexible Output**: TSV, compressed TSV, or Parquet formats
- 🎯 **Expression Filters**: Complex filtering with logical expressions
- 🏷️ **Boolean Flag Support**: INFO field flags (PASS, DB, etc.) extracted as True/False columns
- ⚡ **Memory Optimized**: Two-step workflow for large files (prepare → filter)
- 💾 **Parquet Support**: Pre-process large files for repeated, memory-efficient analysis
- 🔎 **VEP Annotation**: Query Ensembl VEP REST API for individual variants or MNV candidates
- 🧩 **MNV Annotation**: Automatic re-annotation of MNV candidates via VEP with transcript matching

---

## Quick Start

### One-Time Usage (Recommended for Most Users)

Use `uvx` to run PyWombat without installation:

```bash
# Basic filtering
uvx pywombat filter input.tsv -o output

# With filter configuration
uvx pywombat filter input.tsv -F examples/rare_variants_high_impact.yml -o output

# De novo mutation detection
uvx pywombat filter input.tsv --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml -o denovo
```

### For Large Files (>1GB or >50 samples)

Use the two-step workflow for memory-efficient processing:

```bash
# Step 1: Prepare (one-time preprocessing)
uvx pywombat prepare input.tsv.gz -o prepared.parquet

# Step 2: Filter (fast, memory-efficient, can be run multiple times)
uvx pywombat filter prepared.parquet \
  -p pedigree.tsv \
  -F config.yml \
  -o filtered
```

**Benefits:**
- Pre-expands INFO fields once (saves time on repeated filtering)
- Applies filters before melting samples (reduces memory by 95%+)
- Parquet format enables fast columnar access

### Installation for Development/Repeated Use

```bash
# Clone the repository
git clone https://github.com/bourgeron-lab/pywombat.git
cd pywombat

# Install with uv
uv sync

# Run with uv run
uv run wombat filter input.tsv -o output
```

---

## What Does PyWombat Do?

PyWombat transforms bcftools tabulated TSV files into analysis-ready formats by:

1. **Expanding the `(null)` INFO column**: Extracts all `NAME=value` fields (e.g., `DP=30;AF=0.5;AC=2`) and boolean flags (e.g., `PASS`, `DB`) into separate columns
2. **Melting sample columns**: Converts wide-format sample data into long format with one row per variant-sample combination
3. **Extracting genotype data**: Parses `GT:DP:GQ:AD` format into separate columns with calculated VAF
4. **Adding parent data**: Joins father/mother genotypes when pedigree is provided
5. **Applying filters**: Quality filters, expression-based filters, or de novo detection

### Example Transformation

**Input (Wide Format):**

```tsv
#CHROM  POS  REF  ALT  (null)                      Sample1:GT:DP:GQ:AD  Sample2:GT:DP:GQ:AD
chr1    100  A    T    DP=30;AF=0.5;PASS;AC=2      0/1:15:99:5,10      1/1:18:99:0,18
```

**Output (Long Format):**

```tsv
#CHROM  POS  REF  ALT  AC  AF   DP  PASS  sample   sample_gt  sample_dp  sample_gq  sample_ad  sample_vaf
chr1    100  A    T    2   0.5  30  true  Sample1  0/1        15         99         10         0.6667
chr1    100  A    T    2   0.5  30  true  Sample2  1/1        18         99         18         1.0
```

**Generated Columns:**

- INFO fields with `=`: Extracted as separate columns (e.g., `DP`, `AF`, `AC`)
- INFO boolean flags: Extracted as True/False columns (e.g., `PASS`, `DB`, `SOMATIC`)
- `sample`: Sample identifier
- `sample_gt`: Genotype (e.g., 0/1, 1/1)
- `sample_dp`: Read depth (total coverage)
- `sample_gq`: Genotype quality score
- `sample_ad`: Alternate allele depth (from second value in AD field)
- `sample_vaf`: Variant allele frequency (sample_ad / sample_dp)

---

## Commands

PyWombat has three main commands:

### `wombat prepare` - Preprocess Large Files

Converts TSV/TSV.gz to optimized Parquet format with pre-expanded INFO fields:

```bash
# Basic usage
wombat prepare input.tsv.gz -o prepared.parquet

# With verbose output
wombat prepare input.tsv.gz -o prepared.parquet -v

# Adjust chunk size for memory constraints
wombat prepare input.tsv.gz -o prepared.parquet --chunk-size 25000
```

**What it does:**
- Extracts all INFO fields (VEP_*, AF, etc.) as separate columns
- Keeps samples in wide format (not melted yet)
- Writes memory-efficient Parquet format
- Processes in chunks to handle files of any size

**When to use:**
- Files >1GB or >50 samples
- Large families (>10 members)
- Running multiple filter configurations
- Repeated analysis of the same dataset

### `wombat filter` - Process and Filter Data

Transforms and filters variant data (works with both TSV and Parquet input):

```bash
# Basic filtering (TSV input)
wombat filter input.tsv -o output

# From prepared Parquet (faster, more memory-efficient)
wombat filter prepared.parquet -o output

# With filter configuration
wombat filter input.tsv -F config.yml -o output

# With pedigree
wombat filter input.tsv -p pedigree.tsv -o output

# Compressed output
wombat filter input.tsv -o output -f tsv.gz

# Parquet output
wombat filter input.tsv -o output -f parquet

# With verbose output
wombat filter input.tsv -o output -v
```

### `wombat vep` - Annotate a Variant via Ensembl VEP

Query the Ensembl VEP REST API for a single variant on GRCh38:

```bash
# Annotate a variant (shows canonical transcript only)
wombat vep chr1:151776102:GC:AA

# Show all transcripts
wombat vep chr1:151776102:GC:AA --all

# With verbose output
wombat vep chr17:43094464:G:A -v
```

**What it does:**
- Queries the Ensembl VEP REST API (GRCh38) with LOFTEE annotations
- Displays transcript-level annotations: gene symbol, consequence, impact, SIFT, PolyPhen, LOFTEE
- By default shows only the canonical transcript; use `--all` for all transcripts
- Supports SNVs, indels, and MNVs in `chr:pos:ref:alt` format

**Variant format**: `chrN:POS:REF:ALT` (e.g., `chr1:151776102:GC:AA`, `chr17:43094464:G:A`)

### With Pedigree (Trio/Family Analysis)

Add parent genotype information for inheritance analysis:

```bash
wombat filter input.tsv --pedigree pedigree.tsv -o output
```

**Pedigree File Format** (tab-separated):

```tsv
FID sample_id FatherBarcode MotherBarcode Sex Pheno
FAM1 Child1 Father1 Mother1 1 2
FAM1 Father1 0 0 1 1
FAM1 Mother1 0 0 2 1
```

- `FID`: Family ID
- `sample_id`: Sample name (must match VCF)
- `FatherBarcode`: Father's sample name (0 = unknown)
- `MotherBarcode`: Mother's sample name (0 = unknown)
- `Sex`: 1=male, 2=female (or M/F)
- `Pheno`: 1=unaffected, 2=affected

**Output with pedigree includes additional columns:**

- `father_gt`, `father_dp`, `father_gq`, `father_ad`, `father_vaf`
- `mother_gt`, `mother_dp`, `mother_gq`, `mother_ad`, `mother_vaf`

---

## Advanced Filtering

PyWombat supports two types of filtering:

1. **Expression-based filtering**: For rare variants, impact filtering, frequency filtering
2. **De novo mutation detection**: Specialized logic for identifying DNMs in trios/families

### 1. Rare Variant Filtering

Filter for ultra-rare, high-impact variants:

```bash
wombat filter input.tsv \
  -F examples/rare_variants_high_impact.yml \
  -o rare_variants
```

**Configuration** (`rare_variants_high_impact.yml`):

```yaml
quality:
  sample_dp_min: 10           # Minimum read depth
  sample_gq_min: 19           # Minimum genotype quality
  sample_vaf_het_min: 0.25    # Het VAF range: 25-75%
  sample_vaf_het_max: 0.75
  sample_vaf_homalt_min: 0.85 # Hom alt VAF ≥ 85%
  apply_to_parents: true      # Also filter parent genotypes
  filter_no_alt_allele: true  # Exclude 0/0 genotypes

expression: "VEP_CANONICAL = YES & VEP_IMPACT = HIGH & VEP_LoF = HC & VEP_LoF_flags = . & ( fafmax_faf95_max_genomes = null | fafmax_faf95_max_genomes <= 0.001 )"
```

**What this does:**

- Filters for canonical transcripts with HIGH impact
- Loss-of-function (LoF) with high confidence, no flags
- Ultra-rare: ≤0.1% frequency in gnomAD genomes
- Stringent quality: DP≥10, GQ≥19, appropriate VAF

### 2. De Novo Mutation Detection

Identify de novo mutations in trio data:

```bash
wombat filter input.tsv \
  --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml \
  -o denovo
```

**Configuration** (`de_novo_mutations.yml`):

```yaml
quality:
  sample_dp_min: 10
  sample_gq_min: 18
  sample_vaf_min: 0.20

dnm:
  enabled: true
  parent_dp_min: 10          # Parent quality thresholds
  parent_gq_min: 18
  parent_vaf_max: 0.02       # Max 2% VAF in parents
  
  sample_vaf_hemizygous_min: 0.85  # For X male, Y, and hom variants
  
  fafmax_faf95_max_genomes_max: 0.001  # Max 0.1% in gnomAD
  genomes_filters_pass_only: true      # Only PASS variants
  
  par_regions:               # Pseudo-autosomal regions (GRCh38)
    grch38:
      PAR1:
        chrom: X
        start: 10000
        end: 2781479
      PAR2:
        chrom: X
        start: 155701383
        end: 156030895
```

**Sex-Chromosome Aware Logic:**

- **Autosomes & PAR**: Both parents must be 0/0 with VAF<2%
- **X in males (non-PAR)**: Mother must be 0/0, father not informative
- **Y in males**: Father must be 0/0, mother N/A
- **Hemizygous variants**: Require VAF ≥ 85%

### 3. Recessive Disease Analysis (Homozygous Alternative Only)

Filter for only homozygous alternative (1/1) variants in recessive disease studies:

```bash
wombat filter input.tsv \
  -F examples/rare_homalt.yml \
  -o recessive_variants

# With pedigree for trio analysis
wombat filter input.tsv \
  --pedigree pedigree.tsv \
  -F examples/rare_homalt.yml \
  -o recessive_variants
```

**Configuration** (`rare_homalt.yml`):

```yaml
quality:
  sample_dp_min: 10
  sample_gq_min: 19
  sample_vaf_homalt_min: 0.98  # Stringent VAF for true hom calls

  homalt_only: true            # Only keep 1/1 genotypes
  filter_no_alt_allele: true

expression: "VEP_CANONICAL = YES & VEP_IMPACT = HIGH & (fafmax_faf95_max_genomes = null | fafmax_faf95_max_genomes <= 0.001)"
```

**What this does:**

- Filters to **only homozygous alternative (1/1)** genotypes
- Excludes heterozygous (0/1), multi-allelic (2/2, 3/3), and hom ref (0/0)
- Parents are NOT filtered (they remain heterozygous carriers in recessive inheritance)
- Ultra-rare frequency (≤0.1% in gnomAD genomes)
- Canonical transcripts with HIGH impact

**Use cases:**
- Autosomal recessive disease analysis
- Homozygosity mapping in consanguineous families
- Identifying runs of homozygosity (ROH)
- Compound heterozygote analysis (first pass for 1/1, second for 0/1)

---

## Expression-Based Filtering

Create custom filters using logical expressions:

### Available Operators

- **Comparison**: `=`, `!=`, `<`, `>`, `<=`, `>=`
- **Logical**: `&` (AND), `|` (OR)
- **Grouping**: `(`, `)`
- **Null checks**: `= null`, `!= null`
- **String matching**: `contains` (case-insensitive substring search)
- **Empty checks**: `is_empty` (checks for null, empty string, or "."), `not_empty` (negation of is_empty)
- **Variant type**: `is_snv` (both REF and ALT are single nucleotides A/C/G/T)
- **Variant type**: `is_indel` (insertions, deletions, MNVs - not SNVs)

**Note**: Column names can contain dots (e.g., `cadd_v1.7`, `faf95.joint`) and alphanumeric characters with underscores.

### Example Expressions

```yaml
# High or moderate impact with low frequency
expression: "(VEP_IMPACT = HIGH | VEP_IMPACT = MODERATE) & gnomad_AF < 0.001"

# Canonical transcripts with CADD score
expression: "VEP_CANONICAL = YES & CADD_PHRED >= 20"

# Specific consequence types
expression: "VEP_Consequence = frameshift_variant | VEP_Consequence = stop_gained"

# Multiple criteria
expression: "VEP_IMPACT = HIGH & VEP_CANONICAL = YES & gnomad_AF < 0.01 & CADD_PHRED >= 25"

# Pathogenic variants that passed gnomAD filters
expression: "VEP_CLIN_SIG contains 'pathogenic' & genomes_filters is_empty"

# Pathogenic or likely pathogenic variants
expression: "(VEP_CLIN_SIG contains 'pathogenic' | VEP_CLIN_SIG contains 'likely_pathogenic') & genomes_filters is_empty"

# Different quality thresholds for SNVs vs INDELs
expression: "(is_snv & QUAL > 18) | (is_indel & QUAL > 22)"

# High-quality SNVs only
expression: "is_snv & QUAL > 30 & genomes_filters is_empty"

# High-impact INDELs with pathogenic annotation
expression: "is_indel & VEP_IMPACT = HIGH & VEP_CLIN_SIG contains 'pathogenic'"

# Filter for variants with ClinVar annotations present
expression: "VEP_CLIN_SIG not_empty & VEP_CANONICAL = YES"

# Require both ClinVar and LoF annotations
expression: "VEP_CLIN_SIG not_empty & VEP_LoF not_empty & VEP_IMPACT = HIGH"

# Combine is_empty and not_empty for complex logic
# Variants with annotations that passed gnomAD filters
expression: "VEP_CLIN_SIG not_empty & genomes_filters is_empty"

# Column names with dots are supported
expression: "cadd_v1.7 >= 28.3 & REVEL not_empty"
```

---

## MNV Detection (Multi-Nucleotide Variants)

PyWombat can detect MNV candidates - cases where multiple nearby variants on the same haplotype may combine to create a different functional effect than the individual variants alone.

### What are MNVs?

Multi-nucleotide variants occur when two or more nearby variants are inherited together on the same chromosome (in cis). These compound variants can have significantly different functional consequences:

- **SNV pairs**: Two SNVs within 2bp may form a di-nucleotide substitution with different protein impact
- **Indel clusters**: Multiple indels that may restore or maintain reading frame
- **Functional impact**: Combined variants may be synonymous when individual variants are missense, or vice versa

### Use Cases

1. **Di-nucleotide Substitutions**:
   - Individual: chr1:100:A:T (missense) + chr1:101:G:C (missense)
   - Combined MNV: chr1:100:AG:TC (may be synonymous!)
   - Impact: Completely different functional prediction

2. **Separated SNV Pairs**:
   - Individual: chr1:200:C:T (missense) + chr1:202:G:A (missense)
   - Combined MNV: chr1:200:CCG:TCA (different codon change)
   - Impact: May affect protein structure differently

3. **Frameshift Restoration** (or lack thereof):
   - Indel1: chr2:500:ATG:A (-2bp frameshift)
   - Indel2: chr2:505:C:CGAT (+3bp frameshift)
   - Combined: -2+3 = +1bp (still frameshift, not rescued)
   - Impact: Understanding combined effect

### Configuration

Enable MNV detection by adding the `mnv` section to your filter config:

```yaml
mnv:
  candidate: true      # Enable MNV candidate detection
  indel_window: 10     # Window size for indel clustering (default: 10)
  annotate: true       # Re-annotate MNV candidates via Ensembl VEP REST API
  only: false          # Set to true to keep only MNV candidates in output
  # annotate_timeout: 30      # HTTP timeout for VEP API requests (default: 30s)
  # annotate_batch_size: 200  # Max variants per VEP API batch (default/max: 200)
  exclude_mhc:                # Exclude MHC region (noisy, unreliable calling)
    enabled: true
    chrom: "6"                # Matches both "6" and "chr6"
    start: 28510120           # GRCh38 xMHC start
    end: 33480577             # GRCh38 xMHC end
```

### Required Command-Line Options

MNV detection requires two additional files:

```bash
wombat filter variants.parquet \
  -F config_with_mnv.yml \
  --bcf cohort.bcf \      # Tabix-indexed BCF file (same as input source)
  --fasta GRCh38.fa \     # Indexed FASTA reference genome
  -o output.tsv
```

**Files needed**:
- `--bcf`: Tabix-indexed BCF file (`.bcf.csi` or `.bcf.tbi` index must exist)
- `--fasta`: Indexed FASTA reference (`.fai` index must exist)

### Installation

MNV detection requires additional dependencies:

```bash
# Install with MNV support
pip install pywombat[mnv]

# Or with uv
uv pip install pywombat[mnv]
```

**Dependencies added**: `cyvcf2`, `pysam`, `scipy`

### Output Columns

Three new columns are added to the output:

| Column | Type | Description |
|--------|------|-------------|
| `mnv_proba` | float or null | Phasing probability (0.0-1.0) indicating confidence that variants are in-phase. Higher values = stronger evidence. |
| `mnv_inframe` | bool or null | For indels only: `True` if combined indels maintain reading frame (net length is multiple of 3). |
| `mnv_variant` | string | For SNVs only: Reconstructed variant notation (chr:pos:ref:alt) for re-annotation. |

### Example Output

```
#CHROM  POS  REF  ALT  sample  mnv_proba  mnv_inframe  mnv_variant
chr1    100  A    T    S001    0.95          null          chr1:100:AGC:TGC
chr1    102  G    C    S001    0.95          null          chr1:100:AGC:TGC
chr2    500  ATG  A    S002    0.88          true
chr2    505  C    CGAT S002    0.88          true
chr3    200  A    T    S003    null          null
```

**Interpretation**:
- chr1:100-102 (S001): Two SNVs with 95% probability of being in-phase → likely MNV
- chr2:500-505 (S002): Two indels with 88% probability, combined are inframe
- chr3:200 (S003): No nearby variants found → not an MNV candidate

### Detection Logic

**For SNVs** (single nucleotide variants):
1. Search window: pos-2 to pos+2 (5bp total)
2. Find other SNVs carried by the same sample
3. Calculate phasing probability using allelic depth (AD) values
4. Reconstruct combined variant notation using reference genome
5. Output: `mnv_proba` (probability), `mnv_variant` (notation)

**For INDELs** (insertions/deletions):
1. Calculate indel length: `abs(len(REF) - len(ALT))`
2. Search window: `pos-indel_window` to `pos+length+indel_window`
3. Find other indels carried by the same sample
4. Calculate phasing probability using AD values
5. Sum net lengths to determine if inframe (multiple of 3)
6. Output: `mnv_proba` (probability), `mnv_inframe` (bool)

**Phasing Probability Calculation**:
- Compares allelic depth (AD) patterns between variants
- **In-phase (cis)**: AD values similar (ref1≈ref2, alt1≈alt2)
- **Out-of-phase (trans)**: AD values swapped (ref1≈alt2, alt1≈ref2)
- Uses Gaussian model to calculate `P(in-phase) / P(total)`
- Values close to 1.0 = high confidence in-phase
- Values close to 0.0 = likely out-of-phase

### Example Configuration

See [examples/rare_homalt_mnv.yml](examples/rare_homalt_mnv.yml) for a complete example with:
- Quality filters for homozygous alternative variants
- Expression filters for rare, high-impact variants
- MNV detection enabled with detailed comments
- Usage examples and customization tips

### Automatic MNV Annotation via VEP

When `annotate: true` is set in the MNV config, PyWombat automatically:

1. Identifies SNV-based MNV candidates (rows with non-empty `mnv_variant`)
2. Queries the Ensembl VEP REST API in batches (up to 200 variants per request)
3. Matches the VEP response to the same transcript (`VEP_Feature`) as the original variant
4. Appends new rows to the output with the MNV's VEP annotation

**New rows have:**
- `#CHROM`, `POS`, `REF`, `ALT` from the reconstructed MNV variant
- Sample columns (`sample_gt`, `sample_dp`, etc.) copied from the original row
- Parent columns set to null
- `mnv_proba` set to `1 - original_value` (complement)
- `VEP_*` columns populated from the VEP API response

This allows direct comparison of individual vs. combined variant annotations in a single output file.

**Note**: Requires internet access to reach the Ensembl REST API. Only SNV-based MNVs are annotated (indel MNVs produce `mnv_inframe` but no reconstructed variant).

### Downstream Analysis

**Filter by phasing confidence**:
```bash
# Keep only high-confidence MNV candidates
awk '$mnv_proba > 0.8 || $mnv_proba == ""' output.tsv > high_conf_mnvs.tsv
```

**Identify rescued frameshifts**:
```bash
# Find indel pairs that restore reading frame
awk '$mnv_inframe == "true"' output.tsv > rescued_frameshifts.tsv
```

### Limitations

- Requires original BCF and reference FASTA files
- SNV window is fixed at 2bp (based on typical codon spacing)
- Indel window is configurable (default: 10bp)
- Phasing probability is statistical estimate, not definitive proof
- True phasing validation requires read-level analysis (future: `validate: true`)

### Future Features (Planned)

- `validate: true` - Validate phasing using CRAM/BAM read alignments
- Configurable SNV window size
- Multi-sample MNV detection (family-based)

---

## Debug Mode

Inspect specific variants for troubleshooting:

```bash
wombat filter input.tsv \
  -F config.yml \
  --debug chr11:70486013
```

Shows:

- All rows matching the position
- VEP_SYMBOL if available
- All columns mentioned in filter expression
- Useful for understanding why variants pass/fail filters

---

## Output Formats

### TSV (Default)

```bash
wombat filter input.tsv -o output          # Creates output.tsv
wombat filter input.tsv -o output -f tsv   # Same as above
```

### Compressed TSV

```bash
wombat filter input.tsv -o output -f tsv.gz  # Creates output.tsv.gz
```

### Parquet (Fastest for Large Files)

```bash
wombat filter input.tsv -o output -f parquet  # Creates output.parquet
```

**When to use Parquet:**

- Large cohorts (>100 samples)
- Downstream analysis with Polars/Pandas
- Need for fast repeated access
- Storage efficiency

---

## Example Workflows

### 1. Rare Disease Gene Discovery

```bash
# Step 1: Filter for rare, high-impact variants
wombat filter cohort.tsv \
  -F examples/rare_variants_high_impact.yml \
  -o rare_variants

# Step 2: Further filter with gene list (in downstream analysis)
# Use the output TSV with your favorite tools (R, Python, etc.)
```

### 2. Autism Trio Analysis

```bash
# Identify de novo mutations in autism cohort
wombat filter autism_trios.tsv \
  --pedigree autism_pedigree.tsv \
  -F examples/de_novo_mutations.yml \
  -o autism_denovo \
  -v

# Review output for genes in autism risk lists
```

### 3. Large Multi-Family Analysis (Memory-Optimized)

```bash
# Step 1: Prepare once (preprocesses INFO fields)
wombat prepare large_cohort.tsv.gz -o prepared.parquet -v

# Step 2: Filter with different configurations (fast, memory-efficient)
wombat filter prepared.parquet \
  --pedigree families_pedigree.tsv \
  -F examples/rare_variants_high_impact.yml \
  -o families_rare_variants \
  -v

# Step 3: Run additional filters without re-preparing
wombat filter prepared.parquet \
  --pedigree families_pedigree.tsv \
  -F examples/de_novo_mutations.yml \
  -o families_denovo \
  -v
```

### 4. Custom Expression Filter

Create `custom_filter.yml`:

```yaml
quality:
  sample_dp_min: 15
  sample_gq_min: 20
  sample_vaf_het_min: 0.30
  sample_vaf_het_max: 0.70

expression: "VEP_IMPACT = HIGH & (gnomad_AF < 0.0001 | gnomad_AF = null)"
```

Apply:

```bash
wombat filter input.tsv -F custom_filter.yml -o output
```

---

## Input Requirements

### bcftools Tabulated Format

PyWombat expects TSV files created with bcftools:

```bash
# From VCF to tabulated TSV
bcftools query -HH -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT:%DP:%GQ:%AD]\n' \
  input.vcf.gz > input.tsv

# With specific INFO fields
bcftools query -HH \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP;%INFO/AF;%INFO/AC[\t%GT:%DP:%GQ:%AD]\n' \
  input.vcf.gz > input.tsv
```

**Expected format:**

- Tab-separated values
- Header row with column names (use `-HH` for proper headers)
- `(null)` column containing INFO fields
- Sample columns in `GT:DP:GQ:AD` format (or similar)
- Optional FILTER column for quality control

### VEP Annotations

For expression-based filtering on VEP annotations, a two-step process is required:

#### Step 1: Split VEP CSQ Field

**IMPORTANT**: PyWombat requires VEP annotations to be split into individual fields with the `VEP_` prefix:

```bash
# Split VEP CSQ field - creates one row per transcript/consequence
# This prefixes all CSQ columns with "VEP_" (e.g., VEP_SYMBOL, VEP_IMPACT)
bcftools +split-vep -c - -p VEP_ -O b -o annotated.split.bcf input.vcf.gz
```

#### Step 2: Convert to Tabulated Format

```bash
# Extract to TSV with genotype information
bcftools query -HH \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT:%DP:%GQ:%AD]\n' \
  annotated.split.bcf > annotated.tsv
```

#### Complete VEP Workflow

```bash
# 1. Annotate with VEP
vep -i input.vcf.gz \
  --cache --offline \
  --format vcf \
  --vcf \
  --everything \
  --canonical \
  --plugin LoF \
  -o annotated.vcf.gz

# 2. Split VEP CSQ field (REQUIRED for PyWombat)
bcftools +split-vep -c - -p VEP_ -O b -o annotated.split.bcf annotated.vcf.gz

# 3. Convert to tabulated format
bcftools query -HH \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT:%DP:%GQ:%AD]\n' \
  annotated.split.bcf > annotated.tsv

# 4. Process with PyWombat
wombat filter annotated.tsv -F examples/rare_variants_high_impact.yml -o output
```

**Why split-vep is required:**

- Creates one row per transcript/consequence (instead of all in one CSQ field)
- Prefixes VEP fields with `VEP_` making them accessible in expressions
- Enables filtering on `VEP_CANONICAL`, `VEP_IMPACT`, `VEP_LoF`, etc.

#### Pipeline Optimization

For production workflows, these commands can be piped together:

```bash
# Efficient pipeline (single pass through data)
bcftools +split-vep -c - -p VEP_ input.vcf.gz | \
  bcftools query -HH -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT:%DP:%GQ:%AD]\n' | \
  wombat filter - -F config.yml -o output
```

**Note**: For multiple filter configurations, it's more efficient to save the intermediate TSV file rather than regenerating it each time.

### gnomAD Annotations

For frequency filtering, annotate with gnomAD:

```bash
# Using bcftools annotate
bcftools annotate -a gnomad.genomes.vcf.gz \
  -c INFO/AF,INFO/fafmax_faf95_max \
  input.vcf.gz -o annotated.vcf.gz
```

---

## Configuration Examples

See the [`examples/`](examples/) directory for complete configuration files:

- **[`rare_variants_high_impact.yml`](examples/rare_variants_high_impact.yml)**: Ultra-rare, high-impact variants
- **[`de_novo_mutations.yml`](examples/de_novo_mutations.yml)**: De novo mutation detection with sex-chromosome handling

Each configuration file is fully documented with:

- Parameter descriptions
- Use case recommendations
- Customization tips
- Example command lines

---

## Performance Tips

### For Large Files (>1GB or >50 samples)

1. **Use the two-step workflow**: `wombat prepare` → `wombat filter`
   - Reduces memory usage by 95%+ (4.2M variants → ~100 after early filtering)
   - Pre-expands INFO fields once, reuse for multiple filter configurations
   - Example: 38-sample family with 4.2M variants processes in <1 second with ~1.2GB RAM

2. **Parquet format benefits**:
   - Columnar storage enables selective column loading
   - Pre-filtering before melting (expression filters applied before expanding to per-sample rows)
   - **Per-chromosome processing for DNM**: Automatically processes DNM filtering chromosome-by-chromosome
   - 30% smaller file size vs gzipped TSV

3. **De Novo Mutation (DNM) filtering optimization**:
   - Automatically uses per-chromosome processing when DNM mode is enabled
   - Processes one chromosome at a time to reduce peak memory
   - Applies frequency filters before melting to reduce data expansion
   - Example: 38-sample family with 4.2M variants completes in 20 seconds with ~24GB RAM (vs 200GB+ OOM failure)

### For All Files

3. **Pre-filter with bcftools**: Filter by region/gene before PyWombat
4. **Compressed input**: PyWombat handles `.gz` files natively
5. **Use verbose mode** (`-v`): Monitor progress and filtering statistics

### Memory Comparison

**Expression Filtering** (e.g., VEP_IMPACT filters):

| Approach | 38 samples, 4.2M variants | Memory | Time |
|----------|---------------------------|--------|------|
| Direct TSV | ❌ OOM (>200GB) | 200+ GB | Failed |
| TSV with chunking | ⚠️ Slow | ~30GB | ~3 min |
| **Parquet + pre-filter** | ✅ **Optimal** | **~1.2GB** | **<1 sec** |

**De Novo Mutation (DNM) Filtering**:

| Approach | 38 samples, 4.2M variants | Memory | Time | Result |
|----------|---------------------------|--------|------|--------|
| Without optimization | ❌ OOM (>200GB) | 200+ GB | Failed | N/A |
| **Parquet + per-chromosome** | ✅ **Success** | **~24GB** | **20 sec** | **6,788 DNM variants** |

*DNM filtering requires sample-level data (cannot pre-filter before melting), but per-chromosome processing reduces peak memory by 88%.*

---

## Development

### Setup

```bash
# Clone repository
git clone https://github.com/bourgeron-lab/pywombat.git
cd pywombat

# Install dependencies
uv sync

# Run tests (if available)
uv run pytest
```

### Project Structure

```
pywombat/
├── src/pywombat/
│   ├── __init__.py
│   ├── cli.py           # Main CLI implementation
│   ├── mnv.py           # MNV detection module (optional)
│   └── vep.py           # VEP annotation tools
├── examples/            # Configuration examples
│   ├── README.md
│   ├── rare_variants_high_impact.yml
│   ├── rare_homalt_mnv.yml
│   └── de_novo_mutations.yml
├── tests/               # Test files and data
├── pyproject.toml       # Project metadata
└── README.md
```

### Technology Stack

- **[Polars](https://pola.rs/)**: Fast DataFrame operations
- **[Click](https://click.palletsprojects.com/)**: CLI interface
- **[PyYAML](https://pyyaml.org/)**: Configuration parsing
- **[uv](https://github.com/astral-sh/uv)**: Package management

---

## Troubleshooting

### Common Issues

**Issue**: `Column '(null)' not found`

- **Solution**: Ensure input is bcftools tabulated format with INFO column

**Issue**: `No parent genotypes found`

- **Solution**: Check pedigree file format and sample name matching

**Issue**: `DNM filter requires pedigree`

- **Solution**: Add `--pedigree` option when using de novo config

**Issue**: Variants missing from output

- **Solution**: Use `--debug chr:pos` to see why variants are filtered

**Issue**: Memory errors on large files

- **Solution**: Use the two-step workflow: `wombat prepare` then `wombat filter` for 95%+ memory reduction

**Issue**: Command not found after upgrading

- **Solution**: PyWombat now uses subcommands - use `wombat filter` instead of just `wombat`

### Getting Help

1. Check `--help` for command options: `wombat --help` or `wombat filter --help`
2. Review example configurations in [`examples/`](examples/)
3. Use `--debug` mode to inspect specific variants
4. Use `--verbose` to see filtering steps

---

## Citation

If you use PyWombat in your research, please cite:

```
PyWombat: A tool for processing and filtering bcftools tabulated files
Bourgeron Lab, Institut Pasteur
https://github.com/bourgeron-lab/pywombat
```

---

## License

MIT License - see [LICENSE](LICENSE) file for details

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---

## Authors

- **Freddy Cliquet** - [Bourgeron Lab](https://research.pasteur.fr/en/team/human-genetics-and-cognitive-functions/), Institut Pasteur

---

## Acknowledgments

- Bourgeron Lab for project support
- Institut Pasteur for infrastructure
- Polars team for the excellent DataFrame library
