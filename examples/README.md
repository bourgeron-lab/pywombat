# PyWombat Configuration Examples

This directory contains example configuration files demonstrating different filtering strategies for variant analysis.

## Available Configurations

### 1. `rare_high_impact.yml`

**Purpose:** Filter for ultra-rare variants with high functional impact

**Use Case:** Rare disease gene discovery, identifying causal variants in Mendelian disorders

**Key Features:**
- High-impact variants only (VEP: HIGH)
- Loss-of-function (LoF) with high confidence
- Ultra-rare frequency (≤0.1% in gnomAD) or pathogenic annotation
- Variant-type-specific quality (SNV: AQ>15, INDEL: AQ>22)
- Stringent quality filters (DP≥10, GQ≥19)

**Example Usage:**
```bash
# Single sample or cohort analysis
uvx pywombat filter input.tsv -F examples/rare_high_impact.yml -o output

# With prepared Parquet for large files
uvx pywombat filter prepared.parquet -F examples/rare_high_impact.yml -o output
```

**Expected Annotations:**
- VEP annotations (CANONICAL, IMPACT, LoF, LoF_flags, CLIN_SIG)
- gnomAD v4 (fafmax_faf95_max_genomes, nhomalt_genomes, genomes_filters)

---

### 2. `rare_homalt.yml`

**Purpose:** Filter for homozygous alternative (1/1) variants only

**Use Case:** Autosomal recessive disease analysis, homozygosity mapping, consanguineous families

**Key Features:**
- **Only 1/1 genotypes** (excludes heterozygous and multi-allelic)
- Ultra-rare frequency (≤0.1% in gnomAD) or pathogenic annotation
- High-impact variants in canonical transcripts
- Stringent VAF threshold (≥98% for homozygous calls)
- Parents NOT filtered (allows heterozygous carriers)

**Example Usage:**
```bash
# Recessive disease analysis
uvx pywombat filter input.tsv -F examples/rare_homalt.yml -o output

# With pedigree for trio analysis
uvx pywombat filter input.tsv --pedigree pedigree.tsv \
  -F examples/rare_homalt.yml -o recessive_variants
```

**Key Configuration:**
```yaml
quality:
  homalt_only: true  # Only keep 1/1 genotypes
  sample_vaf_homalt_min: 0.98  # Stringent VAF for true hom calls
```

---

### 3. `rare_coding.yml`

**Purpose:** Filter for all coding variants (HIGH, MODERATE, and synonymous)

**Use Case:** Comprehensive coding variant analysis, including silent changes for conservation studies

**Key Features:**
- HIGH and MODERATE impact variants
- Synonymous variants (for reference)
- Canonical transcripts only
- Pass gnomAD filters with variant-type-specific quality

**Example Usage:**
```bash
uvx pywombat filter input.tsv -F examples/rare_coding.yml -o output
```

---

### 4. `rare_missense.yml`

**Purpose:** Filter for rare missense variants with high pathogenicity scores

**Use Case:** Missense variant analysis using computational predictors

**Key Features:**
- Missense variants in canonical transcripts
- Ultra-rare frequency (≤0.1% in gnomAD)
- Multiple pathogenicity scores (MPC, REVEL, CADD, AlphaMissense, etc.)
- ClinVar pathogenic/likely pathogenic annotations

**Pathogenicity Predictors:**
- MPC ≥ 2.6
- AlphaMissense ≥ 0.9
- REVEL ≥ 0.773
- CADD v1.7 ≥ 28.3
- MISTIC ≥ 0.8
- popEVE ≤ -4.6

**Example Usage:**
```bash
uvx pywombat filter input.tsv -F examples/rare_missense.yml -o output
```

---

### 5. `rare_spliceai.yml`

**Purpose:** Filter for variants predicted to affect splicing

**Use Case:** Splicing variant analysis using SpliceAI predictions

**Key Features:**
- Ultra-rare frequency (≤0.1% in gnomAD)
- High SpliceAI score (≥0.8) OR pass quality filters
- Predicts donor/acceptor gain/loss

**Example Usage:**
```bash
uvx pywombat filter input.tsv -F examples/rare_spliceai.yml -o output
```

---

### 6. `de_novo_mutations.yml`

**Purpose:** Identify de novo mutations in trio/family data

**Use Case:** Autism, developmental disorders, sporadic diseases

**Key Features:**
- Sex-chromosome aware (X/Y with PAR region handling)
- Stringent parent filters (VAF<2%, DP≥10, GQ≥18)
- Hemizygous variant detection (VAF≥85%)
- Population frequency filtering (≤0.1% in gnomAD)
- Quality filter enforcement (genomes_filters PASS only)

**Example Usage:**
```bash
# Basic de novo detection with trio
uvx pywombat filter input.tsv --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml -o denovo

# With verbose output to track filtering steps
uvx pywombat filter input.tsv --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml -o denovo --verbose
```

**Required Files:**
1. **Input TSV:** bcftools-formatted tabulated file with all family members
2. **Pedigree file:** Tab-separated file with family relationships and sex information

**Pedigree Format:**
```tsv
FID	sample_id	FatherBarcode	MotherBarcode	Sex	Pheno
FAM1	Proband1	Father1	Mother1	1	2
FAM1	Father1	0	0	1	1
FAM1	Mother1	0	0	2	1
```

- `FID`: Family identifier
- `sample_id`: Sample name (must match VCF sample names)
- `FatherBarcode`: Father's sample name (0 = unknown)
- `MotherBarcode`: Mother's sample name (0 = unknown)
- `Sex`: 1=male, 2=female (or M/F)
- `Pheno`: Phenotype (1=unaffected, 2=affected)

---

## Quality Filter Options

All example configurations support the following quality filter options:

### Standard Options
```yaml
quality:
  sample_dp_min: 10            # Minimum read depth
  sample_gq_min: 19            # Minimum genotype quality
  sample_vaf_het_min: 0.25     # Het VAF minimum (25%)
  sample_vaf_het_max: 0.75     # Het VAF maximum (75%)
  sample_vaf_homalt_min: 0.85  # Hom alt VAF minimum (85%)
  sample_vaf_hom_ref_max: 0.15 # Hom ref VAF maximum (15%)

  apply_to_parents: false      # Apply filters to parents
  filter_no_alt_allele: true   # Exclude 0/0 and ./.
```

### Recessive Analysis Option
```yaml
quality:
  homalt_only: true  # Only keep 1/1 genotypes (NEW in v1.3.1)
```

**When to use `homalt_only`:**
- Autosomal recessive disease analysis
- Homozygosity mapping in consanguineous families
- Identifying regions of identical-by-descent
- Compound heterozygote analysis (first pass for 1/1, second pass for 0/1)

---

## Expression Filter Operators

PyWombat supports powerful expression-based filtering with the following operators:

### Comparison Operators
- `=`, `!=`: Equality and inequality
- `<`, `>`, `<=`, `>=`: Numeric comparisons

### Logical Operators
- `&`: AND
- `|`: OR
- `(`, `)`: Grouping

### Special Operators (NEW in v1.3.0)
- **`contains`**: Case-insensitive substring matching
  - Example: `VEP_CLIN_SIG contains 'pathogenic'`

- **`is_empty`**: Checks for null, empty string, or "."
  - Example: `genomes_filters is_empty`

- **`is_snv`**: TRUE when both REF and ALT are single nucleotides (A/C/G/T)
  - Example: `is_snv & QUAL > 30`

- **`is_indel`**: TRUE for insertions, deletions, and MNVs (not SNVs)
  - Example: `is_indel & VEP_IMPACT = HIGH`

### Null Checks
- `= null`: Check if value is null
- `!= null`: Check if value is not null

### Example Expressions

```yaml
# High or moderate impact with low frequency
expression: "(VEP_IMPACT = HIGH | VEP_IMPACT = MODERATE) & gnomad_AF < 0.001"

# Pathogenic variants that passed filters
expression: "VEP_CLIN_SIG contains 'pathogenic' & genomes_filters is_empty"

# Different quality for SNVs vs INDELs
expression: "(is_snv & QUAL > 18) | (is_indel & QUAL > 22)"

# Canonical transcripts with CADD score
expression: "VEP_CANONICAL = YES & CADD_PHRED >= 20"

# Complex multi-criteria filter
expression: >
  VEP_CANONICAL = YES &
  VEP_IMPACT = HIGH &
  (fafmax_faf95_max_genomes = null | fafmax_faf95_max_genomes <= 0.001) &
  genomes_filters is_empty
```

---

## Customizing Configurations

All YAML configuration files can be modified to suit your specific needs:

### Adjusting Quality Thresholds

```yaml
quality:
  sample_dp_min: 15        # Increase for higher coverage datasets
  sample_gq_min: 20        # Increase for more stringent quality
  sample_vaf_het_min: 0.30 # Tighten het VAF range
  sample_vaf_homalt_min: 0.90  # More stringent hom calls
```

### Modifying Frequency Cutoffs

```yaml
# In rare_high_impact.yml
expression: "... & ( fafmax_faf95_max_genomes = null | fafmax_faf95_max_genomes <= 0.01 )"  # 1% instead of 0.1%

# In de_novo_mutations.yml
dnm:
  fafmax_faf95_max_genomes_max: 0.01  # Allow more common variants
```

### Changing PAR Regions (GRCh37/hg19)

If using GRCh37/hg19 instead of GRCh38, update PAR coordinates in `de_novo_mutations.yml`:

```yaml
par_regions:
  grch37:
    PAR1:
      chrom: X
      start: 60001
      end: 2699520
    PAR2:
      chrom: X
      start: 154931044
      end: 155260560
```

---

## Output Files

All configurations produce filtered TSV files with:
- All original VCF/TSV columns
- Expanded annotation fields (from INFO/FORMAT)
- Sample-specific columns (sample_gt, sample_dp, sample_gq, sample_vaf)
- Parent columns (father_*, mother_*) when using pedigree

**Example output columns:**
```
#CHROM  POS  REF  ALT  sample  sample_gt  sample_dp  sample_gq  sample_vaf  VEP_SYMBOL  VEP_IMPACT  ...
chr1    12345  A    G    Child1  0/1        45         99         0.47        GENE1       HIGH        ...
```

---

## Tips and Best Practices

### For Rare Variant Analysis:
1. Start with the provided thresholds and adjust based on your data quality
2. Review filtered variants for known disease genes
3. Consider adding gene lists or inheritance pattern filters
4. Validate top candidates with orthogonal methods

### For Recessive Disease Analysis:
1. Use `rare_homalt.yml` as a starting point
2. Ensure both parents are in the pedigree file
3. Consider combining 1/1 analysis with compound heterozygote detection
4. Look for runs of homozygosity (ROH) in consanguineous families

### For De Novo Analysis:
1. Always verify both parents are present in the VCF
2. Check for sample swaps using known variants or fingerprinting
3. Validate DNMs with visual inspection (IGV) or Sanger sequencing
4. Be cautious of repetitive regions and segmental duplications
5. Consider parental age effects and known DNM hotspots

### Performance Optimization:
- For large cohorts (>50 samples), use the two-step workflow:
  1. `wombat prepare input.tsv.gz -o prepared.parquet`
  2. `wombat filter prepared.parquet -F config.yml -o output`
- Use `--verbose` to monitor filtering steps
- Pre-filter VCF with bcftools for specific regions/genes if needed

---

## Creating Your Own Configuration

Create a new YAML file with the following structure:

```yaml
# Quality filters (always applied first)
quality:
  sample_dp_min: 10
  sample_gq_min: 18
  sample_vaf_het_min: 0.25
  sample_vaf_het_max: 0.75
  sample_vaf_homalt_min: 0.85
  filter_no_alt_allele: true

  # Optional: only keep homozygous alternative (1/1)
  # homalt_only: false

# Expression-based filter (uses Polars expressions)
expression: "VEP_IMPACT = HIGH & gnomad_AF < 0.01"

# OR: De novo mode (mutually exclusive with expression)
dnm:
  enabled: true
  parent_dp_min: 10
  parent_gq_min: 18
  parent_vaf_max: 0.02
  # ... additional DNM settings
```

---

## Questions or Issues?

For more information:
- See main [README.md](../README.md) for installation and basic usage
- Check [PyWombat on GitHub](https://github.com/bourgeron-lab/pywombat) for documentation
- Report issues on GitHub

---

## Version History

### v1.3.1 (Current)
- Added `homalt_only` quality filter option for recessive analysis
- Enhanced documentation across all example files
- Added comprehensive headers to all configuration files

### v1.3.0
- Added new expression operators: `contains`, `is_empty`, `is_snv`, `is_indel`
- Improved multi-line expression support in YAML

### v1.2.x
- Per-chromosome DNM processing for memory optimization
- Two-step workflow (prepare → filter) for large files
