# PyWombat Configuration Examples

This directory contains example configuration files demonstrating different filtering strategies for variant analysis.

## Available Configurations

### 1. `rare_variants_high_impact.yml`

**Purpose:** Filter for ultra-rare variants with high functional impact

**Use Case:** Rare disease gene discovery, identifying causal variants in Mendelian disorders

**Key Features:**
- High-impact variants only (VEP: HIGH)
- Loss-of-function (LoF) with high confidence
- Ultra-rare frequency (≤0.1% in gnomAD)
- Stringent quality filters (DP≥10, GQ≥19)
- Appropriate VAF thresholds for het/hom calls

**Example Usage:**
```bash
# Single sample or cohort analysis
uvx pywombat input.tsv -F examples/rare_variants_high_impact.yml -o output

# With verbose output
uvx pywombat input.tsv -F examples/rare_variants_high_impact.yml -o output --verbose

# Output as compressed TSV
uvx pywombat input.tsv -F examples/rare_variants_high_impact.yml -o output -f tsv.gz
```

**Expected Annotations:**
- VEP annotations (CANONICAL, IMPACT, LoF, LoF_flags)
- gnomAD v4 frequency annotations (fafmax_faf95_max_genomes)

---

### 2. `de_novo_mutations.yml`

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
uvx pywombat input.tsv --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml -o output

# With verbose output to track filtering steps
uvx pywombat input.tsv --pedigree pedigree.tsv \
  -F examples/de_novo_mutations.yml -o output --verbose
```

**Required Files:**
1. **Input VCF/TSV:** bcftools-formatted tabulated file with all family members
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

**Expected Annotations:**
- Basic VCF genotype fields (GT, DP, GQ, AD)
- gnomAD annotations (fafmax_faf95_max_genomes, genomes_filters)
- Optional: VEP annotations for downstream filtering

---

## Customizing Configurations

All YAML configuration files can be modified to suit your specific needs:

### Adjusting Quality Thresholds

```yaml
quality:
  sample_dp_min: 15        # Increase for higher coverage datasets
  sample_gq_min: 20        # Increase for more stringent quality
  sample_vaf_het_min: 0.30 # Tighten het VAF range
```

### Modifying Frequency Cutoffs

```yaml
# In rare_variants_high_impact.yml
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

### For De Novo Analysis:
1. Always verify both parents are present in the VCF
2. Check for sample swaps using known variants or fingerprinting
3. Validate DNMs with visual inspection (IGV) or Sanger sequencing
4. Be cautious of repetitive regions and segmental duplications
5. Consider parental age effects and known DNM hotspots

### Performance Optimization:
- For large cohorts, use `-f parquet` for faster downstream analysis
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

**Available operators in expressions:**
- Comparison: `=`, `!=`, `<`, `>`, `<=`, `>=`
- Logical: `&` (AND), `|` (OR)
- Grouping: `(` and `)`
- Null checks: `= null`, `!= null`

**Example expressions:**
```yaml
# High impact OR moderate impact with low frequency
expression: "(VEP_IMPACT = HIGH | VEP_IMPACT = MODERATE) & gnomad_AF < 0.001"

# Canonical transcripts only, with CADD score
expression: "VEP_CANONICAL = YES & CADD_PHRED >= 20"

# Specific consequence types
expression: "VEP_Consequence = frameshift_variant | VEP_Consequence = stop_gained"
```

---

## Questions or Issues?

For more information:
- See main [README.md](../README.md) for installation and basic usage
- Check [PyWombat documentation](https://github.com/bourgeron-lab/pywombat)
- Report issues on GitHub
