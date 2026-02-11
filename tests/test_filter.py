"""Tests for the filter command."""

import polars as pl
from click.testing import CliRunner

from pywombat.cli import cli, apply_filters_lazy


class TestFilterCommand:
    """Test the filter subcommand."""

    def test_filter_tsv_input(self, small_tsv, tmp_test_dir):
        """Test filtering TSV input."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["filter", str(small_tsv), "-o", str(output_prefix)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.tsv"
        assert output_path.exists(), "Output TSV file was not created"

    def test_filter_parquet_input(self, small_tsv, tmp_test_dir):
        """Test filtering Parquet input (prepared file)."""
        # First, prepare the Parquet file
        parquet_path = tmp_test_dir / "prepared.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(parquet_path)]
        )
        assert result.exit_code == 0, f"Prepare failed: {result.output}"

        # Now filter the Parquet file
        output_prefix = tmp_test_dir / "filtered"
        result = runner.invoke(
            cli, ["filter", str(parquet_path), "-o", str(output_prefix)]
        )

        assert result.exit_code == 0, f"Filter failed: {result.output}"

        output_path = tmp_test_dir / "filtered.tsv"
        assert output_path.exists(), "Output TSV file was not created"

    def test_filter_with_pedigree(self, small_tsv, small_pedigree, tmp_test_dir):
        """Test filtering with pedigree file."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "filter",
                str(small_tsv),
                "-o",
                str(output_prefix),
                "-p",
                str(small_pedigree),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Check that output contains parent genotype columns
        output_path = tmp_test_dir / "output.tsv"
        df = pl.read_csv(output_path, separator="\t")

        # Should have father_gt and mother_gt columns
        assert "father_gt" in df.columns or "father_id" in df.columns, \
            "Parent columns should be present when pedigree is provided"

    def test_filter_parquet_output_format(self, small_tsv, tmp_test_dir):
        """Test that Parquet output format works."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["filter", str(small_tsv), "-o", str(output_prefix), "-f", "parquet"],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.parquet"
        assert output_path.exists(), "Parquet output file was not created"

        # Verify it's valid Parquet
        df = pl.read_parquet(output_path)
        assert len(df) > 0, "Parquet file should contain data"

    def test_filter_gzipped_output_format(self, small_tsv, tmp_test_dir):
        """Test that gzipped TSV output format works."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["filter", str(small_tsv), "-o", str(output_prefix), "-f", "tsv.gz"],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.tsv.gz"
        assert output_path.exists(), "Gzipped output file was not created"

    def test_filter_melts_samples(self, small_tsv, tmp_test_dir):
        """Test that samples are melted into rows."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["filter", str(small_tsv), "-o", str(output_prefix)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.tsv"
        df = pl.read_csv(output_path, separator="\t")

        # Should have 'sample' column (melted)
        assert "sample" in df.columns, "Sample column should exist after melting"

        # Should have sample genotype columns
        assert "sample_gt" in df.columns, "sample_gt column should exist"
        assert "sample_dp" in df.columns, "sample_dp column should exist"

    def test_filter_extracts_vaf(self, small_tsv, tmp_test_dir):
        """Test that VAF is calculated correctly."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["filter", str(small_tsv), "-o", str(output_prefix)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.tsv"
        df = pl.read_csv(output_path, separator="\t")

        assert "sample_vaf" in df.columns, "sample_vaf column should exist"

    def test_filter_verbose_mode(self, small_tsv, tmp_test_dir):
        """Test verbose output mode."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["filter", str(small_tsv), "-o", str(output_prefix), "-v"]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Verbose output should contain status messages
        assert "Reading" in result.output or "Processing" in result.output

    def test_filter_default_output_name(self, small_tsv, tmp_test_dir):
        """Test that default output name is derived from input."""
        runner = CliRunner()

        # Change to temp dir so output goes there
        with runner.isolated_filesystem(temp_dir=tmp_test_dir):
            # Copy the small_tsv to current dir
            import shutil
            shutil.copy(small_tsv, "test_input.tsv")

            result = runner.invoke(cli, ["filter", "test_input.tsv"])

            assert result.exit_code == 0, f"Command failed: {result.output}"
            # Should create test_input.tsv as output


class TestFilterWithConfig:
    """Test filter command with configuration files."""

    def test_filter_with_quality_config(
        self, small_tsv, filter_config, tmp_test_dir
    ):
        """Test filtering with quality filter config."""
        output_prefix = tmp_test_dir / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "filter",
                str(small_tsv),
                "-o",
                str(output_prefix),
                "-F",
                str(filter_config),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        output_path = tmp_test_dir / "output.tsv"
        assert output_path.exists(), "Output file was not created"


class TestDNMOptimization:
    """Test DNM filtering optimizations."""

    def test_get_unique_chromosomes(self, tmp_test_dir):
        """Test chromosome discovery from Parquet file."""
        import polars as pl
        from pywombat.cli import get_unique_chromosomes

        # Create a test Parquet file with multiple chromosomes
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr2", "chr1", "chrX", "chr22", "chrY", "chr10"],
            "POS": [100, 200, 300, 400, 500, 600, 700],
            "REF": ["A", "C", "G", "T", "A", "C", "G"],
            "ALT": ["T", "G", "C", "A", "T", "G", "C"],
        })

        parquet_path = tmp_test_dir / "test.parquet"
        test_data.write_parquet(parquet_path)

        # Get unique chromosomes
        chroms = get_unique_chromosomes(parquet_path)

        # Should be sorted naturally: 1, 2, 10, 22, X, Y
        assert chroms == ["chr1", "chr2", "chr10", "chr22", "chrX", "chrY"], \
            f"Expected natural chromosome order, got {chroms}"

    def test_apply_dnm_prefilters(self, tmp_test_dir):
        """Test that frequency prefilters reduce variant count."""
        import polars as pl
        from pywombat.cli import apply_dnm_prefilters

        # Create test data with varying frequencies
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [100, 200, 300, 400],
            "fafmax_faf95_max_genomes": [0.0001, 0.002, None, 0.00005],
            "genomes_filters": [".", "FAIL", ".", "."],
        })

        # Config with frequency filter
        config = {
            "dnm": {
                "fafmax_faf95_max_genomes_max": 0.001,
                "genomes_filters_pass_only": True
            }
        }

        # Apply prefilters
        lazy_df = test_data.lazy()
        filtered_lazy = apply_dnm_prefilters(lazy_df, config, verbose=False)
        result = filtered_lazy.collect()

        # Should keep only rows where:
        # - fafmax <= 0.001 OR NULL
        # - genomes_filters == "." OR NULL
        # Row 0: 0.0001 <= 0.001 AND "." -> PASS
        # Row 1: 0.002 > 0.001 -> FAIL
        # Row 2: NULL AND "." -> PASS
        # Row 3: 0.00005 <= 0.001 AND "." -> PASS
        assert result.shape[0] == 3, f"Expected 3 rows after prefilter, got {result.shape[0]}"

    def test_dnm_skip_prefilters(self):
        """Test that skip_prefilters parameter works."""
        import polars as pl
        from pywombat.cli import apply_de_novo_filter

        # Create minimal test data with all required columns
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "sample": ["Sample1", "Sample1"],
            "sample_gt": ["0/1", "0/1"],
            "sample_dp": [30, 35],
            "sample_gq": [50, 50],
            "sample_vaf": [0.45, 0.48],
            "father_id": ["Father1", "Father1"],
            "mother_id": ["Mother1", "Mother1"],
            "father_gt": ["0/0", "0/0"],
            "father_dp": [25, 28],
            "father_gq": [50, 50],
            "father_vaf": [0.01, 0.01],
            "mother_gt": ["0/0", "0/0"],
            "mother_dp": [22, 26],
            "mother_gq": [50, 50],
            "mother_vaf": [0.01, 0.01],
            "sex": ["1", "1"],
            "fafmax_faf95_max_genomes": [0.002, 0.0005],
        })

        # DNM config with frequency filter
        dnm_config = {
            "sample_dp_min": 10,
            "sample_gq_min": 18,
            "sample_vaf_min": 0.2,
            "parent_dp_min": 10,
            "parent_gq_min": 18,
            "parent_vaf_max": 0.02,
            "fafmax_faf95_max_genomes_max": 0.001,
        }

        # Apply with skip_prefilters=True (should NOT filter by frequency)
        result_skip = apply_de_novo_filter(
            test_data, dnm_config, verbose=False,
            skip_prefilters=True
        )

        # Apply with skip_prefilters=False (should filter by frequency)
        result_no_skip = apply_de_novo_filter(
            test_data, dnm_config, verbose=False,
            skip_prefilters=False
        )

        # With skip_prefilters=True, both rows pass (no frequency filter applied)
        # With skip_prefilters=False, only row 2 passes (0.0005 <= 0.001)
        assert result_skip.shape[0] == 2, \
            f"Expected 2 rows with skip_prefilters=True, got {result_skip.shape[0]}"
        assert result_no_skip.shape[0] == 1, \
            f"Expected 1 row with skip_prefilters=False, got {result_no_skip.shape[0]}"


class TestHomaltFilter:
    """Test homalt_only quality filter."""

    def test_homalt_only_filter_excludes_heterozygous(self):
        """Test that homalt_only=true excludes heterozygous (0/1) variants."""
        # Create test data with mixed genotypes
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [100, 200, 300, 400],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
            "sample_gt": ["1/1", "0/1", "1/1", "0/0"],
            "sample_dp": [20, 20, 20, 20],
            "sample_gq": [30, 30, 30, 30],
            "sample_vaf": [0.98, 0.50, 0.99, 0.02],
        })

        # Apply homalt_only filter using lazy path
        result = apply_filters_lazy(
            test_data.lazy(),
            {"quality": {"homalt_only": True}},
            verbose=False,
            pedigree_df=None
        ).collect()

        # Should only keep variants with 1/1 genotypes
        assert result.shape[0] == 2, \
            f"Expected 2 variants (1/1 only), got {result.shape[0]}"
        assert all(gt == "1/1" for gt in result["sample_gt"]), \
            "All genotypes should be 1/1"
        assert set(result["POS"].to_list()) == {100, 300}, \
            "Should keep positions 100 and 300 only"

    def test_homalt_only_with_quality_filters(self):
        """Test homalt_only combined with other quality filters."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [100, 200, 300, 400],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
            "sample_gt": ["1/1", "1/1", "1/1", "0/1"],
            "sample_dp": [5, 15, 20, 25],  # Low depth for pos 100
            "sample_gq": [30, 30, 30, 30],
            "sample_vaf": [0.98, 0.99, 0.98, 0.50],
        })

        # Apply homalt_only with minimum depth filter
        result = apply_filters_lazy(
            test_data.lazy(),
            {"quality": {"homalt_only": True, "sample_dp_min": 10}},
            verbose=False,
            pedigree_df=None
        ).collect()

        # Should keep only 1/1 with dp >= 10
        assert result.shape[0] == 2, \
            f"Expected 2 variants (1/1 and dp >= 10), got {result.shape[0]}"
        assert set(result["POS"].to_list()) == {200, 300}, \
            "Should keep positions 200 and 300 only"
        assert all(gt == "1/1" for gt in result["sample_gt"]), \
            "All genotypes should be 1/1"


class TestMNVOnlyFilter:
    """Test mnv.only filter option."""

    def test_mnv_only_keeps_variants_with_mnv_proba(self):
        """Test that mnv.only=true keeps only variants with non-null mnv_proba."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [100, 101, 200, 300],
            "REF": ["A", "T", "C", "G"],
            "ALT": ["G", "A", "T", "A"],
            "sample": ["S1", "S1", "S1", "S1"],
            "sample_gt": ["0/1", "0/1", "0/1", "0/1"],
            "sample_dp": [20, 20, 20, 20],
            "sample_gq": [30, 30, 30, 30],
            "sample_vaf": [0.50, 0.50, 0.50, 0.50],
            "mnv_proba": [0.95, 0.95, None, 0.87],
            "mnv_variant": ["chr1:100:AT:GA", "chr1:100:AT:GA", "", "chr1:300:GC:AT"],
            "mnv_inframe": [None, None, None, None],
        })

        # Apply MNV-only filter (same logic as cli.py)
        result = test_data.filter(pl.col("mnv_proba").is_not_null())

        assert result.shape[0] == 3, \
            f"Expected 3 variants with mnv_proba, got {result.shape[0]}"
        assert set(result["POS"].to_list()) == {100, 101, 300}, \
            "Should keep positions 100, 101, and 300 (all with mnv_proba values)"

    def test_mnv_only_keeps_both_snv_and_indel_candidates(self):
        """Test that mnv.only keeps both SNV and indel MNV candidates."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1", "chr2", "chr2", "chr3"],
            "POS": [100, 101, 500, 505, 700],
            "REF": ["A", "T", "ATG", "C", "G"],
            "ALT": ["G", "A", "A", "CGAT", "T"],
            "sample": ["S1", "S1", "S1", "S1", "S1"],
            "sample_gt": ["0/1", "0/1", "0/1", "0/1", "0/1"],
            "sample_dp": [20, 20, 20, 20, 20],
            "sample_gq": [30, 30, 30, 30, 30],
            "sample_vaf": [0.50, 0.50, 0.50, 0.50, 0.50],
            "mnv_proba": [0.95, 0.95, 0.88, 0.88, None],
            "mnv_variant": ["chr1:100:AT:GA", "chr1:100:AT:GA", "", "", ""],
            "mnv_inframe": [None, None, True, True, None],
        })

        result = test_data.filter(pl.col("mnv_proba").is_not_null())

        assert result.shape[0] == 4, \
            f"Expected 4 MNV candidates (2 SNV + 2 indel), got {result.shape[0]}"
        # SNV candidates (pos 100, 101) and indel candidates (pos 500, 505)
        assert set(result["POS"].to_list()) == {100, 101, 500, 505}

    def test_mnv_only_returns_empty_when_no_candidates(self):
        """Test that mnv.only returns empty DataFrame when no MNV candidates exist."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "sample": ["S1", "S1"],
            "sample_gt": ["0/1", "0/1"],
            "sample_dp": [20, 20],
            "sample_gq": [30, 30],
            "sample_vaf": [0.50, 0.50],
            "mnv_proba": [None, None],
            "mnv_variant": ["", ""],
            "mnv_inframe": [None, None],
        })

        result = test_data.filter(pl.col("mnv_proba").is_not_null())

        assert result.shape[0] == 0, \
            f"Expected 0 variants (no MNV candidates), got {result.shape[0]}"

    def test_mnv_only_preserves_all_columns(self):
        """Test that mnv.only filter preserves all columns in the output."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "FILTER": ["PASS", "PASS"],
            "VEP_SYMBOL": ["BRCA1", "TP53"],
            "sample": ["S1", "S1"],
            "sample_gt": ["0/1", "0/1"],
            "sample_dp": [20, 20],
            "sample_gq": [30, 30],
            "sample_vaf": [0.50, 0.50],
            "mnv_proba": [0.95, None],
            "mnv_variant": ["chr1:100:AC:GT", ""],
            "mnv_inframe": [None, None],
        })

        result = test_data.filter(pl.col("mnv_proba").is_not_null())

        assert result.shape[0] == 1
        assert result.columns == test_data.columns, \
            "All columns should be preserved"
        assert result["VEP_SYMBOL"][0] == "BRCA1"


class TestExcludeMHCFilter:
    """Test mnv.exclude_mhc filter option."""

    @staticmethod
    def _apply_mhc_filter(df, exclude_mhc):
        """Reproduce the MHC exclusion logic from cli.py."""
        mhc_chrom = exclude_mhc.get("chrom", "6")
        mhc_start = exclude_mhc.get("start", 28510120)
        mhc_end = exclude_mhc.get("end", 33480577)
        chrom_col = pl.col("#CHROM").cast(pl.Utf8)
        is_mhc = (
            (chrom_col == mhc_chrom)
            | (chrom_col == f"chr{mhc_chrom}")
        ) & (pl.col("POS") >= mhc_start) & (pl.col("POS") <= mhc_end)
        return df.filter(~is_mhc)

    def test_excludes_mhc_variants_chr_prefix(self):
        """Test that MHC variants with chr6 prefix are excluded."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr6", "chr6", "chr6"],
            "POS": [100, 30000000, 28510120, 34000000],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
        })
        result = self._apply_mhc_filter(
            test_data, {"enabled": True, "chrom": "6", "start": 28510120, "end": 33480577}
        )
        # chr1:100 (not chr6), chr6:34000000 (outside MHC) should remain
        assert result.shape[0] == 2
        assert set(result["POS"].to_list()) == {100, 34000000}

    def test_excludes_mhc_variants_no_chr_prefix(self):
        """Test that MHC variants without chr prefix are excluded."""
        test_data = pl.DataFrame({
            "#CHROM": ["1", "6", "6", "6"],
            "POS": [100, 30000000, 33480577, 33480578],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
        })
        result = self._apply_mhc_filter(
            test_data, {"enabled": True, "chrom": "6", "start": 28510120, "end": 33480577}
        )
        # 1:100 (not chr6), 6:33480578 (outside MHC) should remain
        assert result.shape[0] == 2
        assert set(result["POS"].to_list()) == {100, 33480578}

    def test_mhc_boundary_inclusive(self):
        """Test that MHC boundaries are inclusive (start and end are excluded)."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr6", "chr6", "chr6", "chr6"],
            "POS": [28510119, 28510120, 33480577, 33480578],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
        })
        result = self._apply_mhc_filter(
            test_data, {"enabled": True, "chrom": "6", "start": 28510120, "end": 33480577}
        )
        # Only positions just outside the boundaries should remain
        assert result.shape[0] == 2
        assert set(result["POS"].to_list()) == {28510119, 33480578}

    def test_no_mhc_variants_unchanged(self):
        """Test that DataFrame without MHC variants is unchanged."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr1", "chr2", "chr3"],
            "POS": [100, 200, 300],
            "REF": ["A", "C", "G"],
            "ALT": ["G", "T", "A"],
        })
        result = self._apply_mhc_filter(
            test_data, {"enabled": True, "chrom": "6", "start": 28510120, "end": 33480577}
        )
        assert result.shape[0] == 3

    def test_default_grch38_coordinates(self):
        """Test that default coordinates match GRCh38 xMHC region."""
        test_data = pl.DataFrame({
            "#CHROM": ["chr6", "chr6"],
            "POS": [30000000, 40000000],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
        })
        # Use defaults (no explicit start/end)
        result = self._apply_mhc_filter(test_data, {"enabled": True, "chrom": "6"})
        # 30M is inside default MHC, 40M is outside
        assert result.shape[0] == 1
        assert result["POS"][0] == 40000000
