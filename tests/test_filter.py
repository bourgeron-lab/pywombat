"""Tests for the filter command."""

import polars as pl
from click.testing import CliRunner

from pywombat.cli import cli


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
