"""Tests for the prepare command."""

import polars as pl
from click.testing import CliRunner

from pywombat.cli import cli


class TestPrepareCommand:
    """Test the prepare subcommand."""

    def test_prepare_creates_parquet(self, small_tsv, tmp_test_dir):
        """Test that prepare command creates a Parquet file."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_path.exists(), "Parquet file was not created"

    def test_prepare_extracts_info_fields(self, small_tsv, tmp_test_dir):
        """Test that INFO fields are extracted as separate columns."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Read the Parquet file and check columns
        df = pl.read_parquet(output_path)

        # Check that INFO fields are now columns
        assert "VEP_SYMBOL" in df.columns, "VEP_SYMBOL column missing"
        assert "VEP_IMPACT" in df.columns, "VEP_IMPACT column missing"
        assert "AF" in df.columns, "AF column missing"

        # Check that (null) column is removed
        assert "(null)" not in df.columns, "(null) column should be removed"

    def test_prepare_preserves_samples(self, small_tsv, tmp_test_dir):
        """Test that sample columns are preserved (not melted)."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Read the Parquet file
        df = pl.read_parquet(output_path)

        # Check that sample columns still exist (wide format)
        sample_cols = [c for c in df.columns if "Sample1" in c or "Sample2" in c]
        assert len(sample_cols) > 0, "Sample columns should be preserved"

        # Check row count matches input (3 variants)
        assert len(df) == 3, f"Expected 3 rows, got {len(df)}"

    def test_prepare_gzipped_input(self, small_tsv_gz, tmp_test_dir):
        """Test that prepare handles gzipped input."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv_gz), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert output_path.exists(), "Parquet file was not created"

    def test_prepare_verbose_output(self, small_tsv, tmp_test_dir):
        """Test that verbose flag produces output."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path), "-v"]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Verbose output should contain progress information
        assert "Preparing" in result.output or "Pass" in result.output

    def test_prepare_adds_parquet_extension(self, small_tsv, tmp_test_dir):
        """Test that .parquet extension is added if missing."""
        output_path = tmp_test_dir / "output"  # No extension

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Should have created output.parquet
        expected_path = tmp_test_dir / "output.parquet"
        assert expected_path.exists(), "Parquet file with .parquet extension was not created"

    def test_prepare_memory_efficient_dtypes(self, small_tsv, tmp_test_dir):
        """Test that memory-efficient dtypes are applied."""
        output_path = tmp_test_dir / "output.parquet"

        runner = CliRunner()
        result = runner.invoke(
            cli, ["prepare", str(small_tsv), "-o", str(output_path)]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

        df = pl.read_parquet(output_path)

        # POS should be UInt32 (or compatible integer type)
        assert df.schema["POS"] in [pl.UInt32, pl.Int64, pl.UInt64], \
            f"POS should be integer type, got {df.schema['POS']}"
