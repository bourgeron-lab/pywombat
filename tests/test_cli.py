"""Tests for the CLI interface."""

from click.testing import CliRunner

from pywombat.cli import cli


class TestCLIStructure:
    """Test the CLI command structure."""

    def test_cli_help(self):
        """Test that CLI help displays correctly."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "Wombat" in result.output
        assert "filter" in result.output
        assert "prepare" in result.output

    def test_filter_help(self):
        """Test that filter subcommand help displays correctly."""
        runner = CliRunner()
        result = runner.invoke(cli, ["filter", "--help"])
        assert result.exit_code == 0
        assert "Process and filter" in result.output
        assert "--pedigree" in result.output
        assert "--filter-config" in result.output

    def test_prepare_help(self):
        """Test that prepare subcommand help displays correctly."""
        runner = CliRunner()
        result = runner.invoke(cli, ["prepare", "--help"])
        assert result.exit_code == 0
        assert "Convert bcftools TSV" in result.output
        assert "--output" in result.output
        assert "--chunk-size" in result.output

    def test_old_syntax_fails(self):
        """Test that old syntax (without subcommand) fails with helpful error."""
        runner = CliRunner()
        # Trying to pass a file directly without subcommand should fail
        result = runner.invoke(cli, ["nonexistent.tsv"])
        # Should show error or usage message since 'nonexistent.tsv' is not a command
        assert result.exit_code != 0

    def test_missing_subcommand(self):
        """Test that running without subcommand shows help."""
        runner = CliRunner()
        result = runner.invoke(cli, [])
        # Should show usage/help
        assert "Usage:" in result.output or "Commands:" in result.output
