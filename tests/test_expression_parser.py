"""Tests for expression parser functionality."""

import polars as pl
import pytest

from pywombat.cli import parse_impact_filter_expression


class TestExpressionParser:
    """Test the expression parser for filter expressions."""

    def test_column_names_with_dots(self):
        """Test that column names with dots are supported."""
        # Create test dataframe with column names containing dots
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2", "chr3"],
            "cadd_v1.7": [25.0, 30.0, 15.0],
            "faf95.joint": [0.001, 0.005, 0.0001],
            "some.nested.col": ["A", "B", "C"],
        })

        # Test comparison operators with dotted column names
        expr = parse_impact_filter_expression("cadd_v1.7 >= 28.3", df)
        result = df.filter(expr)
        assert result.shape[0] == 1
        assert result["cadd_v1.7"][0] == 30.0

        # Test with multiple conditions
        expr = parse_impact_filter_expression(
            "cadd_v1.7 >= 20 & faf95.joint <= 0.001", df
        )
        result = df.filter(expr)
        assert result.shape[0] == 1  # Only chr1 matches both conditions

    def test_not_empty_with_dotted_columns(self):
        """Test not_empty operator with column names containing dots."""
        df = pl.DataFrame({
            "variant_id": ["v1", "v2", "v3", "v4"],
            "score.v1.2": ["0.8", None, "0.9", "."],
        })

        expr = parse_impact_filter_expression("score.v1.2 not_empty", df)
        result = df.filter(expr)
        assert result.shape[0] == 2  # Only v1 and v3
        assert "v1" in result["variant_id"].to_list()
        assert "v3" in result["variant_id"].to_list()

    def test_is_empty_with_dotted_columns(self):
        """Test is_empty operator with column names containing dots."""
        df = pl.DataFrame({
            "variant_id": ["v1", "v2", "v3", "v4"],
            "anno.field": ["value", None, ".", ""],
        })

        expr = parse_impact_filter_expression("anno.field is_empty", df)
        result = df.filter(expr)
        assert result.shape[0] == 3  # v2, v3, v4
        assert "v1" not in result["variant_id"].to_list()

    def test_contains_with_dotted_columns(self):
        """Test contains operator with column names containing dots."""
        df = pl.DataFrame({
            "variant_id": ["v1", "v2", "v3"],
            "pred.consequence": ["missense_variant", "frameshift", "missense_splice"],
        })

        expr = parse_impact_filter_expression("pred.consequence contains 'missense'", df)
        result = df.filter(expr)
        assert result.shape[0] == 2
        assert set(result["variant_id"].to_list()) == {"v1", "v3"}

    def test_complex_expression_with_dotted_columns(self):
        """Test complex expressions with multiple dotted column names."""
        df = pl.DataFrame({
            "variant_id": ["v1", "v2", "v3", "v4"],
            "cadd_v1.7": [28.0, 30.0, 25.0, 35.0],
            "mpc.v2": [2.5, 3.0, 1.5, 3.5],
            "revel.score": [0.7, 0.8, 0.5, 0.9],
            "filters": ["PASS", "FAIL", "PASS", "PASS"],
        })

        expr = parse_impact_filter_expression(
            "(cadd_v1.7 >= 28.3 & mpc.v2 >= 2.6) | "
            "(revel.score >= 0.773 & filters = PASS)",
            df,
        )
        result = df.filter(expr)
        assert result.shape[0] == 2  # v2 and v4 match first condition
        assert set(result["variant_id"].to_list()) == {"v2", "v4"}

    def test_not_empty_combined_with_comparison(self):
        """Test combining not_empty with comparison operators on dotted columns."""
        df = pl.DataFrame({
            "variant_id": ["v1", "v2", "v3", "v4"],
            "pathogenicity.score": ["0.8", None, "0.95", "."],
        })

        expr = parse_impact_filter_expression(
            "pathogenicity.score not_empty & pathogenicity.score >= 0.9", df
        )
        result = df.filter(expr)
        assert result.shape[0] == 1
        assert result["variant_id"][0] == "v3"

    def test_error_on_nonexistent_dotted_column(self):
        """Test that proper error is raised for non-existent dotted columns."""
        df = pl.DataFrame({
            "CHROM": ["chr1"],
            "existing.column": [1.0],
        })

        # Should raise ValueError for non-existent column
        with pytest.raises(ValueError, match="not found in dataframe"):
            parse_impact_filter_expression("nonexistent.column >= 1.0", df)

    def test_single_dot_in_column_name(self):
        """Test column names with single dot."""
        df = pl.DataFrame({
            "score.1": [1.0, 2.0, 3.0],
            "score.2": [4.0, 5.0, 6.0],
        })

        expr = parse_impact_filter_expression("score.1 > 1.5 & score.2 < 5.5", df)
        result = df.filter(expr)
        assert result.shape[0] == 1

    def test_multiple_dots_in_column_name(self):
        """Test column names with multiple dots."""
        df = pl.DataFrame({
            "anno.db.v3.2.score": [10.0, 20.0, 30.0],
        })

        expr = parse_impact_filter_expression("anno.db.v3.2.score >= 20", df)
        result = df.filter(expr)
        assert result.shape[0] == 2
