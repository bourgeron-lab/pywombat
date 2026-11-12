"""Tests for pedigree reading, with and without a header line."""

import polars as pl
import pytest
from click.testing import CliRunner

from pywombat.cli import cli, read_pedigree

# A trio with the child listed first, as a pedigree subset per family often is
CHILD = "FAM1\tSample3\tSample1\tSample2\t1\t2"
FATHER = "FAM1\tSample1\t0\t0\t1\t1"
MOTHER = "FAM1\tSample2\t0\t0\t2\t1"


def _write(tmp_path, lines, name="ped.tsv"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def _by_sample(ped: pl.DataFrame) -> dict:
    return {row["sample_id"]: row for row in ped.iter_rows(named=True)}


class TestReadPedigree:
    def test_headerless_keeps_first_member(self, tmp_path):
        """A headerless pedigree's first line is a member, not column names."""
        ped = read_pedigree(_write(tmp_path, [CHILD, MOTHER, FATHER]))

        rows = _by_sample(ped)
        assert set(rows) == {"Sample1", "Sample2", "Sample3"}
        assert rows["Sample3"]["father_id"] == "Sample1"
        assert rows["Sample3"]["mother_id"] == "Sample2"
        assert str(rows["Sample3"]["sex"]) == "1"
        assert rows["Sample1"]["father_id"] is None
        assert rows["Sample2"]["mother_id"] is None

    @pytest.mark.parametrize(
        "header",
        [
            "FID\tsample_id\tFatherBarcode\tMotherBarcode\tSex\tPheno",
            "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE",
            "#FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE",
        ],
    )
    def test_header_is_not_a_member(self, tmp_path, header):
        ped = read_pedigree(_write(tmp_path, [header, CHILD, MOTHER, FATHER]))

        rows = _by_sample(ped)
        assert set(rows) == {"Sample1", "Sample2", "Sample3"}
        assert rows["Sample3"]["father_id"] == "Sample1"
        assert rows["Sample3"]["mother_id"] == "Sample2"
        assert str(rows["Sample3"]["sex"]) == "1"

    def test_headerless_numeric_ids_stay_strings(self, tmp_path):
        """IDs that look like numbers must still match the sample names as strings."""
        ped = read_pedigree(
            _write(
                tmp_path,
                [
                    "42\t00103\t00101\t00102\t2\t2",
                    "42\t00101\t0\t0\t1\t1",
                    "42\t00102\t0\t0\t2\t1",
                ],
            )
        )

        rows = _by_sample(ped)
        assert set(rows) == {"00101", "00102", "00103"}
        assert rows["00103"]["father_id"] == "00101"
        assert rows["00103"]["mother_id"] == "00102"

    def test_headerless_extra_columns(self, tmp_path):
        ped = read_pedigree(
            _write(tmp_path, [CHILD + "\textra", MOTHER + "\textra", FATHER + "\textra"])
        )

        assert _by_sample(ped)["Sample3"]["father_id"] == "Sample1"


class TestFilterWithHeaderlessPedigree:
    def test_first_listed_child_gets_parents(self, small_tsv, tmp_test_dir):
        """The child on a headerless pedigree's first line keeps its parents' genotypes."""
        ped = _write(tmp_test_dir, [CHILD, MOTHER, FATHER])
        output_prefix = tmp_test_dir / "output"

        result = CliRunner().invoke(
            cli, ["filter", str(small_tsv), "-o", str(output_prefix), "-p", str(ped)]
        )
        assert result.exit_code == 0, f"Command failed: {result.output}"

        df = pl.read_csv(
            tmp_test_dir / "output.tsv", separator="\t", infer_schema=False
        )
        child = df.filter(pl.col("sample") == "Sample3")
        assert child.height > 0
        assert child["father"].to_list() == ["Sample1"] * child.height
        assert child["mother"].to_list() == ["Sample2"] * child.height
        assert child["father_gt"].null_count() == 0
