"""Tests for ACMG variant-level evidence annotation."""

import polars as pl
import pytest
from pathlib import Path

from pywombat.acmg import (
    annotate_acmg,
    load_ruleset,
    _freq_status,
    _pp3_bp4_missense,
    _splice_evidence,
    _bp7,
    _pm4_bp3,
    _parse_intron_offset,
    _has_consequence,
)


@pytest.fixture
def default_ruleset():
    return load_ruleset("acmg-svi-2024.1")


class TestLoadRuleset:
    def test_load_builtin(self):
        rs = load_ruleset("acmg-svi-2024.1")
        assert rs["name"] == "acmg-svi-2024.1"
        assert rs["framework"] == "points"
        assert rs["points"]["supporting"] == 1

    def test_load_by_path(self):
        path = Path(__file__).parent.parent / "src" / "pywombat" / "rulesets" / "acmg-svi-2024.1.yml"
        rs = load_ruleset(path)
        assert rs["name"] == "acmg-svi-2024.1"

    def test_missing_raises(self):
        with pytest.raises(FileNotFoundError):
            load_ruleset("nonexistent-ruleset")


class TestFrequency:
    def test_ba1(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.06],
            "AC_joint": [100],
            "AN_joint": [50000],
        })
        result = _freq_status(df, default_ruleset)
        assert result[0] == "BA1"

    def test_ba1_boundary(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.05],
            "AC_joint": [100],
            "AN_joint": [50000],
        })
        result = _freq_status(df, default_ruleset)
        assert result[0] == "rare"

    def test_absent_high_coverage(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.0],
            "AC_joint": [0],
            "AN_joint": [30000],
        })
        result = _freq_status(df, default_ruleset)
        assert result[0] == "absent"

    def test_absent_low_coverage(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.0],
            "AC_joint": [0],
            "AN_joint": [5000],
        })
        result = _freq_status(df, default_ruleset)
        assert result[0] == "absent_low_coverage"

    def test_null_freq_absent_low_coverage(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [None],
            "AC_joint": [None],
            "AN_joint": [None],
        }).cast({"fafmax_faf95_max_joint": pl.Float64})
        result = _freq_status(df, default_ruleset)
        assert result[0] == "absent_low_coverage"

    def test_rare(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
        })
        result = _freq_status(df, default_ruleset)
        assert result[0] == "rare"

    def test_missing_column(self, default_ruleset):
        df = pl.DataFrame({"other_col": [1]})
        result = _freq_status(df, default_ruleset)
        assert result[0] is None


class TestPP3BP4Missense:
    def test_pp3_strong(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.95],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "PP3_Strong" in result[0]

    def test_pp3_moderate(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.80],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "PP3_Moderate" in result[0]

    def test_pp3_supporting(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.65],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "PP3_Supporting" in result[0]

    def test_bp4_supporting(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.28],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "BP4_Supporting" in result[0]

    def test_bp4_moderate(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.15],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "BP4_Moderate" in result[0]

    def test_bp4_strong(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.01],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "BP4_Strong" in result[0]

    def test_bp4_very_strong(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.002],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "BP4_VeryStrong" in result[0]

    def test_non_missense_returns_none(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.95],
            "VEP_Consequence": ["synonymous_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert result[0] is None

    def test_null_score_returns_none(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [None],
            "VEP_Consequence": ["missense_variant"],
        }).cast({"REVEL": pl.Float64})
        result = _pp3_bp4_missense(df, default_ruleset)
        assert result[0] is None

    def test_indeterminate_score(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.5],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert result[0] is None

    def test_boundary_pp3_supporting(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.644],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "PP3_Supporting" in result[0]

    def test_boundary_bp4_supporting(self, default_ruleset):
        df = pl.DataFrame({
            "REVEL": [0.290],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _pp3_bp4_missense(df, default_ruleset)
        assert "BP4_Supporting" in result[0]


class TestSpliceEvidence:
    def test_pp3_splice(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.5],
            "VEP_Consequence": ["intron_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert "PP3_Supporting:SpliceAI" in result[0]

    def test_bp4_splice(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.05],
            "VEP_Consequence": ["intron_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert "BP4_Supporting:SpliceAI" in result[0]

    def test_canonical_splice_site_excluded(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.9],
            "VEP_Consequence": ["splice_donor_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert result[0] is None

    def test_canonical_acceptor_excluded(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.9],
            "VEP_Consequence": ["splice_acceptor_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert result[0] is None

    def test_null_score(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [None],
            "VEP_Consequence": ["intron_variant"],
        }).cast({"SpliceAI_maxDS": pl.Float64})
        result = _splice_evidence(df, default_ruleset)
        assert result[0] is None

    def test_indeterminate_score(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.15],
            "VEP_Consequence": ["intron_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert result[0] is None

    def test_boundary_pp3(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.2],
            "VEP_Consequence": ["intron_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert "PP3_Supporting:SpliceAI" in result[0]

    def test_boundary_bp4(self, default_ruleset):
        df = pl.DataFrame({
            "SpliceAI_maxDS": [0.1],
            "VEP_Consequence": ["intron_variant"],
        })
        result = _splice_evidence(df, default_ruleset)
        assert "BP4_Supporting:SpliceAI" in result[0]


class TestBP7:
    def test_synonymous_low_splice(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["synonymous_variant"],
            "SpliceAI_maxDS": [0.05],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] == "BP7"

    def test_synonymous_high_splice_no_bp7(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["synonymous_variant"],
            "SpliceAI_maxDS": [0.3],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] is None

    def test_deep_intronic_bp7(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["intron_variant"],
            "SpliceAI_maxDS": [0.05],
            "VEP_HGVSc": ["ENST00000454634.6:c.123+8A>G"],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] == "BP7"

    def test_near_splice_intronic_no_bp7(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["intron_variant"],
            "SpliceAI_maxDS": [0.05],
            "VEP_HGVSc": ["ENST00000454634.6:c.123+3A>G"],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] is None

    def test_acceptor_deep(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["intron_variant"],
            "SpliceAI_maxDS": [0.05],
            "VEP_HGVSc": ["ENST00000454634.6:c.124-25A>G"],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] == "BP7"

    def test_acceptor_near_splice_no_bp7(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["intron_variant"],
            "SpliceAI_maxDS": [0.05],
            "VEP_HGVSc": ["ENST00000454634.6:c.124-10A>G"],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] is None

    def test_missense_not_eligible(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
            "SpliceAI_maxDS": [0.05],
        })
        result = _bp7(df, default_ruleset)
        assert result[0] is None


class TestPM4BP3:
    def test_inframe_insertion(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["inframe_insertion"],
            "LCR": [None],
        }).cast({"LCR": pl.Utf8})
        result = _pm4_bp3(df)
        assert result[0] == "PM4_candidate"

    def test_inframe_deletion(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["inframe_deletion"],
            "LCR": [None],
        }).cast({"LCR": pl.Utf8})
        result = _pm4_bp3(df)
        assert result[0] == "PM4_candidate"

    def test_stop_lost(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["stop_lost"],
            "LCR": [None],
        }).cast({"LCR": pl.Utf8})
        result = _pm4_bp3(df)
        assert result[0] == "PM4_candidate"

    def test_inframe_in_lcr_is_bp3(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["inframe_insertion"],
            "LCR": ["true"],
        })
        result = _pm4_bp3(df)
        assert result[0] == "BP3_candidate"

    def test_missense_not_eligible(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
            "LCR": [None],
        }).cast({"LCR": pl.Utf8})
        result = _pm4_bp3(df)
        assert result[0] is None

    def test_no_lcr_column(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["inframe_insertion"],
        })
        result = _pm4_bp3(df)
        assert result[0] == "PM4_candidate"


class TestParseIntronOffset:
    def test_donor_positive(self):
        assert _parse_intron_offset("ENST00000454634.6:c.123+7A>G") == 7

    def test_acceptor_negative(self):
        assert _parse_intron_offset("ENST00000454634.6:c.124-21A>G") == -21

    def test_no_offset(self):
        assert _parse_intron_offset("ENST00000454634.6:c.123A>G") is None

    def test_none(self):
        assert _parse_intron_offset(None) is None


class TestHasConsequence:
    def test_single_match(self):
        assert _has_consequence("missense_variant", {"missense_variant"})

    def test_ampersand_separated(self):
        assert _has_consequence("missense_variant&splice_region_variant", {"missense_variant"})

    def test_no_match(self):
        assert not _has_consequence("synonymous_variant", {"missense_variant"})


class TestPP3OnceRule:
    def test_missense_wins_over_splice(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [0.95],
            "SpliceAI_maxDS": [0.5],
            "VEP_Consequence": ["missense_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
        }).cast({"VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert criteria.count("PP3") == 1
        assert "PP3_Strong" in criteria

    def test_splice_wins_when_no_missense_pp3(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [0.5],
            "SpliceAI_maxDS": [0.5],
            "VEP_Consequence": ["missense_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
        }).cast({"VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert "PP3_Supporting" in criteria


class TestAnnotateACMG:
    def test_full_annotation(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [0.80],
            "SpliceAI_maxDS": [0.05],
            "VEP_Consequence": ["missense_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
        }).cast({"VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8})

        result = annotate_acmg(df, default_ruleset)

        assert "acmg_ruleset" in result.columns
        assert "acmg_freq" in result.columns
        assert "acmg_pp3_bp4" in result.columns
        assert "acmg_splice" in result.columns
        assert "acmg_bp7" in result.columns
        assert "acmg_pm4_bp3" in result.columns
        assert "acmg_criteria" in result.columns
        assert "acmg_points_auto" in result.columns
        assert "acmg_missing" in result.columns
        assert result["acmg_ruleset"][0].startswith("acmg-svi-2024.1/pywombat-")
        assert result["acmg_freq"][0] == "rare"
        assert "PP3_Moderate" in result["acmg_pp3_bp4"][0]

    def test_ba1_gives_negative_points(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.1],
            "AC_joint": [1000],
            "AN_joint": [50000],
            "REVEL": [None],
            "SpliceAI_maxDS": [None],
            "VEP_Consequence": ["missense_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
        }).cast({"REVEL": pl.Float64, "SpliceAI_maxDS": pl.Float64, "VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8})

        result = annotate_acmg(df, default_ruleset)
        assert result["acmg_points_auto"][0] == -8
        assert "BA1" in result["acmg_criteria"][0]

    def test_missing_columns_reported(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
        })
        result = annotate_acmg(df, default_ruleset)
        missing = result["acmg_missing"][0]
        assert "REVEL" in missing
        assert "SpliceAI_maxDS" in missing

    def test_no_criteria_gives_none(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [0.5],
            "SpliceAI_maxDS": [0.15],
            "VEP_Consequence": ["missense_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
        }).cast({"VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        assert result["acmg_criteria"][0] is None
        assert result["acmg_points_auto"][0] is None
