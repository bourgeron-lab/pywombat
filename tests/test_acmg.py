"""Tests for ACMG variant-level evidence annotation."""

import polars as pl
import pytest
from pathlib import Path

from click.testing import CliRunner

from pywombat.cli import cli
from pywombat.acmg import (
    annotate_acmg,
    load_ruleset,
    _freq_status,
    _pp3_bp4_missense,
    _splice_evidence,
    _bp7,
    _pm4_bp3,
    _pvs1,
    _denovo_status,
    _phase_pair,
    _find_partners,
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


class TestPVS1:
    def test_frameshift_hc_no_flags(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["3/10"],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate:exon=3/10"

    def test_stop_gained_hc(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["stop_gained"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["5/12"],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate:exon=5/12"

    def test_splice_donor_hc(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["splice_donor_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate"

    def test_splice_acceptor_hc(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["splice_acceptor_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["2/8"],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate:exon=2/8"

    def test_hc_with_flags_is_review(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": ["SINGLE_EXON"],
            "VEP_EXON": ["1/1"],
        }).cast({"VEP_LoF_filter": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:lof_flags=SINGLE_EXON:exon=1/1"

    def test_lc_with_filter(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["stop_gained"],
            "VEP_LoF": ["LC"],
            "VEP_LoF_filter": ["END_TRUNC"],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["10/10"],
        }).cast({"VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:lof_lc(END_TRUNC):exon=10/10"

    def test_lc_no_filter(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_LoF": ["LC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:lof_lc"

    def test_os_loftee(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["splice_donor_variant"],
            "VEP_LoF": ["OS"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:lof_os"

    def test_start_lost_always_review(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["start_lost"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["1/10"],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:start_lost:exon=1/10"

    def test_no_loftee_data(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_EXON": ["4/10"],
        })
        result = _pvs1(df)
        assert result[0] == "review:no_loftee:exon=4/10"

    def test_loftee_dot_is_no_loftee(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["stop_gained"],
            "VEP_LoF": ["."],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "review:no_loftee"

    def test_missense_not_eligible(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] is None

    def test_synonymous_not_eligible(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["synonymous_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] is None

    def test_multi_consequence_with_null(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant&splice_region_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["3/10"],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate:exon=3/10"

    def test_no_exon_no_suffix(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate"

    def test_dot_exon_no_suffix(self):
        df = pl.DataFrame({
            "VEP_Consequence": ["stop_gained"],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": ["."],
            "VEP_EXON": ["."],
        }).cast({"VEP_LoF_filter": pl.Utf8})
        result = _pvs1(df)
        assert result[0] == "PVS1_candidate"

    def test_null_consequence(self):
        df = pl.DataFrame({
            "VEP_Consequence": [None],
            "VEP_LoF": ["HC"],
        }).cast({"VEP_Consequence": pl.Utf8})
        result = _pvs1(df)
        assert result[0] is None

    def test_no_consequence_column(self):
        df = pl.DataFrame({"other_col": [1]})
        result = _pvs1(df)
        assert result[0] is None


class TestPVS1InCriteria:
    def test_pvs1_candidate_in_criteria(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.0],
            "AC_joint": [0],
            "AN_joint": [30000],
            "REVEL": [None],
            "SpliceAI_maxDS": [None],
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
            "VEP_EXON": ["3/10"],
        }).cast({"REVEL": pl.Float64, "SpliceAI_maxDS": pl.Float64,
                  "VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8,
                  "VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert "PVS1_candidate" in criteria

    def test_pvs1_candidate_not_scored(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [None],
            "SpliceAI_maxDS": [None],
            "VEP_Consequence": ["frameshift_variant"],
            "VEP_HGVSc": [None],
            "LCR": [None],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"REVEL": pl.Float64, "SpliceAI_maxDS": pl.Float64,
                  "VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8,
                  "VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        assert result["acmg_points_auto"][0] is None

    def test_review_not_in_criteria(self, default_ruleset):
        df = pl.DataFrame({
            "fafmax_faf95_max_joint": [0.001],
            "AC_joint": [5],
            "AN_joint": [50000],
            "REVEL": [None],
            "SpliceAI_maxDS": [None],
            "VEP_Consequence": ["start_lost"],
            "VEP_HGVSc": [None],
            "LCR": [None],
            "VEP_LoF": ["HC"],
            "VEP_LoF_filter": [None],
            "VEP_LoF_flags": [None],
        }).cast({"REVEL": pl.Float64, "SpliceAI_maxDS": pl.Float64,
                  "VEP_HGVSc": pl.Utf8, "LCR": pl.Utf8,
                  "VEP_LoF_filter": pl.Utf8, "VEP_LoF_flags": pl.Utf8})
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert criteria is None or "PVS1" not in criteria

    def test_lof_missing_in_acmg_missing(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
        })
        result = annotate_acmg(df, default_ruleset)
        missing = result["acmg_missing"][0]
        assert "VEP_LoF" in missing
        assert "VEP_EXON" in missing


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

    def test_acmg_denovo_column_present(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
        })
        result = annotate_acmg(df, default_ruleset)
        assert "acmg_denovo" in result.columns


def _trio_row(chrom="chr1", pos=100, father_gt="0/0", mother_gt="0/0",
              father_dp="30", mother_dp="30", father_gq="99", mother_gq="99",
              father_vaf=0.0, mother_vaf=0.0, sex="1"):
    """Build a minimal melted row with parent columns for de novo tests."""
    return pl.DataFrame({
        "#CHROM": [chrom],
        "POS": [pos],
        "REF": ["A"],
        "ALT": ["T"],
        "sample": ["Child1"],
        "sample_gt": ["0/1"],
        "sample_dp": ["30"],
        "sample_gq": ["99"],
        "sample_ad": ["15"],
        "sample_vaf": [0.5],
        "father": ["Father1"],
        "mother": ["Mother1"],
        "father_gt": [father_gt],
        "father_dp": [father_dp],
        "father_gq": [father_gq],
        "father_ad": ["30"],
        "father_vaf": [father_vaf],
        "mother_gt": [mother_gt],
        "mother_dp": [mother_dp],
        "mother_gq": [mother_gq],
        "mother_ad": ["30"],
        "mother_vaf": [mother_vaf],
        "sex": [sex],
        "VEP_Consequence": ["missense_variant"],
    })


class TestDenovoStatus:
    def test_de_novo_autosome(self, default_ruleset):
        df = _trio_row(chrom="chr1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "de_novo"

    def test_inherited_paternal(self, default_ruleset):
        df = _trio_row(father_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "inherited_paternal"

    def test_inherited_maternal(self, default_ruleset):
        df = _trio_row(mother_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "inherited_maternal"

    def test_inherited_both(self, default_ruleset):
        df = _trio_row(father_gt="0/1", mother_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "inherited_both"

    def test_no_parent_columns(self, default_ruleset):
        df = pl.DataFrame({
            "#CHROM": ["chr1"],
            "POS": [100],
            "VEP_Consequence": ["missense_variant"],
        })
        result = _denovo_status(df, default_ruleset)
        assert result[0] is None

    def test_no_parents_in_pedigree(self, default_ruleset):
        df = _trio_row()
        df = df.with_columns([
            pl.lit(None).cast(pl.Utf8).alias("father"),
            pl.lit(None).cast(pl.Utf8).alias("mother"),
        ])
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "not_testable:no_parents"

    def test_parent_low_dp(self, default_ruleset):
        df = _trio_row(father_dp="5")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "not_testable:parent_qc"

    def test_parent_low_gq(self, default_ruleset):
        df = _trio_row(mother_gq="10")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "not_testable:parent_qc"

    def test_parent_missing_gt(self, default_ruleset):
        df = _trio_row(father_gt="./.")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "not_testable:parent_gt"

    def test_parent_high_vaf(self, default_ruleset):
        df = _trio_row(father_vaf=0.05)
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "not_testable:parent_vaf"

    def test_x_linked_male_de_novo(self, default_ruleset):
        df = _trio_row(chrom="chrX", pos=50000000, sex="1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "de_novo"

    def test_x_linked_male_inherited_maternal(self, default_ruleset):
        df = _trio_row(chrom="chrX", pos=50000000, sex="1", mother_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "inherited_maternal"

    def test_x_linked_male_father_ignored(self, default_ruleset):
        """Father's X genotype is not informative for males outside PAR."""
        df = _trio_row(chrom="chrX", pos=50000000, sex="1", father_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "de_novo"

    def test_x_par_tests_both_parents(self, default_ruleset):
        dnm_cfg = {
            "par_regions": {
                "grch38": {
                    "PAR1": {"chrom": "X", "start": 10000, "end": 2781479},
                }
            }
        }
        df = _trio_row(chrom="chrX", pos=100000, sex="1", father_gt="0/1")
        result = _denovo_status(df, default_ruleset, dnm_cfg)
        assert result[0] == "inherited_paternal"

    def test_y_de_novo(self, default_ruleset):
        df = _trio_row(chrom="chrY", pos=10000000, sex="1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "de_novo"

    def test_y_inherited_paternal(self, default_ruleset):
        df = _trio_row(chrom="chrY", pos=10000000, sex="1", father_gt="0/1")
        result = _denovo_status(df, default_ruleset)
        assert result[0] == "inherited_paternal"

    def test_y_excluded_region(self, default_ruleset):
        dnm_cfg = {
            "y_excluded_regions": {
                "grch38": {
                    "Y_Ampliconic1": {"chrom": "Y", "start": 6235111, "end": 6532906},
                }
            }
        }
        df = _trio_row(chrom="chrY", pos=6300000, sex="1")
        result = _denovo_status(df, default_ruleset, dnm_cfg)
        assert result[0] == "not_testable:y_excluded"

    def test_de_novo_in_criteria(self, default_ruleset):
        df = _trio_row()
        df = df.with_columns([
            pl.lit(0.001).alias("fafmax_faf95_max_joint"),
            pl.lit(5).alias("AC_joint"),
            pl.lit(50000).alias("AN_joint"),
            pl.lit(0.80).alias("REVEL"),
            pl.lit(0.05).alias("SpliceAI_maxDS"),
            pl.lit(None).cast(pl.Utf8).alias("VEP_HGVSc"),
            pl.lit(None).cast(pl.Utf8).alias("LCR"),
        ])
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert "de_novo" in criteria

    def test_inherited_not_in_criteria(self, default_ruleset):
        df = _trio_row(father_gt="0/1")
        df = df.with_columns([
            pl.lit(0.001).alias("fafmax_faf95_max_joint"),
            pl.lit(5).alias("AC_joint"),
            pl.lit(50000).alias("AN_joint"),
            pl.lit(0.80).alias("REVEL"),
            pl.lit(0.05).alias("SpliceAI_maxDS"),
            pl.lit(None).cast(pl.Utf8).alias("VEP_HGVSc"),
            pl.lit(None).cast(pl.Utf8).alias("LCR"),
        ])
        result = annotate_acmg(df, default_ruleset)
        criteria = result["acmg_criteria"][0]
        assert "de_novo" not in (criteria or "")

    def test_dnm_config_overrides_ruleset(self, default_ruleset):
        dnm_cfg = {"parent_vaf_max": 0.001}
        df = _trio_row(father_vaf=0.015)
        result = _denovo_status(df, default_ruleset, dnm_cfg)
        assert result[0] == "not_testable:parent_vaf"
        result_default = _denovo_status(df, default_ruleset)
        assert result_default[0] == "de_novo"


class TestPhasePair:
    def test_paternal_maternal_is_trans(self):
        assert _phase_pair("inherited_paternal", "inherited_maternal") == "trans"

    def test_maternal_paternal_is_trans(self):
        assert _phase_pair("inherited_maternal", "inherited_paternal") == "trans"

    def test_paternal_paternal_is_cis(self):
        assert _phase_pair("inherited_paternal", "inherited_paternal") == "cis"

    def test_maternal_maternal_is_cis(self):
        assert _phase_pair("inherited_maternal", "inherited_maternal") == "cis"

    def test_de_novo_inherited_paternal_is_trans(self):
        assert _phase_pair("de_novo", "inherited_paternal") == "trans"

    def test_inherited_maternal_de_novo_is_trans(self):
        assert _phase_pair("inherited_maternal", "de_novo") == "trans"

    def test_de_novo_de_novo_is_unknown(self):
        assert _phase_pair("de_novo", "de_novo") == "unknown"

    def test_inherited_both_is_unknown(self):
        assert _phase_pair("inherited_both", "inherited_paternal") == "unknown"
        assert _phase_pair("inherited_maternal", "inherited_both") == "unknown"

    def test_not_testable_is_unknown(self):
        assert _phase_pair("not_testable:parent_qc", "inherited_paternal") == "unknown"

    def test_none_is_unknown(self):
        assert _phase_pair(None, "inherited_paternal") == "unknown"
        assert _phase_pair("de_novo", None) == "unknown"
        assert _phase_pair(None, None) == "unknown"

    def test_de_novo_inherited_both_is_unknown(self):
        assert _phase_pair("de_novo", "inherited_both") == "unknown"


def _partner_df(rows):
    """Build a DataFrame from a list of dicts for partner tests."""
    if not rows:
        return pl.DataFrame()
    all_keys = set()
    for r in rows:
        all_keys.update(r.keys())
    data = {k: [r.get(k) for r in rows] for k in all_keys}
    return pl.DataFrame(data)


class TestFindPartners:
    def test_two_variants_same_gene_trans(self, default_ruleset):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] == "chr1:200:G:C(trans)"
        assert result[1] == "chr1:100:A:T(trans)"

    def test_two_variants_same_gene_cis(self, default_ruleset):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] == "chr1:200:G:C(cis)"
        assert result[1] == "chr1:100:A:T(cis)"

    def test_singleton_no_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "de_novo", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] is None

    def test_different_genes_no_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "TP53", "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] is None

    def test_different_samples_no_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child2",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] is None

    def test_ba1_excluded_as_partner(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "BA1"},
        ])
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] == "chr1:100:A:T(trans)"

    def test_no_vep_symbol_no_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": None, "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": None, "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] is None

    def test_dot_symbol_no_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": ".", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": ".", "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] is None

    def test_three_variants_multiple_partners(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_maternal", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "300", "REF": "C", "ALT": "G",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "de_novo", "acmg_freq": "absent"},
        ])
        result = _find_partners(df)
        parts_0 = result[0].split(";")
        assert len(parts_0) == 2
        assert "chr1:200:G:C(trans)" in parts_0
        assert "chr1:300:C:G(trans)" in parts_0

    def test_no_denovo_column_gives_unknown_phase(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1", "acmg_freq": "rare"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1", "acmg_freq": "rare"},
        ])
        result = _find_partners(df)
        assert result[0] == "chr1:200:G:C(unknown)"

    def test_missing_required_columns(self):
        df = pl.DataFrame({"other_col": [1, 2]})
        result = _find_partners(df)
        assert result[0] is None
        assert result[1] is None

    def test_no_freq_column_still_works(self):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_paternal"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "acmg_denovo": "inherited_maternal"},
        ])
        result = _find_partners(df)
        assert result[0] == "chr1:200:G:C(trans)"


class TestAnnotateACMGPartners:
    def test_partners_column_present(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
        })
        result = annotate_acmg(df, default_ruleset)
        assert "acmg_partners" in result.columns

    def test_partners_in_full_annotation(self, default_ruleset):
        df = _partner_df([
            {"#CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "fafmax_faf95_max_joint": 0.001, "AC_joint": 5, "AN_joint": 50000,
             "REVEL": 0.80, "SpliceAI_maxDS": 0.05,
             "VEP_Consequence": "missense_variant", "VEP_HGVSc": None,
             "LCR": None,
             "father": "Father1", "mother": "Mother1",
             "father_gt": "0/1", "father_dp": "30", "father_gq": "99",
             "father_ad": "30", "father_vaf": 0.5,
             "mother_gt": "0/0", "mother_dp": "30", "mother_gq": "99",
             "mother_ad": "30", "mother_vaf": 0.0,
             "sex": "1"},
            {"#CHROM": "chr1", "POS": "200", "REF": "G", "ALT": "C",
             "VEP_SYMBOL": "BRCA1", "sample": "Child1",
             "fafmax_faf95_max_joint": 0.001, "AC_joint": 5, "AN_joint": 50000,
             "REVEL": 0.65, "SpliceAI_maxDS": 0.05,
             "VEP_Consequence": "missense_variant", "VEP_HGVSc": None,
             "LCR": None,
             "father": "Father1", "mother": "Mother1",
             "father_gt": "0/0", "father_dp": "30", "father_gq": "99",
             "father_ad": "30", "father_vaf": 0.0,
             "mother_gt": "0/1", "mother_dp": "30", "mother_gq": "99",
             "mother_ad": "30", "mother_vaf": 0.5,
             "sex": "1"},
        ])
        result = annotate_acmg(df, default_ruleset)
        assert result["acmg_partners"][0] == "chr1:200:G:C(trans)"
        assert result["acmg_partners"][1] == "chr1:100:A:T(trans)"

    def test_vep_symbol_in_missing_when_absent(self, default_ruleset):
        df = pl.DataFrame({
            "VEP_Consequence": ["missense_variant"],
        })
        result = annotate_acmg(df, default_ruleset)
        assert "VEP_SYMBOL" in result["acmg_missing"][0]


def _write_filtered_tsv(path, rows=None):
    """Write a minimal filtered-output TSV (melted format) for acmg command tests."""
    if rows is None:
        rows = [
            "chr1\t100\tA\tT\t.\tSample1\t0/1\t20\t99\t5\t0.25\tmissense_variant\tGENE1\tYES\t0.75\t0.3",
        ]
    header = (
        "#CHROM\tPOS\tREF\tALT\tFILTER\tsample\tsample_gt\tsample_dp\tsample_gq\t"
        "sample_ad\tsample_vaf\tVEP_Consequence\tVEP_SYMBOL\tVEP_CANONICAL\tREVEL\tSpliceAI_maxDS"
    )
    path.write_text("\n".join([header] + rows) + "\n")
    return path


class TestAcmgCommand:
    def test_basic_annotation(self, tmp_path):
        tsv = _write_filtered_tsv(tmp_path / "input.tsv")
        out = str(tmp_path / "output")
        result = CliRunner().invoke(cli, ["acmg", str(tsv), "-o", out])
        assert result.exit_code == 0, f"Command failed: {result.output}"

        df = pl.read_csv(tmp_path / "output.tsv", separator="\t", infer_schema=False)
        assert "acmg_ruleset" in df.columns
        assert "acmg_freq" in df.columns
        assert "acmg_criteria" in df.columns
        assert df["acmg_ruleset"][0].startswith("acmg-svi-2024.1/pywombat-")

    def test_strips_existing_acmg_columns(self, tmp_path):
        header = (
            "#CHROM\tPOS\tREF\tALT\tFILTER\tsample\tsample_gt\tsample_dp\tsample_gq\t"
            "sample_ad\tsample_vaf\tVEP_Consequence\tVEP_SYMBOL\tVEP_CANONICAL\tREVEL\t"
            "SpliceAI_maxDS\tacmg_ruleset\tacmg_freq\tacmg_criteria"
        )
        row = (
            "chr1\t100\tA\tT\t.\tSample1\t0/1\t20\t99\t5\t0.25\tmissense_variant\t"
            "GENE1\tYES\t0.75\t0.3\told-ruleset/pywombat-1.0.0\trare\tPP3_Moderate"
        )
        tsv = tmp_path / "input.tsv"
        tsv.write_text(header + "\n" + row + "\n")

        out = str(tmp_path / "output")
        result = CliRunner().invoke(cli, ["acmg", str(tsv), "-o", out, "-v"])
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Stripping 3 existing acmg_* columns" in result.output

        df = pl.read_csv(tmp_path / "output.tsv", separator="\t", infer_schema=False)
        assert df["acmg_ruleset"][0].startswith("acmg-svi-2024.1/pywombat-")

    def test_parquet_input(self, tmp_path):
        tsv = _write_filtered_tsv(tmp_path / "input.tsv")
        df = pl.read_csv(tsv, separator="\t")
        pq = tmp_path / "input.parquet"
        df.write_parquet(pq)

        out = str(tmp_path / "output")
        result = CliRunner().invoke(cli, ["acmg", str(pq), "-o", out])
        assert result.exit_code == 0, f"Command failed: {result.output}"

        result_df = pl.read_csv(tmp_path / "output.tsv", separator="\t", infer_schema=False)
        assert "acmg_ruleset" in result_df.columns

    def test_default_output_name(self, tmp_path, monkeypatch):
        tsv = _write_filtered_tsv(tmp_path / "mydata.tsv")
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(cli, ["acmg", str(tsv)])
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert (tmp_path / "mydata.acmg.tsv").exists()

    def test_output_format_parquet(self, tmp_path):
        tsv = _write_filtered_tsv(tmp_path / "input.tsv")
        out = str(tmp_path / "output")
        result = CliRunner().invoke(cli, ["acmg", str(tsv), "-o", out, "-f", "parquet"])
        assert result.exit_code == 0, f"Command failed: {result.output}"

        result_df = pl.read_parquet(tmp_path / "output.parquet")
        assert "acmg_ruleset" in result_df.columns

    def test_verbose_output(self, tmp_path):
        tsv = _write_filtered_tsv(tmp_path / "input.tsv")
        out = str(tmp_path / "output")
        result = CliRunner().invoke(cli, ["acmg", str(tsv), "-o", out, "-v"])
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Reading input file" in result.output
        assert "Annotating ACMG evidence" in result.output

    def test_acmg_help(self):
        result = CliRunner().invoke(cli, ["acmg", "--help"])
        assert result.exit_code == 0
        assert "Re-annotate ACMG evidence" in result.output
        assert "--ruleset" in result.output
        assert "--pedigree" in result.output

    def test_with_dnm_config(self, tmp_path):
        rows = [
            "chr1\t100\tA\tT\t.\tChild1\t0/1\t20\t99\t5\t0.25\tmissense_variant\tGENE1\tYES\t0.75\t0.3",
        ]
        tsv = _write_filtered_tsv(tmp_path / "input.tsv", rows)

        config = tmp_path / "config.yml"
        config.write_text(
            "dnm:\n"
            "  enabled: true\n"
            "  parent_dp_min: 10\n"
            "  parent_gq_min: 18\n"
            "  parent_vaf_max: 0.02\n"
        )

        out = str(tmp_path / "output")
        result = CliRunner().invoke(
            cli, ["acmg", str(tsv), "-o", out, "--dnm-config", str(config)]
        )
        assert result.exit_code == 0, f"Command failed: {result.output}"
        df = pl.read_csv(tmp_path / "output.tsv", separator="\t", infer_schema=False)
        assert "acmg_denovo" in df.columns
