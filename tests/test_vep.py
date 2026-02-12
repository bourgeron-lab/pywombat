"""Tests for VEP annotation module."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest
from click.testing import CliRunner

from pywombat.cli import cli, _normalize_annotate_mode
from pywombat.vep import (
    VEP_COLUMN_MAPPING,
    VEP_FILE_COLUMN_MAPPING,
    _find_transcript_annotation,
    _format_variant_for_vep,
    annotate_mnv_variants,
    create_mnv_placeholders_and_vcf,
    extract_annotations,
    format_annotations,
    merge_vep_annotations,
    parse_variant,
    parse_vep_annotation_file,
    query_vep,
)


# --- Fixtures ---


@pytest.fixture
def sample_vep_response():
    """Return a realistic VEP API response dict with multiple transcripts."""
    return {
        "input": "1 151776102 . GC AA . . .",
        "most_severe_consequence": "missense_variant",
        "assembly_name": "GRCh38",
        "transcript_consequences": [
            {
                "gene_symbol": "GENE1",
                "gene_id": "ENSG00000000001",
                "transcript_id": "ENST00000000001",
                "consequence_terms": ["missense_variant"],
                "impact": "MODERATE",
                "canonical": 1,
                "biotype": "protein_coding",
                "hgvsc": "ENST00000000001.3:c.100G>A",
                "hgvsp": "ENSP00000000001:p.Ala34Thr",
                "protein_position": "34/500",
                "amino_acids": "A/T",
                "codons": "Gcc/Acc",
                "exon": "3/12",
                "mane_select": "NM_001234.5",
                "sift_prediction": "deleterious",
                "sift_score": 0.01,
                "polyphen_prediction": "probably_damaging",
                "polyphen_score": 0.998,
            },
            {
                "gene_symbol": "GENE1",
                "gene_id": "ENSG00000000001",
                "transcript_id": "ENST00000000002",
                "consequence_terms": ["downstream_gene_variant"],
                "impact": "MODIFIER",
                "biotype": "retained_intron",
            },
            {
                "gene_symbol": "GENE1",
                "gene_id": "ENSG00000000001",
                "transcript_id": "ENST00000000003",
                "consequence_terms": ["missense_variant", "splice_region_variant"],
                "impact": "MODERATE",
                "biotype": "protein_coding",
                "hgvsc": "ENST00000000003.1:c.50G>A",
            },
        ],
    }


@pytest.fixture
def lof_vep_response():
    """Return a VEP response with LOFTEE annotations."""
    return {
        "input": "17 41245466 . G A . . .",
        "most_severe_consequence": "stop_gained",
        "transcript_consequences": [
            {
                "gene_symbol": "BRCA1",
                "gene_id": "ENSG00000012048",
                "transcript_id": "ENST00000357654",
                "consequence_terms": ["stop_gained"],
                "impact": "HIGH",
                "canonical": 1,
                "biotype": "protein_coding",
                "lof": "HC",
                "lof_filter": None,
                "lof_flags": "SINGLE_EXON",
                "lof_info": "PERCENTILE:0.05",
            },
        ],
    }


@pytest.fixture
def intergenic_vep_response():
    """Return a VEP response for an intergenic variant."""
    return {
        "input": "1 100000 . A T . . .",
        "most_severe_consequence": "intergenic_variant",
        "intergenic_consequences": [
            {
                "consequence_terms": ["intergenic_variant"],
                "impact": "MODIFIER",
            }
        ],
    }


# --- TestParseVariant ---


class TestParseVariant:
    """Test variant string parsing."""

    def test_standard_format(self):
        result = parse_variant("chr1:151776102:GC:AA")
        assert result["chrom"] == "1"
        assert result["pos"] == 151776102
        assert result["ref"] == "GC"
        assert result["alt"] == "AA"
        assert result["chrom_raw"] == "chr1"

    def test_without_chr_prefix(self):
        result = parse_variant("1:151776102:GC:AA")
        assert result["chrom"] == "1"
        assert result["chrom_raw"] == "1"

    def test_chrx(self):
        result = parse_variant("chrX:12345:A:T")
        assert result["chrom"] == "X"

    def test_chry(self):
        result = parse_variant("chrY:12345:A:T")
        assert result["chrom"] == "Y"

    def test_snv(self):
        result = parse_variant("chr1:100:A:T")
        assert result["ref"] == "A"
        assert result["alt"] == "T"

    def test_lowercase_normalized(self):
        result = parse_variant("chr1:100:a:t")
        assert result["ref"] == "A"
        assert result["alt"] == "T"

    def test_invalid_too_few_fields(self):
        with pytest.raises(ValueError, match="Invalid variant format"):
            parse_variant("chr1:100:A")

    def test_invalid_too_many_fields(self):
        with pytest.raises(ValueError, match="Invalid variant format"):
            parse_variant("chr1:100:A:T:extra")

    def test_invalid_position(self):
        with pytest.raises(ValueError, match="Invalid position"):
            parse_variant("chr1:abc:A:T")

    def test_negative_position(self):
        with pytest.raises(ValueError, match="must be positive"):
            parse_variant("chr1:-5:A:T")

    def test_invalid_ref_bases(self):
        with pytest.raises(ValueError, match="Invalid REF"):
            parse_variant("chr1:100:Z:T")

    def test_invalid_alt_bases(self):
        with pytest.raises(ValueError, match="Invalid ALT"):
            parse_variant("chr1:100:A:Z")


# --- TestFormatVariantForVep ---


class TestFormatVariantForVep:
    """Test VEP API format conversion."""

    def test_snv(self):
        variant = {"chrom": "1", "pos": 100, "ref": "A", "alt": "T"}
        assert _format_variant_for_vep(variant) == "1 100 . A T . . ."

    def test_mnv(self):
        variant = {"chrom": "1", "pos": 151776102, "ref": "GC", "alt": "AA"}
        assert _format_variant_for_vep(variant) == "1 151776102 . GC AA . . ."

    def test_chrx(self):
        variant = {"chrom": "X", "pos": 12345, "ref": "A", "alt": "G"}
        assert _format_variant_for_vep(variant) == "X 12345 . A G . . ."


# --- TestExtractAnnotations ---


class TestExtractAnnotations:
    """Test VEP response annotation extraction."""

    def test_canonical_only(self, sample_vep_response):
        result = extract_annotations(sample_vep_response, canonical_only=True)
        assert len(result) == 1
        assert result[0]["transcript_id"] == "ENST00000000001"
        assert result[0]["canonical"] is True

    def test_all_transcripts(self, sample_vep_response):
        result = extract_annotations(sample_vep_response, canonical_only=False)
        assert len(result) == 3

    def test_fields_extracted(self, sample_vep_response):
        result = extract_annotations(sample_vep_response, canonical_only=True)
        ann = result[0]
        assert ann["gene_symbol"] == "GENE1"
        assert ann["consequence"] == "missense_variant"
        assert ann["impact"] == "MODERATE"
        assert ann["hgvsc"] == "ENST00000000001.3:c.100G>A"
        assert ann["hgvsp"] == "ENSP00000000001:p.Ala34Thr"
        assert ann["amino_acids"] == "A/T"
        assert ann["exon"] == "3/12"
        assert ann["mane_select"] == "NM_001234.5"
        assert ann["sift_prediction"] == "deleterious"
        assert ann["sift_score"] == 0.01
        assert ann["polyphen_prediction"] == "probably_damaging"
        assert ann["polyphen_score"] == 0.998

    def test_consequence_terms_joined(self, sample_vep_response):
        result = extract_annotations(sample_vep_response, canonical_only=False)
        # Third transcript has two consequence terms
        third = [a for a in result if a["transcript_id"] == "ENST00000000003"][0]
        assert third["consequence"] == "missense_variant&splice_region_variant"

    def test_missing_fields_are_none(self, sample_vep_response):
        result = extract_annotations(sample_vep_response, canonical_only=False)
        # Second transcript has minimal fields
        second = [a for a in result if a["transcript_id"] == "ENST00000000002"][0]
        assert second["hgvsc"] is None
        assert second["hgvsp"] is None
        assert second["lof"] is None
        assert second["sift_prediction"] is None

    def test_lof_fields_extracted(self, lof_vep_response):
        result = extract_annotations(lof_vep_response, canonical_only=True)
        assert len(result) == 1
        ann = result[0]
        assert ann["lof"] == "HC"
        assert ann["lof_filter"] is None
        assert ann["lof_flags"] == "SINGLE_EXON"
        assert ann["lof_info"] == "PERCENTILE:0.05"

    def test_intergenic_returns_empty(self, intergenic_vep_response):
        result = extract_annotations(intergenic_vep_response, canonical_only=True)
        assert result == []

    def test_no_canonical_returns_empty(self):
        response = {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST00000000099",
                    "consequence_terms": ["intron_variant"],
                    "impact": "MODIFIER",
                }
            ]
        }
        result = extract_annotations(response, canonical_only=True)
        assert result == []


# --- TestFormatAnnotations ---


class TestFormatAnnotations:
    """Test CLI output formatting."""

    def test_header_present(self, sample_vep_response):
        annotations = extract_annotations(sample_vep_response, canonical_only=True)
        output = format_annotations(annotations, "chr1:151776102:GC:AA")
        assert "Variant: chr1:151776102:GC:AA" in output

    def test_canonical_label(self, sample_vep_response):
        annotations = extract_annotations(sample_vep_response, canonical_only=True)
        output = format_annotations(annotations, "chr1:151776102:GC:AA")
        assert "(canonical)" in output
        assert "ENST00000000001" in output

    def test_fields_displayed(self, sample_vep_response):
        annotations = extract_annotations(sample_vep_response, canonical_only=True)
        output = format_annotations(annotations, "chr1:151776102:GC:AA")
        assert "GENE1" in output
        assert "missense_variant" in output
        assert "MODERATE" in output
        assert "SIFT" in output
        assert "PolyPhen" in output

    def test_none_fields_omitted(self):
        annotations = [
            {
                "gene_symbol": "TEST",
                "gene_id": "ENSG000",
                "transcript_id": "ENST000",
                "consequence": "intron_variant",
                "impact": "MODIFIER",
                "biotype": "protein_coding",
                "canonical": False,
                "hgvsc": None,
                "hgvsp": None,
                "protein_position": None,
                "amino_acids": None,
                "codons": None,
                "exon": None,
                "intron": "3/5",
                "mane_select": None,
                "mane_plus_clinical": None,
                "lof": None,
                "lof_filter": None,
                "lof_flags": None,
                "lof_info": None,
                "sift_prediction": None,
                "sift_score": None,
                "polyphen_prediction": None,
                "polyphen_score": None,
            }
        ]
        output = format_annotations(annotations, "chr1:100:A:T")
        assert "HGVS.c" not in output
        assert "HGVS.p" not in output
        assert "LoF" not in output
        assert "SIFT" not in output
        assert "Intron" in output

    def test_empty_annotations(self):
        output = format_annotations([], "chr1:100:A:T")
        assert "No transcript annotations found" in output

    def test_lof_fields_displayed(self, lof_vep_response):
        annotations = extract_annotations(lof_vep_response, canonical_only=True)
        output = format_annotations(annotations, "chr17:41245466:G:A")
        assert "LoF" in output
        assert "HC" in output
        assert "SINGLE_EXON" in output


# --- TestQueryVep ---


class TestQueryVep:
    """Test VEP API queries (mocked)."""

    def _mock_urlopen(self, response_data, status=200):
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(response_data).encode("utf-8")
        mock_response.__enter__ = MagicMock(return_value=mock_response)
        mock_response.__exit__ = MagicMock(return_value=False)
        return mock_response

    @patch("pywombat.vep.urllib.request.urlopen")
    def test_successful_query(self, mock_urlopen, sample_vep_response):
        mock_urlopen.return_value = self._mock_urlopen([sample_vep_response])
        variants = [parse_variant("chr1:151776102:GC:AA")]
        result = query_vep(variants)
        assert len(result) == 1
        mock_urlopen.assert_called_once()

    @patch("pywombat.vep.urllib.request.urlopen")
    def test_request_body_contains_variant(self, mock_urlopen, sample_vep_response):
        mock_urlopen.return_value = self._mock_urlopen([sample_vep_response])
        variants = [parse_variant("chr1:151776102:GC:AA")]
        query_vep(variants)

        call_args = mock_urlopen.call_args
        request = call_args[0][0]
        body = json.loads(request.data.decode("utf-8"))
        assert "1 151776102 . GC AA . . ." in body["variants"]

    @patch("pywombat.vep.urllib.request.urlopen")
    def test_default_params_sent(self, mock_urlopen, sample_vep_response):
        mock_urlopen.return_value = self._mock_urlopen([sample_vep_response])
        variants = [parse_variant("chr1:100:A:T")]
        query_vep(variants)

        call_args = mock_urlopen.call_args
        request = call_args[0][0]
        body = json.loads(request.data.decode("utf-8"))
        assert body["canonical"] == 1
        assert body["LoF"] == 1
        assert body["hgvs"] == 1
        assert body["mane"] == 1

    @patch("pywombat.vep.urllib.request.urlopen")
    def test_http_error_raises_runtime_error(self, mock_urlopen):
        import urllib.error

        mock_urlopen.side_effect = urllib.error.HTTPError(
            url="https://rest.ensembl.org/vep/homo_sapiens/region",
            code=400,
            msg="Bad Request",
            hdrs={},
            fp=MagicMock(read=MagicMock(return_value=b"Invalid variant")),
        )
        variants = [parse_variant("chr1:100:A:T")]
        with pytest.raises(RuntimeError, match="VEP API error 400"):
            query_vep(variants)

    @patch("pywombat.vep.urllib.request.urlopen")
    def test_url_error_raises_connection_error(self, mock_urlopen):
        import urllib.error

        mock_urlopen.side_effect = urllib.error.URLError("Network unreachable")
        variants = [parse_variant("chr1:100:A:T")]
        with pytest.raises(ConnectionError, match="Cannot reach VEP API"):
            query_vep(variants)


# --- TestVepCLICommand ---


class TestVepCLICommand:
    """Test the wombat vep CLI command."""

    def test_vep_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "--help"])
        assert result.exit_code == 0
        assert "Annotate a variant" in result.output
        assert "--all" in result.output

    def test_vep_invalid_format(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "invalid"])
        assert result.exit_code != 0

    @patch("pywombat.vep.query_vep")
    def test_vep_canonical_only(self, mock_query, sample_vep_response):
        mock_query.return_value = [sample_vep_response]
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "chr1:151776102:GC:AA"])
        assert result.exit_code == 0
        assert "ENST00000000001" in result.output
        assert "(canonical)" in result.output
        # Non-canonical transcripts should not appear
        assert "ENST00000000002" not in result.output

    @patch("pywombat.vep.query_vep")
    def test_vep_all_transcripts(self, mock_query, sample_vep_response):
        mock_query.return_value = [sample_vep_response]
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "chr1:151776102:GC:AA", "--all"])
        assert result.exit_code == 0
        assert "ENST00000000001" in result.output
        assert "ENST00000000002" in result.output
        assert "ENST00000000003" in result.output

    @patch("pywombat.vep.query_vep")
    def test_vep_verbose(self, mock_query, sample_vep_response):
        mock_query.return_value = [sample_vep_response]
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "chr1:151776102:GC:AA", "-v"])
        assert result.exit_code == 0
        assert "Querying Ensembl VEP API" in result.output

    @patch("pywombat.vep.query_vep")
    def test_vep_connection_error(self, mock_query):
        mock_query.side_effect = ConnectionError("Cannot reach VEP API")
        runner = CliRunner()
        result = runner.invoke(cli, ["vep", "chr1:100:A:T"])
        assert result.exit_code != 0


# --- TestVepColumnMapping ---


class TestVepColumnMapping:
    """Test the VEP column mapping constant."""

    def test_mapping_has_expected_keys(self):
        expected_keys = {
            "gene_symbol", "gene_id", "transcript_id", "consequence",
            "impact", "biotype", "hgvsc", "hgvsp", "lof", "canonical",
        }
        assert expected_keys.issubset(VEP_COLUMN_MAPPING.keys())

    def test_mapping_values_are_vep_prefixed(self):
        for col_name in VEP_COLUMN_MAPPING.values():
            assert col_name.startswith("VEP_"), f"{col_name} doesn't start with VEP_"


# --- TestFindTranscriptAnnotation ---


class TestFindTranscriptAnnotation:
    """Test transcript lookup in VEP responses."""

    def test_exact_match(self, sample_vep_response):
        ann = _find_transcript_annotation(sample_vep_response, "ENST00000000001")
        assert ann is not None
        assert ann["transcript_id"] == "ENST00000000001"

    def test_version_stripped_match(self, sample_vep_response):
        """VEP_Feature may have version (ENST...001.3), API returns without."""
        ann = _find_transcript_annotation(sample_vep_response, "ENST00000000001.3")
        assert ann is not None
        assert ann["transcript_id"] == "ENST00000000001"

    def test_not_found_returns_none(self, sample_vep_response):
        ann = _find_transcript_annotation(sample_vep_response, "ENST99999999999")
        assert ann is None

    def test_intergenic_returns_none(self, intergenic_vep_response):
        ann = _find_transcript_annotation(intergenic_vep_response, "ENST00000000001")
        assert ann is None


# --- Fixtures for MNV annotation tests ---


@pytest.fixture
def mnv_dataframe():
    """Create a DataFrame simulating MNV detection output with VEP columns."""
    return pl.DataFrame(
        {
            "#CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "FILTER": [".", "."],
            "VEP_SYMBOL": ["GENE1", "GENE2"],
            "VEP_IMPACT": ["HIGH", "MODERATE"],
            "VEP_Consequence": ["stop_gained", "missense_variant"],
            "VEP_Feature": ["ENST00000000001", "ENST00000000002"],
            "VEP_CANONICAL": ["YES", "YES"],
            "VEP_HGVSc": ["ENST00000000001:c.100A>T", None],
            "VEP_LoF": ["HC", None],
            "sample": ["Sample1", "Sample1"],
            "sample_gt": ["0/1", "0/1"],
            "sample_dp": [30, 25],
            "sample_gq": [99, 85],
            "sample_ad": [15, 12],
            "sample_vaf": [0.5, 0.48],
            "father_gt": ["0/0", "0/0"],
            "father_dp": [25, 20],
            "mother_gt": ["0/0", "0/1"],
            "mother_dp": [28, 22],
            "mnv_proba": [0.95, None],
            "mnv_inframe": [None, None],
            "mnv_variant": ["chr1:100:AG:TC", ""],
        }
    )


@pytest.fixture
def mnv_vep_response():
    """VEP response for the MNV variant chr1:100:AG:TC."""
    return {
        "input": "1 100 . AG TC . . .",
        "most_severe_consequence": "synonymous_variant",
        "transcript_consequences": [
            {
                "gene_symbol": "GENE1",
                "gene_id": "ENSG00000000001",
                "transcript_id": "ENST00000000001",
                "consequence_terms": ["synonymous_variant"],
                "impact": "LOW",
                "canonical": 1,
                "biotype": "protein_coding",
                "hgvsc": "ENST00000000001.3:c.100_101delinsTC",
                "hgvsp": "ENSP00000000001:p.Ala34=",
                "amino_acids": "A",
                "codons": "gcT/gcC",
                "exon": "3/12",
                "lof": None,
            },
            {
                "gene_symbol": "GENE1",
                "gene_id": "ENSG00000000001",
                "transcript_id": "ENST00000000099",
                "consequence_terms": ["downstream_gene_variant"],
                "impact": "MODIFIER",
                "biotype": "retained_intron",
            },
        ],
    }


# --- TestAnnotateMnvVariants ---


class TestAnnotateMnvVariants:
    """Test MNV VEP annotation integration."""

    @patch("pywombat.vep.query_vep")
    def test_basic_annotation_creates_new_row(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)

        # Original 2 rows + 1 new MNV row
        assert result.shape[0] == 3

    @patch("pywombat.vep.query_vep")
    def test_new_row_has_mnv_coordinates(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)
        new_row = result.row(2, named=True)

        assert new_row["#CHROM"] == "chr1"
        assert new_row["POS"] == 100
        assert new_row["REF"] == "AG"
        assert new_row["ALT"] == "TC"

    @patch("pywombat.vep.query_vep")
    def test_sample_columns_copied(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)
        new_row = result.row(2, named=True)

        assert new_row["sample"] == "Sample1"
        assert new_row["sample_gt"] == "0/1"
        assert new_row["sample_dp"] == 30
        assert new_row["sample_vaf"] == 0.5

    @patch("pywombat.vep.query_vep")
    def test_parent_columns_null(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)
        new_row = result.row(2, named=True)

        assert new_row["father_gt"] is None
        assert new_row["father_dp"] is None
        assert new_row["mother_gt"] is None
        assert new_row["mother_dp"] is None

    @patch("pywombat.vep.query_vep")
    def test_mnv_proba_inverted(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)
        new_row = result.row(2, named=True)

        assert new_row["mnv_proba"] == pytest.approx(0.05)  # 1.0 - 0.95
        assert new_row["mnv_variant"] is None
        assert new_row["mnv_inframe"] is None

    @patch("pywombat.vep.query_vep")
    def test_vep_columns_populated(
        self, mock_query, mnv_dataframe, mnv_vep_response
    ):
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(mnv_dataframe)
        new_row = result.row(2, named=True)

        assert new_row["VEP_SYMBOL"] == "GENE1"
        assert new_row["VEP_IMPACT"] == "LOW"
        assert new_row["VEP_Consequence"] == "synonymous_variant"
        assert new_row["VEP_Feature"] == "ENST00000000001"
        assert new_row["VEP_CANONICAL"] == "YES"

    @patch("pywombat.vep.query_vep")
    def test_only_existing_columns_populated(self, mock_query, mnv_vep_response):
        """DataFrame without VEP_LoF should not get LoF in new rows."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "VEP_SYMBOL": ["GENE1"],
                "VEP_IMPACT": ["HIGH"],
                "VEP_Feature": ["ENST00000000001"],
                "sample": ["S1"],
                "sample_gt": ["0/1"],
                "sample_dp": [30],
                "sample_gq": [99],
                "sample_ad": [15],
                "sample_vaf": [0.5],
                "mnv_proba": [0.9],
                "mnv_inframe": [None],
                "mnv_variant": ["chr1:100:AG:TC"],
            }
        )
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(df)

        # VEP_LoF column doesn't exist in df, so it shouldn't appear
        assert "VEP_LoF" not in result.columns

    @patch("pywombat.vep.query_vep")
    def test_empty_mnv_variant_skipped(self, mock_query, mnv_dataframe):
        """Rows with empty mnv_variant should not trigger VEP queries."""
        # Set all mnv_variant to empty
        df = mnv_dataframe.with_columns(pl.lit("").alias("mnv_variant"))
        result = annotate_mnv_variants(df)

        assert result.shape[0] == 2  # No new rows
        mock_query.assert_not_called()

    @patch("pywombat.vep.query_vep")
    def test_null_vep_feature_skipped(self, mock_query, mnv_vep_response):
        """Rows with null VEP_Feature should be skipped."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "VEP_Feature": [None],
                "VEP_IMPACT": ["HIGH"],
                "sample": ["S1"],
                "sample_gt": ["0/1"],
                "sample_dp": [30],
                "sample_gq": [99],
                "sample_ad": [15],
                "sample_vaf": [0.5],
                "mnv_proba": [0.9],
                "mnv_inframe": [None],
                "mnv_variant": ["chr1:100:AG:TC"],
            }
        )
        result = annotate_mnv_variants(df)
        assert result.shape[0] == 1  # No new rows
        mock_query.assert_not_called()

    @patch("pywombat.vep.query_vep")
    def test_transcript_not_found_skips_row(
        self, mock_query, mnv_dataframe
    ):
        """If VEP doesn't return the target transcript, skip."""
        wrong_transcript_response = {
            "input": "1 100 . AG TC . . .",
            "transcript_consequences": [
                {
                    "gene_symbol": "OTHER",
                    "transcript_id": "ENST99999999999",
                    "consequence_terms": ["intron_variant"],
                    "impact": "MODIFIER",
                },
            ],
        }
        mock_query.return_value = [wrong_transcript_response]
        result = annotate_mnv_variants(mnv_dataframe)

        assert result.shape[0] == 2  # No new rows added

    @patch("pywombat.vep.query_vep")
    def test_api_failure_returns_original(self, mock_query, mnv_dataframe):
        """VEP API failure should return original DataFrame unchanged."""
        mock_query.side_effect = ConnectionError("Network error")
        result = annotate_mnv_variants(mnv_dataframe)

        assert result.shape[0] == 2  # Original rows only

    @patch("pywombat.vep.query_vep")
    def test_deduplication_queries_once(self, mock_query, mnv_vep_response):
        """Same mnv_variant for two samples should query VEP only once."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1", "chr1"],
                "POS": [100, 100],
                "REF": ["A", "A"],
                "ALT": ["T", "T"],
                "VEP_Feature": ["ENST00000000001", "ENST00000000001"],
                "VEP_IMPACT": ["HIGH", "HIGH"],
                "VEP_SYMBOL": ["GENE1", "GENE1"],
                "VEP_Consequence": ["stop_gained", "stop_gained"],
                "VEP_CANONICAL": ["YES", "YES"],
                "sample": ["Sample1", "Sample2"],
                "sample_gt": ["0/1", "0/1"],
                "sample_dp": [30, 28],
                "sample_gq": [99, 95],
                "sample_ad": [15, 14],
                "sample_vaf": [0.5, 0.5],
                "mnv_proba": [0.95, 0.90],
                "mnv_inframe": [None, None],
                "mnv_variant": ["chr1:100:AG:TC", "chr1:100:AG:TC"],
            }
        )
        mock_query.return_value = [mnv_vep_response]
        result = annotate_mnv_variants(df)

        # Should add 2 new rows (one per sample), but only 1 API call
        assert result.shape[0] == 4  # 2 original + 2 new
        mock_query.assert_called_once()

    @patch("pywombat.vep.query_vep")
    def test_batching(self, mock_query, mnv_vep_response):
        """More than batch_size variants should result in multiple API calls."""
        # Create 3 unique MNV variants
        rows = {
            "#CHROM": [f"chr1"] * 3,
            "POS": [100, 200, 300],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "A"],
            "VEP_Feature": ["ENST00000000001"] * 3,
            "VEP_IMPACT": ["HIGH"] * 3,
            "sample": ["S1"] * 3,
            "sample_gt": ["0/1"] * 3,
            "sample_dp": [30] * 3,
            "sample_gq": [99] * 3,
            "sample_ad": [15] * 3,
            "sample_vaf": [0.5] * 3,
            "mnv_proba": [0.95, 0.90, 0.85],
            "mnv_inframe": [None] * 3,
            "mnv_variant": ["chr1:100:AG:TC", "chr1:200:GC:CA", "chr1:300:CA:AG"],
        }
        df = pl.DataFrame(rows)

        # Make VEP responses for all 3
        resp1 = {**mnv_vep_response, "input": "1 100 . AG TC . . ."}
        resp2 = {**mnv_vep_response, "input": "1 200 . GC CA . . ."}
        resp3 = {**mnv_vep_response, "input": "1 300 . CA AG . . ."}

        mock_query.side_effect = [[resp1, resp2], [resp3]]
        result = annotate_mnv_variants(df, batch_size=2)

        # 2 API calls (batch of 2 + batch of 1)
        assert mock_query.call_count == 2

    @patch("pywombat.vep.query_vep")
    def test_colon_colon_colon_filtered(self, mock_query, mnv_dataframe):
        """Rows with mnv_variant=':::' (same-position edge case) should be skipped."""
        df = mnv_dataframe.with_columns(pl.lit(":::").alias("mnv_variant"))
        result = annotate_mnv_variants(df)

        assert result.shape[0] == 2  # No new rows
        mock_query.assert_not_called()

    def test_missing_vep_feature_column(self):
        """DataFrame without VEP_Feature column returns unchanged."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "sample": ["S1"],
                "mnv_proba": [0.95],
                "mnv_variant": ["chr1:100:AG:TC"],
            }
        )
        result = annotate_mnv_variants(df)
        assert result.shape[0] == 1


# --- TestNormalizeAnnotateMode ---


class TestNormalizeAnnotateMode:
    """Test _normalize_annotate_mode helper."""

    def test_true_returns_vep_api(self):
        assert _normalize_annotate_mode(True) == "vep_api"

    def test_false_returns_false(self):
        assert _normalize_annotate_mode(False) == "false"

    def test_none_returns_false(self):
        assert _normalize_annotate_mode(None) == "false"

    def test_string_true_returns_vep_api(self):
        assert _normalize_annotate_mode("true") == "vep_api"

    def test_string_false_returns_false(self):
        assert _normalize_annotate_mode("false") == "false"

    def test_string_vep_api(self):
        assert _normalize_annotate_mode("vep_api") == "vep_api"

    def test_string_external(self):
        assert _normalize_annotate_mode("external") == "external"

    def test_invalid_raises(self):
        with pytest.raises(Exception):
            _normalize_annotate_mode("invalid_value")


# --- TestCreateMnvPlaceholdersAndVcf ---


class TestCreateMnvPlaceholdersAndVcf:
    """Test create_mnv_placeholders_and_vcf."""

    @pytest.fixture
    def mnv_df(self):
        """DataFrame with MNV candidates for testing."""
        return pl.DataFrame(
            {
                "#CHROM": ["chr1", "chr1", "chr2"],
                "POS": [100, 200, 300],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "A"],
                "FILTER": [".", ".", "."],
                "VEP_SYMBOL": ["GENE1", "GENE2", "GENE3"],
                "VEP_IMPACT": ["HIGH", "MODERATE", "LOW"],
                "VEP_Consequence": ["stop_gained", "missense_variant", "synonymous_variant"],
                "VEP_Feature": ["ENST00000000001", "ENST00000000002", "ENST00000000003"],
                "VEP_CANONICAL": ["YES", "YES", "YES"],
                "sample": ["S1", "S1", "S1"],
                "sample_gt": ["0/1", "0/1", "0/1"],
                "sample_dp": [30, 25, 20],
                "sample_gq": [99, 85, 70],
                "sample_ad": [15, 12, 10],
                "sample_vaf": [0.5, 0.48, 0.5],
                "father_gt": ["0/0", "0/0", "0/0"],
                "father_dp": [25, 20, 18],
                "mother_gt": ["0/0", "0/0", "0/0"],
                "mother_dp": [28, 22, 20],
                "mnv_proba": [0.95, 0.80, None],
                "mnv_inframe": [None, None, None],
                "mnv_variant": ["chr1:100:AG:TC", "chr1:200:GC:CA", ""],
            }
        )

    def test_placeholder_rows_created(self, mnv_df, tmp_path):
        """Should create placeholder rows for eligible MNV candidates."""
        result, vcf_path = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        # 3 original + 2 placeholders (row 3 has no mnv_proba)
        assert result.shape[0] == 5

    def test_vcf_created_with_variants(self, mnv_df, tmp_path):
        """Should create VCF file with unique MNV variants."""
        _, vcf_path = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        assert vcf_path.exists()
        content = vcf_path.read_text()
        lines = content.strip().split("\n")
        assert lines[0].startswith("#CHROM")
        # 2 unique MNV variants
        assert len(lines) == 3  # header + 2 variants

    def test_vcf_variants_have_correct_format(self, mnv_df, tmp_path):
        """VCF variant lines should have ID = chr:pos:ref:alt."""
        _, vcf_path = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        lines = vcf_path.read_text().strip().split("\n")
        # First data line
        fields = lines[1].split("\t")
        assert fields[0] == "chr1"  # CHROM
        assert fields[1] == "100"  # POS
        assert fields[2] == "chr1:100:AG:TC"  # ID
        assert fields[3] == "AG"  # REF
        assert fields[4] == "TC"  # ALT

    def test_placeholder_mnv_proba_inverted(self, mnv_df, tmp_path):
        """Placeholder row mnv_proba should be 1.0 - original."""
        result, _ = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        placeholder_rows = result.filter(
            pl.col("VEP_Consequence").is_null()
        )
        probas = placeholder_rows["mnv_proba"].to_list()
        assert pytest.approx(0.05) in probas  # 1.0 - 0.95
        assert pytest.approx(0.20) in probas  # 1.0 - 0.80

    def test_placeholder_vep_feature_copied(self, mnv_df, tmp_path):
        """Placeholder row should have VEP_Feature copied from original."""
        result, _ = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        placeholder_rows = result.filter(
            pl.col("VEP_Consequence").is_null()
        )
        features = placeholder_rows["VEP_Feature"].to_list()
        assert "ENST00000000001" in features
        assert "ENST00000000002" in features

    def test_placeholder_vep_consequence_null(self, mnv_df, tmp_path):
        """Placeholder rows should have null VEP_Consequence."""
        result, _ = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        placeholder_rows = result.filter(
            pl.col("VEP_Consequence").is_null()
        )
        assert placeholder_rows.shape[0] == 2

    def test_placeholder_sample_columns_copied(self, mnv_df, tmp_path):
        """Sample columns should be copied from original row."""
        result, _ = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        placeholder_rows = result.filter(
            pl.col("VEP_Consequence").is_null()
        )
        first_placeholder = placeholder_rows.row(0, named=True)
        assert first_placeholder["sample"] == "S1"
        assert first_placeholder["sample_gt"] == "0/1"

    def test_placeholder_parent_columns_null(self, mnv_df, tmp_path):
        """Parent columns should be null in placeholder rows."""
        result, _ = create_mnv_placeholders_and_vcf(
            mnv_df, output_prefix=str(tmp_path / "out")
        )
        placeholder_rows = result.filter(
            pl.col("VEP_Consequence").is_null()
        )
        first_placeholder = placeholder_rows.row(0, named=True)
        assert first_placeholder["father_gt"] is None
        assert first_placeholder["mother_gt"] is None

    def test_vcf_deduplicated(self, tmp_path):
        """Same MNV variant from different samples should appear once in VCF."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1", "chr1"],
                "POS": [100, 100],
                "REF": ["A", "A"],
                "ALT": ["T", "T"],
                "VEP_Feature": ["ENST00000000001", "ENST00000000001"],
                "VEP_Consequence": ["stop_gained", "stop_gained"],
                "sample": ["S1", "S2"],
                "sample_gt": ["0/1", "0/1"],
                "sample_dp": [30, 28],
                "sample_gq": [99, 95],
                "sample_ad": [15, 14],
                "sample_vaf": [0.5, 0.5],
                "mnv_proba": [0.95, 0.90],
                "mnv_inframe": [None, None],
                "mnv_variant": ["chr1:100:AG:TC", "chr1:100:AG:TC"],
            }
        )
        _, vcf_path = create_mnv_placeholders_and_vcf(
            df, output_prefix=str(tmp_path / "out")
        )
        lines = vcf_path.read_text().strip().split("\n")
        # Header + 1 unique variant (deduplicated)
        assert len(lines) == 2

    def test_empty_mnv_returns_unchanged(self, tmp_path):
        """DataFrame with no eligible MNV rows returns unchanged."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1"],
                "POS": [100],
                "REF": ["A"],
                "ALT": ["T"],
                "VEP_Feature": ["ENST00000000001"],
                "VEP_Consequence": ["stop_gained"],
                "sample": ["S1"],
                "sample_gt": ["0/1"],
                "sample_dp": [30],
                "sample_gq": [99],
                "sample_ad": [15],
                "sample_vaf": [0.5],
                "mnv_proba": [None],
                "mnv_inframe": [None],
                "mnv_variant": [""],
            }
        )
        result, vcf_path = create_mnv_placeholders_and_vcf(
            df, output_prefix=str(tmp_path / "out")
        )
        assert result.shape[0] == 1  # Unchanged
        assert vcf_path.read_text() == ""  # Empty VCF


# --- TestParseVepAnnotationFile ---


class TestParseVepAnnotationFile:
    """Test parse_vep_annotation_file."""

    def test_parses_example_file(self):
        """Should parse the example VEP annotated file."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        assert len(lookup) > 0

    def test_skips_metadata_comments(self):
        """Lines starting with ## should be skipped."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        # Should have entries for the data rows only
        # The example has 11 data rows (lines 78-88)
        assert len(lookup) == 11

    def test_header_hash_stripped(self):
        """First column should have # stripped from header."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        # Keys should be (Uploaded_variation, transcript_base)
        keys = list(lookup.keys())
        assert all(isinstance(k, tuple) and len(k) == 2 for k in keys)
        # First variant should be chr1:65568:A:C
        variant_ids = [k[0] for k in keys]
        assert "chr1:65568:A:C" in variant_ids

    def test_dash_converted_to_none(self):
        """VEP '-' values should be converted to None."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        # The downstream_gene_variant for ENST00000492842 should have null HGVSc
        key = ("chr1:65568:A:C", "ENST00000492842")
        assert key in lookup
        assert lookup[key].get("VEP_HGVSc") is None

    def test_versioned_feature_preserved(self):
        """Feature value should preserve version (e.g. ENST00000641515.2)."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        # ENST00000641515.2 (canonical transcript of OR4F5)
        key = ("chr1:65568:A:C", "ENST00000641515")
        assert key in lookup
        assert lookup[key]["VEP_Feature"] == "ENST00000641515.2"

    def test_consequence_populated(self):
        """Should correctly populate consequence from VEP file."""
        example_path = Path(__file__).parent / "example.annotated.tsv"
        if not example_path.exists():
            pytest.skip("example.annotated.tsv not found")

        lookup = parse_vep_annotation_file(example_path)
        key = ("chr1:65568:A:C", "ENST00000641515")
        assert lookup[key]["VEP_Consequence"] == "missense_variant"
        assert lookup[key]["VEP_IMPACT"] == "MODERATE"
        assert lookup[key]["VEP_SYMBOL"] == "OR4F5"


# --- TestMergeVepAnnotations ---


class TestMergeVepAnnotations:
    """Test merge_vep_annotations."""

    @pytest.fixture
    def placeholder_tsv(self, tmp_path):
        """Create a TSV file with placeholder MNV rows."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr1", "chr1", "chr1"],
                "POS": [50000, 65568, 230710048],
                "REF": ["A", "A", "A"],
                "ALT": ["T", "C", "G"],
                "VEP_Consequence": ["missense_variant", None, None],
                "VEP_IMPACT": ["MODERATE", None, None],
                "VEP_SYMBOL": ["GENE1", None, None],
                "VEP_Feature": [
                    "ENST00000999999",
                    "ENST00000641515",
                    "ENST00000366667",
                ],
                "VEP_Gene": ["ENSG00000999999", None, None],
                "VEP_BIOTYPE": ["protein_coding", None, None],
                "VEP_HGVSc": ["ENST00000999999:c.1A>T", None, None],
                "VEP_HGVSp": [None, None, None],
                "VEP_CANONICAL": ["YES", None, None],
                "sample": ["S1", "S1", "S1"],
                "mnv_proba": [0.95, 0.05, 0.20],
            }
        )
        tsv_path = tmp_path / "filtered.tsv"
        df.write_csv(tsv_path, separator="\t")
        return tsv_path

    def test_basic_merge_fills_vep_columns(self, placeholder_tsv):
        """Placeholder rows should get VEP columns filled."""
        annotation_path = Path(__file__).parent / "example.annotated.tsv"
        if not annotation_path.exists():
            pytest.skip("example.annotated.tsv not found")

        result = merge_vep_annotations(placeholder_tsv, annotation_path)

        # Row 1 (chr1:65568:A:C, ENST00000641515) should now have VEP data
        row1 = result.row(1, named=True)
        assert row1["VEP_Consequence"] == "missense_variant"
        assert row1["VEP_SYMBOL"] == "OR4F5"
        assert row1["VEP_IMPACT"] == "MODERATE"

    def test_non_placeholder_rows_untouched(self, placeholder_tsv):
        """Non-placeholder rows should not be modified."""
        annotation_path = Path(__file__).parent / "example.annotated.tsv"
        if not annotation_path.exists():
            pytest.skip("example.annotated.tsv not found")

        result = merge_vep_annotations(placeholder_tsv, annotation_path)

        # Row 0 (original variant, not a placeholder) should be unchanged
        row0 = result.row(0, named=True)
        assert row0["VEP_Consequence"] == "missense_variant"
        assert row0["VEP_SYMBOL"] == "GENE1"

    def test_unmatched_placeholder_stays_null(self, tmp_path):
        """Placeholder with no matching VEP annotation stays null."""
        df = pl.DataFrame(
            {
                "#CHROM": ["chr99"],
                "POS": [999999],
                "REF": ["A"],
                "ALT": ["T"],
                "VEP_Consequence": [None],
                "VEP_IMPACT": [None],
                "VEP_Feature": ["ENST00000000000"],
                "sample": ["S1"],
                "mnv_proba": [0.05],
            }
        )
        tsv_path = tmp_path / "filtered.tsv"
        df.write_csv(tsv_path, separator="\t")

        annotation_path = Path(__file__).parent / "example.annotated.tsv"
        if not annotation_path.exists():
            pytest.skip("example.annotated.tsv not found")

        result = merge_vep_annotations(tsv_path, annotation_path)
        row = result.row(0, named=True)
        assert row["VEP_Consequence"] is None
        assert row["VEP_IMPACT"] is None

    def test_agt_variant_merged(self, placeholder_tsv):
        """chr1:230710048:A:G for ENST00000366667 should merge AGT annotations."""
        annotation_path = Path(__file__).parent / "example.annotated.tsv"
        if not annotation_path.exists():
            pytest.skip("example.annotated.tsv not found")

        result = merge_vep_annotations(placeholder_tsv, annotation_path)

        # Row 2 (chr1:230710048:A:G, ENST00000366667)
        row2 = result.row(2, named=True)
        assert row2["VEP_Consequence"] == "missense_variant"
        assert row2["VEP_SYMBOL"] == "AGT"
        assert row2["VEP_IMPACT"] == "MODERATE"

    def test_merge_annotations_cli(self, placeholder_tsv, tmp_path):
        """Test the merge-annotations CLI command."""
        annotation_path = Path(__file__).parent / "example.annotated.tsv"
        if not annotation_path.exists():
            pytest.skip("example.annotated.tsv not found")

        output_prefix = str(tmp_path / "merged")
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "merge-annotations",
                "--tsv", str(placeholder_tsv),
                "--annotations", str(annotation_path),
                "-o", output_prefix,
                "-v",
            ],
        )
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert (tmp_path / "merged.tsv").exists()

        # Verify the merged file has filled VEP data
        merged_df = pl.read_csv(tmp_path / "merged.tsv", separator="\t")
        row1 = merged_df.row(1, named=True)
        assert row1["VEP_SYMBOL"] == "OR4F5"
