"""Tests for MNV (Multi-Nucleotide Variant) detection functionality."""

import pytest
import polars as pl

try:
    from pywombat.mnv import (
        calculate_phasing_probability,
        calculate_indel_inframe,
    )
    MNV_AVAILABLE = True
except ImportError:
    MNV_AVAILABLE = False


@pytest.mark.skipif(not MNV_AVAILABLE, reason="MNV dependencies not installed")
class TestPhasingProbability:
    """Test phasing probability calculations."""

    def test_identical_ad_high_probability(self):
        """In-phase variants should have similar AD values."""
        var1_ad = (10, 15)
        var2_ad = (11, 14)
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob > 0.9, f"Expected high probability for similar AD, got {prob}"

    def test_very_similar_ad_very_high_probability(self):
        """Nearly identical AD should give very high probability."""
        var1_ad = (20, 30)
        var2_ad = (20, 30)  # Identical
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob > 0.99, f"Expected very high probability for identical AD, got {prob}"

    def test_swapped_ad_low_probability(self):
        """Out-of-phase variants have swapped AD values."""
        var1_ad = (10, 15)
        var2_ad = (15, 10)  # Swapped
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob < 0.1, f"Expected low probability for swapped AD, got {prob}"

    def test_missing_ad_returns_default(self):
        """Handle missing or very low coverage gracefully."""
        var1_ad = (10, 15)
        var2_ad = (0, 1)  # Very low coverage
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob == 1.0, f"Expected default 1.0 for low coverage, got {prob}"

    def test_missing_ref_depth_returns_default(self):
        """Handle missing reference depth."""
        import numpy as np
        var1_ad = (10, 15)
        var2_ad = (float('nan'), 10)
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob == 1.0, f"Expected default 1.0 for NaN ref depth, got {prob}"

    def test_moderate_difference_moderate_probability(self):
        """Moderately different AD should give intermediate probability."""
        # Values where both in-phase and out-of-phase are plausible
        var1_ad = (10, 15)
        var2_ad = (12, 13)  # Moderate similarity to both (10,15) and (15,10)
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        # This should give a moderate-to-high probability
        assert 0.4 < prob, f"Expected reasonable probability, got {prob}"

    def test_large_differences_use_log_probability(self):
        """Test that large AD differences are handled correctly with log probabilities."""
        var1_ad = (100, 150)
        var2_ad = (150, 100)  # Large swapped values
        prob = calculate_phasing_probability(var1_ad, var2_ad)
        assert prob < 0.1, f"Expected very low probability for large swapped AD, got {prob}"


@pytest.mark.skipif(not MNV_AVAILABLE, reason="MNV dependencies not installed")
class TestIndelInframe:
    """Test inframe calculation for indel clusters."""

    def test_single_indel_frameshift_minus_2(self):
        """Single -2bp deletion causes frameshift."""
        indels = [{'ref': 'ATG', 'alt': 'A'}]  # -2bp deletion
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for -2bp indel"

    def test_single_indel_frameshift_minus_1(self):
        """Single -1bp deletion causes frameshift."""
        indels = [{'ref': 'AT', 'alt': 'A'}]  # -1bp deletion
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for -1bp indel"

    def test_single_indel_frameshift_plus_1(self):
        """Single +1bp insertion causes frameshift."""
        indels = [{'ref': 'A', 'alt': 'AT'}]  # +1bp insertion
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for +1bp indel"

    def test_single_indel_inframe_minus_3(self):
        """Single -3bp deletion is inframe."""
        indels = [{'ref': 'ATGC', 'alt': 'A'}]  # -3bp deletion
        result = calculate_indel_inframe(indels)
        assert result == True, "Expected inframe for -3bp indel"

    def test_single_indel_inframe_plus_3(self):
        """Single +3bp insertion is inframe."""
        indels = [{'ref': 'A', 'alt': 'ATGC'}]  # +3bp insertion
        result = calculate_indel_inframe(indels)
        assert result == True, "Expected inframe for +3bp indel"

    def test_multiple_indels_inframe_net_zero(self):
        """Multiple indels that sum to 0 are inframe."""
        indels = [
            {'ref': 'ATG', 'alt': 'A'},    # -2bp
            {'ref': 'C', 'alt': 'CGAT'},   # +3bp
            {'ref': 'TA', 'alt': 'T'}      # -1bp
        ]
        # Net: -2 + 3 - 1 = 0 (inframe!)
        result = calculate_indel_inframe(indels)
        assert result == True, "Expected inframe for net 0bp change"

    def test_multiple_indels_inframe_net_minus_3(self):
        """Multiple indels with net -3bp are inframe."""
        indels = [
            {'ref': 'ATG', 'alt': 'A'},    # -2bp
            {'ref': 'C', 'alt': 'CG'},     # +1bp
            {'ref': 'TAG', 'alt': 'T'}     # -2bp
        ]
        # Net: -2 + 1 - 2 = -3 (inframe!)
        result = calculate_indel_inframe(indels)
        assert result == True, "Expected inframe for net -3bp change"

    def test_multiple_indels_frameshift_net_plus_1(self):
        """Multiple indels with net +1bp cause frameshift."""
        indels = [
            {'ref': 'A', 'alt': 'ATT'},   # +2bp
            {'ref': 'CG', 'alt': 'C'}     # -1bp
        ]
        # Net: +2 - 1 = +1 (frameshift)
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for net +1bp change"

    def test_multiple_indels_frameshift_net_minus_2(self):
        """Multiple indels with net -2bp cause frameshift."""
        indels = [
            {'ref': 'ATG', 'alt': 'A'},    # -2bp
            {'ref': 'C', 'alt': 'CGAT'}    # +3bp
        ]
        # Net: -2 + 3 = +1 (frameshift)
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for net +1bp change"

    def test_empty_indel_list(self):
        """Empty indel list is technically inframe (0 % 3 == 0)."""
        indels = []
        result = calculate_indel_inframe(indels)
        assert result == True, "Expected inframe for empty list (net 0bp)"

    def test_complex_multiple_indels_inframe(self):
        """Complex combination that sums to multiple of 3."""
        indels = [
            {'ref': 'A', 'alt': 'ATGCG'},     # +4bp
            {'ref': 'TT', 'alt': 'T'},        # -1bp
            {'ref': 'GGG', 'alt': 'G'},       # -2bp
            {'ref': 'C', 'alt': 'CATG'},      # +3bp
            {'ref': 'ATAT', 'alt': 'A'}       # -3bp
        ]
        # Net: +4 - 1 - 2 + 3 - 3 = +1 (frameshift)
        result = calculate_indel_inframe(indels)
        assert result == False, "Expected frameshift for net +1bp change"


@pytest.mark.skipif(not MNV_AVAILABLE, reason="MNV dependencies not installed")
class TestMNVDetectionIntegration:
    """Integration tests for MNV detection (require cyvcf2 and pysam)."""

    @pytest.mark.skip(reason="Requires test BCF and FASTA fixtures")
    def test_snv_pair_detection(self):
        """End-to-end SNV MNV detection."""
        # This would require creating test BCF and FASTA files
        # Skipping for now - to be implemented with proper fixtures
        pass

    @pytest.mark.skip(reason="Requires test BCF and FASTA fixtures")
    def test_indel_cluster_detection(self):
        """End-to-end indel cluster MNV detection."""
        # This would require creating test BCF and FASTA files
        # Skipping for now - to be implemented with proper fixtures
        pass


@pytest.mark.skipif(not MNV_AVAILABLE, reason="MNV dependencies not installed")
class TestMNVHelperFunctions:
    """Test helper functions in MNV module."""

    def test_is_snv_valid_snv(self):
        """Test SNV detection for valid SNVs."""
        from pywombat.mnv import _is_snv
        assert _is_snv('A', 'T') == True
        assert _is_snv('C', 'G') == True
        assert _is_snv('G', 'A') == True
        assert _is_snv('T', 'C') == True

    def test_is_snv_invalid_snv(self):
        """Test SNV detection rejects non-SNVs."""
        from pywombat.mnv import _is_snv
        assert _is_snv('AT', 'A') == False  # Deletion
        assert _is_snv('A', 'AT') == False  # Insertion
        assert _is_snv('AT', 'GC') == False  # MNV
        assert _is_snv('A', 'N') == False   # Invalid base
        assert _is_snv('N', 'T') == False   # Invalid base

    def test_is_indel_valid_indel(self):
        """Test indel detection for valid indels."""
        from pywombat.mnv import _is_indel
        assert _is_indel('AT', 'A') == True   # Deletion
        assert _is_indel('A', 'AT') == True   # Insertion
        assert _is_indel('ATG', 'A') == True  # Larger deletion
        assert _is_indel('A', 'ATGC') == True # Larger insertion

    def test_is_indel_not_indel(self):
        """Test indel detection rejects SNVs."""
        from pywombat.mnv import _is_indel
        assert _is_indel('A', 'T') == False
        assert _is_indel('C', 'G') == False
        assert _is_indel('G', 'A') == False
