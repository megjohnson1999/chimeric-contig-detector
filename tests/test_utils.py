"""
Tests for utility functions.
"""

import pytest
import tempfile
import os
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

from chimeric_detective.utils import (
    calculate_gc_content, calculate_kmer_frequencies, calculate_kmer_distance,
    run_command, merge_overlapping_intervals, setup_logging
)


class TestGCContent:
    """Test GC content calculation functions."""
    
    def test_calculate_gc_content_basic(self):
        """Test basic GC content calculation."""
        sequence = "ATCGATCGATCGATCG"  # 50% GC
        gc_content = calculate_gc_content(sequence, window_size=8)
        
        # Each window should have roughly 50% GC
        for gc in gc_content:
            assert 0.0 <= gc <= 1.0
            assert abs(gc - 0.5) <= 0.1  # Allow some variance due to windowing
    
    def test_calculate_gc_content_edge_cases(self):
        """Test GC content calculation with edge cases."""
        # Empty sequence
        gc_content = calculate_gc_content("", window_size=4)
        assert gc_content == []
        
        # Very short sequence
        gc_content = calculate_gc_content("AT", window_size=4)
        assert len(gc_content) == 1
        assert gc_content[0] == 0.0  # No GC in "AT"
        
        # All GC sequence
        gc_content = calculate_gc_content("GGGGCCCC", window_size=4)
        for gc in gc_content:
            assert gc == 1.0
        
        # No GC sequence
        gc_content = calculate_gc_content("AAAATTTT", window_size=4)
        for gc in gc_content:
            assert gc == 0.0
    
    def test_calculate_gc_content_different_window_sizes(self):
        """Test GC content with different window sizes."""
        sequence = "ATCGATCGATCGATCG" * 10  # 160bp sequence
        
        # Small window
        gc_small = calculate_gc_content(sequence, window_size=4)
        # Large window
        gc_large = calculate_gc_content(sequence, window_size=40)
        
        assert len(gc_small) > len(gc_large)
        
        # Both should have valid GC content
        for gc in gc_small + gc_large:
            assert 0.0 <= gc <= 1.0


class TestKmerFunctions:
    """Test k-mer related functions."""
    
    def test_calculate_kmer_frequencies_basic(self):
        """Test basic k-mer frequency calculation."""
        sequence = "AAATTTCCCGGG"
        kmers = calculate_kmer_frequencies(sequence, k=3)
        
        # Should have overlapping 3-mers
        expected_kmers = ["AAA", "AAT", "ATT", "TTT", "TTC", "TCC", "CCC", "CCG", "CGG", "GGG"]
        
        assert len(kmers) <= len(expected_kmers)  # Some might not be present due to filtering
        for kmer, count in kmers.items():
            assert len(kmer) == 3
            assert count > 0
            assert "N" not in kmer
    
    def test_calculate_kmer_frequencies_with_n(self):
        """Test k-mer frequency calculation with N bases."""
        sequence = "AAANNNTTTCCC"
        kmers = calculate_kmer_frequencies(sequence, k=3)
        
        # K-mers containing N should be filtered out
        for kmer in kmers.keys():
            assert "N" not in kmer
    
    def test_calculate_kmer_frequencies_edge_cases(self):
        """Test k-mer frequency calculation edge cases."""
        # Empty sequence
        kmers = calculate_kmer_frequencies("", k=3)
        assert kmers == {}
        
        # Sequence shorter than k
        kmers = calculate_kmer_frequencies("AT", k=3)
        assert kmers == {}
        
        # Sequence exactly k length
        kmers = calculate_kmer_frequencies("ATG", k=3)
        assert kmers == {"ATG": 1}
    
    def test_calculate_kmer_distance_basic(self):
        """Test k-mer distance calculation."""
        kmers1 = {"AAA": 10, "TTT": 10, "GGG": 5, "CCC": 5}
        kmers2 = {"AAA": 5, "TTT": 5, "GGG": 10, "CCC": 10}
        
        distance = calculate_kmer_distance(kmers1, kmers2)
        
        assert 0.0 <= distance <= 1.0
        assert distance > 0.0  # Should be different distributions
    
    def test_calculate_kmer_distance_edge_cases(self):
        """Test k-mer distance calculation edge cases."""
        # Identical distributions
        kmers1 = {"AAA": 10, "TTT": 10}
        kmers2 = {"AAA": 10, "TTT": 10}
        
        distance = calculate_kmer_distance(kmers1, kmers2)
        assert distance == 0.0
        
        # Empty distributions
        distance = calculate_kmer_distance({}, {})
        assert distance == 0.0
        
        # One empty distribution
        distance = calculate_kmer_distance({"AAA": 10}, {})
        assert distance >= 0.0


class TestUtilityFunctions:
    """Test other utility functions."""
    
    def test_merge_overlapping_intervals(self):
        """Test interval merging function."""
        # Non-overlapping intervals with large gaps (default min_gap=100)
        intervals = [(1, 5), (200, 205), (400, 405)]
        merged = merge_overlapping_intervals(intervals)
        assert merged == intervals
        
        # Overlapping intervals
        intervals = [(1, 5), (3, 8), (200, 205), (202, 208)]
        merged = merge_overlapping_intervals(intervals)
        expected = [(1, 8), (200, 208)]
        assert merged == expected
        
        # Adjacent intervals with gap threshold
        intervals = [(1, 5), (7, 10)]
        merged = merge_overlapping_intervals(intervals, min_gap=5)
        expected = [(1, 10)]  # Should merge due to small gap
        assert merged == expected
    
    def test_merge_overlapping_intervals_edge_cases(self):
        """Test interval merging edge cases."""
        # Empty list
        merged = merge_overlapping_intervals([])
        assert merged == []
        
        # Single interval
        merged = merge_overlapping_intervals([(1, 5)])
        assert merged == [(1, 5)]
        
        # Identical intervals
        merged = merge_overlapping_intervals([(1, 5), (1, 5)])
        assert merged == [(1, 5)]
    
    @patch('subprocess.run')
    def test_run_command_success(self, mock_run):
        """Test successful command execution."""
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="test output",
            stderr=""
        )
        
        result = run_command(["echo", "test"])
        
        assert result.returncode == 0
        assert result.stdout == "test output"
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_run_command_failure(self, mock_run):
        """Test failed command execution."""
        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="error message"
        )
        
        result = run_command(["false"])
        
        assert result.returncode == 1
        assert result.stderr == "error message"
    
    @patch('subprocess.run')
    def test_run_command_exception_handling(self, mock_run):
        """Test command exception handling."""
        mock_run.side_effect = FileNotFoundError("Command not found")
        
        with pytest.raises(FileNotFoundError):
            run_command(["nonexistent_command"])
    
    def test_setup_logging(self):
        """Test logging setup."""
        # Test different log levels
        logger = setup_logging("DEBUG")
        assert logger.level <= 10  # DEBUG level
        
        logger = setup_logging("INFO") 
        assert logger.level <= 20  # INFO level
        
        logger = setup_logging("ERROR")
        assert logger.level <= 40  # ERROR level


class TestGCContentComprehensive:
    """Comprehensive GC content tests."""
    
    @pytest.mark.parametrize("sequence,window_size,expected_first_gc", [
        ("ATCG", 4, 0.5),       # Window: ATCG = 2/4 = 0.5
        ("AAAA", 4, 0.0),       # Window: AAAA = 0/4 = 0.0
        ("GGGG", 4, 1.0),       # Window: GGGG = 4/4 = 1.0
        ("ATCGATCGATCG", 6, 0.33), # Window: ATCGAT = 2/6 = 0.33
        ("GGGGGGGGAAAA", 8, 1.0), # Window: GGGGGGGG = 8/8 = 1.0
    ])
    def test_gc_content_parametrized(self, sequence, window_size, expected_first_gc):
        """Test GC content calculation with parametrized inputs."""
        gc_content = calculate_gc_content(sequence, window_size=window_size)
        
        # Should have at least one window
        assert len(gc_content) >= 1
        # Check first window
        assert abs(gc_content[0] - expected_first_gc) < 0.01


class TestKmerDistanceComprehensive:
    """Comprehensive k-mer distance tests."""
    
    def test_kmer_distance_symmetric(self):
        """Test that k-mer distance is symmetric."""
        kmers1 = {"AAA": 10, "TTT": 5, "GGG": 3}
        kmers2 = {"AAA": 5, "TTT": 10, "CCC": 3}
        
        dist1 = calculate_kmer_distance(kmers1, kmers2)
        dist2 = calculate_kmer_distance(kmers2, kmers1)
        
        assert abs(dist1 - dist2) < 1e-10  # Should be identical
    
    def test_kmer_distance_bounds(self):
        """Test that k-mer distance is properly bounded."""
        # Completely different distributions
        kmers1 = {"AAA": 100}
        kmers2 = {"TTT": 100}
        
        distance = calculate_kmer_distance(kmers1, kmers2)
        assert 0.0 <= distance <= 1.0
        
        # Should be close to maximum distance
        assert distance > 0.5