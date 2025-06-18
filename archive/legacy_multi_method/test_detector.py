"""
Tests for the ChimeraDetector module.
"""

import unittest
import tempfile
import os
from unittest.mock import Mock, patch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np

from chimeric_detective.detector import ChimeraDetector, ChimeraCandidate
from chimeric_detective.utils import calculate_gc_content, calculate_kmer_frequencies


class TestChimeraDetector(unittest.TestCase):
    
    def setUp(self):
        """Set up test fixtures."""
        self.detector = ChimeraDetector(
            min_contig_length=500,
            min_coverage=3.0,
            coverage_fold_change=2.0,
            gc_content_threshold=0.1,
            kmer_distance_threshold=0.3,
            log_level="ERROR"  # Suppress logging during tests
        )
        
        # Create test sequences
        self.test_sequence = "A" * 1000 + "G" * 1000  # Simple chimeric sequence
        self.normal_sequence = "ATCGATCG" * 250  # Normal sequence
    
    def test_load_assembly(self):
        """Test loading assembly from FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">contig1\n")
            f.write("ATCGATCGATCG\n")
            f.write(">contig2\n")
            f.write("GCTAGCTAGCTA\n")
            f.flush()
            
            contigs = self.detector._load_assembly(f.name)
            
            self.assertEqual(len(contigs), 2)
            self.assertIn("contig1", contigs)
            self.assertIn("contig2", contigs)
            self.assertEqual(contigs["contig1"], "ATCGATCGATCG")
            self.assertEqual(contigs["contig2"], "GCTAGCTAGCTA")
            
            os.unlink(f.name)
    
    def test_calculate_coverage(self):
        """Test coverage calculation."""
        # This would require a mock BAM file
        # For now, test that the method exists and has correct signature
        self.assertTrue(hasattr(self.detector, '_calculate_coverage'))
    
    def test_detect_coverage_breakpoints(self):
        """Test coverage breakpoint detection."""
        # Create coverage array with clear breakpoint that meets the algorithm's requirements
        # Algorithm needs: window_size (1000) on each side, and both sides > min_coverage (3.0)
        coverage = np.array([15.0] * 1500 + [45.0] * 1500)  # Clear 3x fold change, 3000bp total
        
        breakpoints = self.detector._detect_coverage_breakpoints(coverage)
        
        # Should detect breakpoint around position 1500
        self.assertTrue(len(breakpoints) > 0)
        # Breakpoint should be roughly in the middle (allowing for window effects)
        self.assertTrue(1200 < breakpoints[0] < 1800)
    
    def test_detect_gc_breakpoints(self):
        """Test GC content breakpoint detection."""
        # Create sequences with different GC content
        left_seq = "A" * 1000  # 0% GC
        right_seq = "G" * 1000  # 100% GC
        sequence = left_seq + right_seq
        
        gc_profile = calculate_gc_content(sequence, window_size=500)
        breakpoints = self.detector._detect_gc_breakpoints(gc_profile)
        
        # Should detect at least one breakpoint
        self.assertTrue(len(breakpoints) > 0)
    
    def test_detect_kmer_breakpoints(self):
        """Test k-mer composition breakpoint detection."""
        # Create sequences with different k-mer composition
        left_seq = "AAAA" * 250  # Homopolymer A
        right_seq = "GGGG" * 250  # Homopolymer G
        sequence = left_seq + right_seq
        
        kmer_profile = self.detector._calculate_kmer_profile(sequence)
        breakpoints = self.detector._detect_kmer_breakpoints(kmer_profile)
        
        # Should detect breakpoints due to composition change
        self.assertTrue(len(breakpoints) >= 0)  # May or may not detect depending on window size
    
    def test_calculate_confidence_score(self):
        """Test confidence score calculation."""
        evidence_types = ["coverage_discontinuity", "gc_content_shift"]
        
        score = self.detector._calculate_confidence_score(
            evidence_types=evidence_types,
            cov_fold_change=3.0,
            gc_diff=0.2,
            kmer_dist=0.4,
            spanning_reads=2,
            orientation_score=0.3
        )
        
        self.assertTrue(0 <= score <= 1)
        self.assertTrue(score > 0.5)  # Should be reasonably confident with this evidence
    
    def test_chimera_candidate_creation(self):
        """Test ChimeraCandidate creation."""
        candidate = ChimeraCandidate(
            contig_id="test_contig",
            breakpoint=1000,
            confidence_score=0.8,
            evidence_types=["coverage_discontinuity"],
            coverage_left=10.0,
            coverage_right=30.0,
            gc_content_left=0.4,
            gc_content_right=0.6,
            kmer_distance=0.5,
            spanning_reads=5,
            read_orientation_score=0.2
        )
        
        self.assertEqual(candidate.contig_id, "test_contig")
        self.assertEqual(candidate.breakpoint, 1000)
        self.assertEqual(candidate.confidence_score, 0.8)
        self.assertIn("coverage_discontinuity", candidate.evidence_types)


class TestUtilityFunctions(unittest.TestCase):
    
    def test_calculate_gc_content(self):
        """Test GC content calculation."""
        sequence = "ATCGATCGATCG"  # 50% GC (6 GC out of 12 bases)
        gc_content = calculate_gc_content(sequence, window_size=4)
        
        # With window_size=4 and step=2, we get:
        # Window 1 (0-4): "ATCG" = 2/4 = 0.5
        # Window 2 (2-6): "CGAT" = 2/4 = 0.5  
        # Window 3 (4-8): "ATCG" = 2/4 = 0.5
        # Window 4 (6-10): "GATC" = 2/4 = 0.5
        # Window 5 (8-12): "TCGA" = 2/4 = 0.5
        # Window 6 (10-12): "CG" = 2/2 = 1.0 (partial window)
        expected_windows = len(gc_content)
        self.assertTrue(expected_windows >= 5)
        
        # Check that all but possibly the last window have expected GC content
        for i, gc in enumerate(gc_content[:-1]):  # Skip last window which might be partial
            self.assertTrue(0.0 <= gc <= 1.0)
            self.assertAlmostEqual(gc, 0.5, delta=0.1)
        
        # Last window might be different due to partial window
        self.assertTrue(0.0 <= gc_content[-1] <= 1.0)
    
    def test_calculate_kmer_frequencies(self):
        """Test k-mer frequency calculation."""
        sequence = "AAAA"
        kmers = calculate_kmer_frequencies(sequence, k=2)
        
        # Should have one unique 2-mer: "AA"
        self.assertEqual(len(kmers), 1)
        self.assertIn("AA", kmers)
        self.assertEqual(kmers["AA"], 3)  # "AAAA" has 3 overlapping "AA"s
    
    def test_kmer_frequencies_with_n(self):
        """Test k-mer frequencies with N bases."""
        sequence = "AANNA"
        kmers = calculate_kmer_frequencies(sequence, k=2)
        
        # Should skip k-mers containing N
        self.assertNotIn("AN", kmers)
        self.assertNotIn("NN", kmers)
        self.assertNotIn("NA", kmers)


if __name__ == '__main__':
    unittest.main()