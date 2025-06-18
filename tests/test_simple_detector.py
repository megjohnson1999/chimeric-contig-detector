"""
Tests for the SimpleChimeraDetector module.
"""

import unittest
import tempfile
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from chimeric_detective.detector_simple import SimpleChimeraDetector


class TestSimpleChimeraDetector(unittest.TestCase):
    """Test the SimpleChimeraDetector class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.detector = SimpleChimeraDetector(
            gc_content_threshold=0.1,
            window_size=500,
            step_size=250
        )
    
    def test_detector_initialization(self):
        """Test detector initialization."""
        self.assertEqual(self.detector.gc_content_threshold, 0.1)
        self.assertEqual(self.detector.window_size, 500)
        self.assertEqual(self.detector.step_size, 250)
    
    def test_empty_assembly(self):
        """Test detector with empty assembly file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            # Write empty file
            pass
        
        try:
            candidates = self.detector.detect_chimeras(f.name)
            self.assertEqual(len(candidates), 0)
        finally:
            os.unlink(f.name)
    
    def test_short_contigs(self):
        """Test detector with contigs too short to analyze."""
        # Create assembly with short contigs
        records = [
            SeqRecord(Seq("ATCG" * 50), id="short1"),  # 200bp
            SeqRecord(Seq("GCTA" * 50), id="short2"),  # 200bp
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            SeqIO.write(records, f, "fasta")
            f.flush()
            
            try:
                candidates = self.detector.detect_chimeras(f.name)
                self.assertEqual(len(candidates), 0)
            finally:
                os.unlink(f.name)
    
    def test_uniform_gc_contig(self):
        """Test detector with uniform GC content (no chimeras expected)."""
        # Create contig with uniform GC content (~50%)
        uniform_seq = "ATGC" * 500  # 2000bp, 50% GC
        records = [SeqRecord(Seq(uniform_seq), id="uniform")]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            SeqIO.write(records, f, "fasta")
            f.flush()
            
            try:
                candidates = self.detector.detect_chimeras(f.name)
                # Should detect no significant GC changes
                self.assertEqual(len(candidates), 0)
            finally:
                os.unlink(f.name)
    
    def test_gc_shift_chimera(self):
        """Test detector with clear GC content shift (chimera expected)."""
        # Create chimeric sequence: low GC + high GC
        low_gc = "ATAT" * 300   # 1200bp, 0% GC
        high_gc = "GCGC" * 300  # 1200bp, 100% GC
        chimeric_seq = low_gc + high_gc
        
        records = [SeqRecord(Seq(chimeric_seq), id="chimera")]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            SeqIO.write(records, f, "fasta")
            f.flush()
            
            try:
                candidates = self.detector.detect_chimeras(f.name)
                # Should detect the GC shift
                self.assertGreater(len(candidates), 0)
                
                # Check candidate properties
                candidate = candidates[0]
                self.assertEqual(candidate.contig_id, "chimera")
                self.assertGreaterEqual(candidate.gc_difference, 0.5)  # Large difference
                self.assertGreater(candidate.confidence_score, 0.5)   # High confidence
                
            finally:
                os.unlink(f.name)
    
    def test_adaptive_window_parameters(self):
        """Test adaptive window parameter calculation."""
        # Test with different sequence lengths
        short_seq = "ATGC" * 250  # 1000bp
        long_seq = "ATGC" * 2500  # 10000bp
        
        short_params = self.detector._calculate_adaptive_window_parameters(short_seq)
        long_params = self.detector._calculate_adaptive_window_parameters(long_seq)
        
        # Longer sequences should get larger windows
        self.assertGreaterEqual(long_params['window'], short_params['window'])
        self.assertGreaterEqual(long_params['step'], short_params['step'])
    
    def test_gc_profile_calculation(self):
        """Test GC profile calculation."""
        # Create sequence with known GC content
        test_seq = "AAAA" + "GGGG"  # 0% GC + 100% GC
        
        profile = self.detector._calculate_gc_profile(test_seq, window_size=4, step_size=1)
        
        # Should have profiles for overlapping windows
        self.assertGreater(len(profile), 0)
        # First window should be 0% GC, some later window should be higher
        self.assertEqual(profile[0], 0.0)
        self.assertGreater(max(profile), 0.5)


if __name__ == '__main__':
    unittest.main()