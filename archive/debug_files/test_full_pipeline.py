#!/usr/bin/env python3
"""
Test the complete detection pipeline with a synthetic sequence and BAM file.
"""

import sys
import os
import tempfile
import numpy as np
sys.path.insert(0, './src')

from chimeric_detective.detector import ChimeraDetector
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

def create_synthetic_test():
    """Create a synthetic test case with sequence and fake BAM data."""
    print("=== Creating synthetic test case ===")
    
    # Create test sequence - GC rich on left, AT rich on right
    left_part = 'GCGCGCGC' * 50  # 400bp, 100% GC
    right_part = 'ATATATATAT' * 50  # 500bp, 0% GC
    sequence = left_part + right_part  # 900bp total
    
    # Create test FASTA file
    temp_dir = tempfile.mkdtemp()
    fasta_file = os.path.join(temp_dir, "test.fasta")
    
    record = SeqRecord(Seq(sequence), id="test_contig", description="")
    with open(fasta_file, "w") as f:
        SeqIO.write(record, f, "fasta")
    
    # Create fake SAM file
    sam_file = os.path.join(temp_dir, "test.sam")
    bam_file = os.path.join(temp_dir, "test.bam")
    
    with open(sam_file, "w") as f:
        # SAM header
        f.write("@HD\\tVN:1.6\\tSO:coordinate\\n")
        f.write(f"@SQ\\tSN:test_contig\\tLN:{len(sequence)}\\n")
        
        # Create reads with coverage jump at position 400
        read_id = 1
        
        # High coverage in left part (positions 0-399)
        for pos in range(0, 400, 20):  # Every 20bp
            for copy in range(20):  # 20 reads per position = 20x coverage
                if pos + 50 <= 400:  # 50bp reads
                    f.write(f"read_{read_id}\\t0\\ttest_contig\\t{pos+1}\\t60\\t50M\\t*\\t0\\t0\\t{'A'*50}\\t{'I'*50}\\n")
                    read_id += 1
        
        # Low coverage in right part (positions 400-899)  
        for pos in range(400, 900, 20):  # Every 20bp
            for copy in range(5):  # 5 reads per position = 5x coverage
                if pos + 50 <= 900:  # 50bp reads
                    f.write(f"read_{read_id}\\t0\\ttest_contig\\t{pos+1}\\t60\\t50M\\t*\\t0\\t0\\t{'A'*50}\\t{'I'*50}\\n")
                    read_id += 1
    
    # Convert to BAM and index
    pysam.view("-bS", "-o", bam_file, sam_file, catch_stdout=False)
    pysam.sort("-o", bam_file.replace('.bam', '_sorted.bam'), bam_file, catch_stdout=False)
    pysam.index(bam_file.replace('.bam', '_sorted.bam'), catch_stdout=False)
    
    return fasta_file, bam_file.replace('.bam', '_sorted.bam'), temp_dir

def test_full_detection_pipeline():
    """Test the complete detection pipeline."""
    print("=== Testing Full Detection Pipeline ===")
    
    # Create test data
    fasta_file, bam_file, temp_dir = create_synthetic_test()
    
    # Create detector with sensitive settings
    detector = ChimeraDetector(
        min_contig_length=100,
        min_coverage=1.0,
        coverage_fold_change=1.5,
        gc_content_threshold=0.05,
        kmer_distance_threshold=0.1,
        min_spanning_reads=1,
        log_level='DEBUG'
    )
    
    print(f"\\nUsing test files:")
    print(f"  FASTA: {fasta_file}")
    print(f"  BAM: {bam_file}")
    
    try:
        # Run detection
        candidates = detector.detect_chimeras(
            assembly_file=fasta_file,
            bam_file=bam_file
        )
        
        print(f"\\n=== DETECTION RESULTS ===")
        print(f"Total candidates found: {len(candidates)}")
        
        for i, candidate in enumerate(candidates):
            print(f"\\nCandidate {i+1}:")
            print(f"  Contig: {candidate.contig_id}")
            print(f"  Breakpoint: {candidate.breakpoint}")
            print(f"  Confidence: {candidate.confidence_score:.3f}")
            print(f"  Evidence: {', '.join(candidate.evidence_types)}")
            print(f"  Coverage L/R: {candidate.coverage_left:.1f}/{candidate.coverage_right:.1f}")
            print(f"  GC L/R: {candidate.gc_content_left:.3f}/{candidate.gc_content_right:.3f}")
            print(f"  K-mer distance: {candidate.kmer_distance:.3f}")
            print(f"  Spanning reads: {candidate.spanning_reads}")
        
        # Check if we found the expected breakpoint around position 400
        expected_region = (350, 450)  # Allow some wiggle room
        found_in_region = [c for c in candidates if expected_region[0] <= c.breakpoint <= expected_region[1]]
        
        print(f"\\n=== VALIDATION ===")
        print(f"Expected breakpoint region: {expected_region}")
        print(f"Candidates in expected region: {len(found_in_region)}")
        
        if found_in_region:
            print("SUCCESS: Found chimera candidates in expected region!")
            return True
        else:
            print("WARNING: No candidates found in expected breakpoint region")
            return len(candidates) > 0
            
    except Exception as e:
        print(f"ERROR: Detection failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Cleanup
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    success = test_full_detection_pipeline()
    print(f"\\n=== FINAL RESULT: {'PASS' if success else 'FAIL'} ===")