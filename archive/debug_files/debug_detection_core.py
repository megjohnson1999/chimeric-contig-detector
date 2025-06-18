#!/usr/bin/env python3
"""
Debug script to test just the core detection methods without needing alignment.
"""

import sys
import traceback
import numpy as np
from pathlib import Path

# Add src to path so we can import
sys.path.insert(0, str(Path(__file__).parent / "src"))

from chimeric_detective.detector import ChimeraDetector

def test_core_detection_methods():
    """Test the core detection methods that are causing range() errors."""
    print("Testing core detection methods...")
    
    detector = ChimeraDetector(
        min_contig_length=50,
        min_coverage=1.0,
        window_size=20,
        step_size=10,
        log_level="DEBUG"
    )
    
    # Test sequences of different lengths
    test_sequences = [
        ("very_short", "ATCGATCGATCGATCG"),  # 16 bp
        ("short", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),  # 50 bp
        ("medium", "ATCGATCGATCGATCG" * 10),  # 160 bp
        ("long", "ATCGATCGATCGATCG" * 50),   # 800 bp
    ]
    
    for name, sequence in test_sequences:
        print(f"\n--- Testing {name} sequence (length: {len(sequence)}) ---")
        
        try:
            # Test coverage breakpoint detection with mock coverage
            print("Testing coverage detection...")
            mock_coverage = np.random.uniform(5, 15, len(sequence))  # Random coverage
            # Add a coverage drop in the middle to simulate a breakpoint
            mid = len(sequence) // 2
            mock_coverage[mid:] *= 0.3  # Reduce coverage in second half
            
            coverage_breakpoints = detector._detect_coverage_breakpoints(mock_coverage, window_size=20)
            print(f"  ✓ Coverage detection: found {len(coverage_breakpoints)} breakpoints")
            
        except Exception as e:
            print(f"  ❌ Coverage detection failed: {e}")
            print(f"  Error in line: {traceback.extract_tb(e.__traceback__)[-1].line}")
            
        try:
            # Test GC content detection
            print("Testing GC detection...")
            gc_profile = []
            window_size = 20
            step_size = 10
            
            # Create GC profile manually to avoid windowing issues
            for i in range(0, max(1, len(sequence) - window_size + 1), step_size):
                window_seq = sequence[i:i+window_size]
                if len(window_seq) >= window_size:
                    gc_content = (window_seq.count('G') + window_seq.count('C')) / len(window_seq)
                    gc_profile.append(gc_content)
            
            if len(gc_profile) > 1:
                gc_breakpoints = detector._detect_gc_breakpoints(gc_profile)
                print(f"  ✓ GC detection: found {len(gc_breakpoints)} breakpoints")
            else:
                print(f"  ⚠️  GC profile too short: {len(gc_profile)} windows")
            
        except Exception as e:
            print(f"  ❌ GC detection failed: {e}")
            print(f"  Error in line: {traceback.extract_tb(e.__traceback__)[-1].line}")
            
        try:
            # Test adaptive GC detection 
            print("Testing adaptive GC detection...")
            step_size = max(1, len(sequence) // 20)  # Adaptive step
            adaptive_gc_breakpoints = detector._detect_gc_breakpoints_adaptive(gc_profile, sequence, step_size)
            print(f"  ✓ Adaptive GC detection: found {len(adaptive_gc_breakpoints)} breakpoints")
            
        except Exception as e:
            print(f"  ❌ Adaptive GC detection failed: {e}")
            print(f"  Error in line: {traceback.extract_tb(e.__traceback__)[-1].line}")
            
        try:
            # Test k-mer detection
            print("Testing k-mer detection...")
            from chimeric_detective.utils import calculate_kmer_frequencies
            
            kmer_profile = []
            for i in range(0, max(1, len(sequence) - window_size + 1), step_size):
                window_seq = sequence[i:i+window_size]
                if len(window_seq) >= window_size:
                    kmers = calculate_kmer_frequencies(window_seq, k=4)
                    kmer_profile.append(kmers)
            
            if len(kmer_profile) > 1:
                kmer_breakpoints = detector._detect_kmer_breakpoints(kmer_profile)
                print(f"  ✓ K-mer detection: found {len(kmer_breakpoints)} breakpoints")
            else:
                print(f"  ⚠️  K-mer profile too short: {len(kmer_profile)} windows")
                
        except Exception as e:
            print(f"  ❌ K-mer detection failed: {e}")
            print(f"  Error in line: {traceback.extract_tb(e.__traceback__)[-1].line}")
            
        try:
            # Test breakpoint refinement if we have any breakpoints
            print("Testing breakpoint refinement...")
            test_position = len(sequence) // 2  # Middle of sequence
            refined_pos = detector._refine_breakpoint(sequence, test_position, window_size // 4)
            print(f"  ✓ Breakpoint refinement: {test_position} -> {refined_pos}")
            
        except Exception as e:
            print(f"  ❌ Breakpoint refinement failed: {e}")
            print(f"  Error in line: {traceback.extract_tb(e.__traceback__)[-1].line}")

if __name__ == "__main__":
    test_core_detection_methods()