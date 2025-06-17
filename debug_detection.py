#!/usr/bin/env python3
"""
Debug chimera detection to find why it's detecting 0 candidates.
"""

import sys
import os
sys.path.insert(0, './src')

from chimeric_detective.detector import ChimeraDetector
import numpy as np

def test_basic_detection():
    """Test the detection logic with a simple synthetic case."""
    print("=== DEBUG: Testing basic detection logic ===")
    
    # Create detector with debug logging
    detector = ChimeraDetector(
        min_contig_length=100,
        min_coverage=1.0,
        coverage_fold_change=1.5,
        gc_content_threshold=0.05,
        kmer_distance_threshold=0.1,
        log_level='DEBUG'
    )
    
    # Create a simple test sequence - GC rich on left, AT rich on right
    sequence = 'GCGCGCGC' * 50 + 'ATATATATAT' * 50  # 400bp + 500bp = 900bp total
    print(f"Test sequence length: {len(sequence)}bp")
    
    # Test the individual components
    print("\n=== Testing validation ===")
    valid = detector._is_valid_contig_for_analysis("test_contig", sequence)
    print(f"Contig valid for analysis: {valid}")
    
    print("\n=== Testing window parameters ===")
    window_params = detector._calculate_adaptive_window_parameters(sequence)
    print(f"Window parameters: {window_params}")
    
    print("\n=== Testing GC profile calculation ===")
    gc_profile = detector._calculate_adaptive_gc_profile(sequence, window_params['window'], window_params['step'])
    print(f"GC profile length: {len(gc_profile)}")
    if gc_profile:
        print(f"GC profile values: {gc_profile[:5]}...{gc_profile[-5:] if len(gc_profile) > 5 else gc_profile}")
    
    print("\n=== Testing k-mer profile calculation ===")
    kmer_profile = detector._calculate_adaptive_kmer_profile(sequence, window_params['window'], window_params['step'])
    print(f"K-mer profile length: {len(kmer_profile)}")
    
    print("\n=== Testing GC breakpoint detection ===")
    gc_breakpoints = detector._detect_gc_breakpoints_adaptive(gc_profile, sequence, window_params['step'])
    print(f"GC breakpoints found: {len(gc_breakpoints)} - {gc_breakpoints}")
    
    print("\n=== Testing k-mer breakpoint detection ===")
    kmer_breakpoints = detector._detect_kmer_breakpoints_adaptive(kmer_profile, sequence, window_params['step'])
    print(f"K-mer breakpoints found: {len(kmer_breakpoints)} - {kmer_breakpoints}")
    
    # Create fake coverage with a breakpoint
    print("\n=== Testing coverage breakpoint detection ===")
    coverage = np.ones(len(sequence)) * 10  # Base coverage of 10
    coverage[400:] *= 3  # 3x coverage jump at position 400
    coverage_breakpoints = detector._detect_coverage_breakpoints_adaptive(coverage, window_params['window'])
    print(f"Coverage breakpoints found: {len(coverage_breakpoints)} - {coverage_breakpoints}")
    
    # Test the combined breakpoint detection
    print("\n=== Testing combined breakpoint detection ===")
    profiles = {
        'coverage': coverage,
        'gc': gc_profile,
        'kmer': kmer_profile
    }
    all_breakpoints = detector._detect_all_breakpoints(sequence, profiles, window_params)
    print(f"Total breakpoints found: {len(all_breakpoints)} - {sorted(all_breakpoints)}")
    
    return len(all_breakpoints) > 0

if __name__ == "__main__":
    success = test_basic_detection()
    print(f"\n=== RESULT: {'PASS' if success else 'FAIL'} ===")