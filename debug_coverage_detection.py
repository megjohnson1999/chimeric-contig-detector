#!/usr/bin/env python3
"""
Debug coverage breakpoint detection specifically.
"""

import sys
import numpy as np
sys.path.insert(0, './src')

from chimeric_detective.detector import ChimeraDetector

def debug_coverage_detection():
    """Debug why coverage breakpoint detection isn't working."""
    print("=== DEBUGGING COVERAGE BREAKPOINT DETECTION ===")
    
    detector = ChimeraDetector(
        min_coverage=1.0,
        coverage_fold_change=1.5,
        log_level='DEBUG'
    )
    
    # Create clear coverage jump: 10x -> 30x at position 400
    sequence_length = 900
    coverage = np.ones(sequence_length) * 10  # Base coverage of 10
    coverage[400:] *= 3  # 3x coverage jump at position 400 (10 -> 30)
    
    print(f"Coverage profile: {coverage[:5]}...{coverage[395:405]}...{coverage[-5:]}")
    print(f"Coverage before breakpoint (395-399): {coverage[395:399]}")
    print(f"Coverage after breakpoint (400-404): {coverage[400:404]}")
    
    # Test coverage detection
    window_size = 200
    print(f"Using window size: {window_size}")
    print(f"Detector settings: min_coverage={detector.min_coverage}, fold_change={detector.coverage_fold_change}")
    print(f"detector.window_size={detector.window_size}")
    
    breakpoints = detector._detect_coverage_breakpoints(coverage, window_size)
    print(f"Detected breakpoints: {breakpoints}")
    
    # Debug the smoothing step
    window = np.ones(window_size) / window_size
    smoothed_coverage = np.convolve(coverage, window, mode='same')
    print(f"\\nSmoothed coverage around breakpoint:")
    print(f"  Positions 395-405: {smoothed_coverage[395:405]}")
    
    # Debug the detection loop
    print(f"\\nDetection loop debug:")
    print(f"  Loop range: {window_size} to {len(smoothed_coverage) - window_size}")
    print(f"  Sequence length: {len(smoothed_coverage)}")
    
    for i in range(max(390, window_size), min(410, len(smoothed_coverage) - window_size)):
        left_cov = np.mean(smoothed_coverage[i-window_size:i])
        right_cov = np.mean(smoothed_coverage[i:i+window_size])
        
        if left_cov > detector.min_coverage and right_cov > detector.min_coverage:
            fold_change = max(left_cov, right_cov) / min(left_cov, right_cov)
            print(f"  Position {i}: left={left_cov:.1f}, right={right_cov:.1f}, fold_change={fold_change:.2f} (threshold={detector.coverage_fold_change})")
            if fold_change >= detector.coverage_fold_change:
                print(f"    -> BREAKPOINT DETECTED at {i}")
        else:
            print(f"  Position {i}: left={left_cov:.1f}, right={right_cov:.1f} (coverage too low)")

if __name__ == "__main__":
    debug_coverage_detection()