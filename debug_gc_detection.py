#!/usr/bin/env python3
"""
Debug GC breakpoint detection specifically.
"""

import sys
sys.path.insert(0, './src')

from chimeric_detective.detector import ChimeraDetector

def debug_gc_detection():
    """Debug why GC breakpoint detection isn't working."""
    print("=== DEBUGGING GC BREAKPOINT DETECTION ===")
    
    detector = ChimeraDetector(
        gc_content_threshold=0.05,  # Very low threshold
        log_level='DEBUG'
    )
    
    # Create clear GC difference: GC rich -> AT rich
    sequence = 'GCGCGCGC' * 50 + 'ATATATATAT' * 50  # 400bp + 500bp
    window_size = 200
    step_size = 25
    
    print(f"Sequence length: {len(sequence)}")
    print(f"First 50bp GC content: {(sequence[:50].count('G') + sequence[:50].count('C'))/50}")
    print(f"Last 50bp GC content: {(sequence[-50:].count('G') + sequence[-50:].count('C'))/50}")
    
    # Calculate GC profile
    gc_profile = detector._calculate_adaptive_gc_profile(sequence, window_size, step_size)
    print(f"\nGC profile: {gc_profile}")
    
    # Manually check GC differences
    print(f"\nGC profile differences:")
    for i in range(len(gc_profile) - 1):
        diff = abs(gc_profile[i+1] - gc_profile[i])
        print(f"  Windows {i}->{i+1}: {gc_profile[i]:.3f} -> {gc_profile[i+1]:.3f}, diff = {diff:.3f}")
    
    # Test the detection logic
    print(f"\nThreshold: {detector.gc_content_threshold}")
    breakpoints = detector._detect_gc_breakpoints_adaptive(gc_profile, sequence, step_size)
    print(f"Detected breakpoints: {breakpoints}")
    
    # Manual implementation of detection to debug
    print(f"\nManual detection debug:")
    gc_diffs = []
    for i in range(len(gc_profile) - 1):
        diff = abs(gc_profile[i+1] - gc_profile[i])
        gc_diffs.append(diff)
        print(f"  Diff {i}: {diff:.3f} (threshold: {detector.gc_content_threshold})")
    
    # Check peak detection logic
    print(f"\nPeak detection debug:")
    for i, diff in enumerate(gc_diffs):
        if diff >= detector.gc_content_threshold:
            print(f"  Candidate {i}: diff={diff:.3f}")
            
            # Check if it's a peak
            is_peak = True
            if i > 0 and gc_diffs[i-1] >= diff:
                is_peak = False
                print(f"    Not peak: previous diff {gc_diffs[i-1]:.3f} >= {diff:.3f}")
            if i < len(gc_diffs) - 1 and gc_diffs[i+1] >= diff:
                is_peak = False
                print(f"    Not peak: next diff {gc_diffs[i+1]:.3f} >= {diff:.3f}")
            
            if is_peak:
                position = (i + 1) * step_size
                print(f"    PEAK DETECTED at window boundary position {position}")
                if 0 < position < len(sequence):
                    print(f"    Position is valid (in range 0-{len(sequence)})")
                else:
                    print(f"    Position is INVALID (outside range 0-{len(sequence)})")

if __name__ == "__main__":
    debug_gc_detection()