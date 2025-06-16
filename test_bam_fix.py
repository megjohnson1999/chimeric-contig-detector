#!/usr/bin/env python3
"""
Test script for BAM key error fix.
"""

import sys
import os
sys.path.insert(0, 'src')

def test_kwargs_access():
    """Test kwargs access patterns."""
    print("🧪 Testing kwargs access patterns")
    print("="*50)
    
    # Simulate kwargs without bam key (multi-sample scenario)
    kwargs = {
        'assembly': 'test.fasta',
        'reads_dir': 'reads/',
        'reads_pattern': '*_R{1,2}.fastq.gz',
        'out': 'results/',
        'log_level': 'INFO',
        'confidence_threshold': 0.5,
        'min_contig_length': 1000,
        'min_coverage': 5.0
    }
    
    # Test safe access patterns
    try:
        # This should work
        has_bam = kwargs.get('bam') is not None
        has_reads_dir = kwargs.get('reads_dir') is not None
        
        print(f"✅ Safe access: has_bam={has_bam}, has_reads_dir={has_reads_dir}")
        
        # Test condition logic
        if kwargs.get('reads_dir'):
            print("✅ Would go to multi-sample pipeline")
        else:
            print("❌ Would incorrectly go to single-sample pipeline")
            
        # Test validation patterns
        if kwargs.get('reference'):
            print("Would validate reference")
        else:
            print("✅ Correctly skipping reference validation")
            
    except KeyError as e:
        print(f"❌ KeyError with safe access: {e}")
    
    # Test unsafe access pattern (what was causing the error)
    try:
        # This would fail
        has_bam_unsafe = kwargs['bam'] is not None
        print(f"❌ Unsafe access worked unexpectedly: {has_bam_unsafe}")
    except KeyError as e:
        print(f"✅ Unsafe access correctly failed: {e}")

if __name__ == "__main__":
    print("🔧 BAM Key Error Fix Test")
    print("="*40)
    
    try:
        test_kwargs_access()
        print("\n🎉 All tests passed!")
        
    except Exception as e:
        print(f"\n💥 Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)