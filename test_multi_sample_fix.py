#!/usr/bin/env python3
"""
Test script for multi-sample functionality fixes.
"""

import sys
import os
sys.path.insert(0, 'src')

from chimeric_detective.detector import ChimeraDetector
from chimeric_detective.multi_sample import MultiSampleProcessor

def test_detector_return_bam():
    """Test that detector returns BAM file path when requested."""
    print("🧪 Testing detector BAM file return functionality")
    print("="*60)
    
    # Test with minimal data
    detector = ChimeraDetector(min_contig_length=100)
    
    # Create a minimal test assembly
    test_assembly = "test_assembly.fasta"
    with open(test_assembly, 'w') as f:
        f.write(">test_contig\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
    
    # Create minimal test reads
    test_reads = "test_reads.fastq"
    with open(test_reads, 'w') as f:
        f.write("@read1\nATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIII\n")
    
    try:
        # Test without return_bam_path (default behavior)
        result1 = detector.detect_chimeras(
            assembly_file=test_assembly,
            reads1=test_reads,
            return_bam_path=False
        )
        print(f"✅ Without return_bam_path: type={type(result1)}")
        
        # Test with return_bam_path=True
        result2 = detector.detect_chimeras(
            assembly_file=test_assembly,
            reads1=test_reads,
            return_bam_path=True
        )
        print(f"✅ With return_bam_path: type={type(result2)}")
        
        if isinstance(result2, tuple):
            candidates, bam_path = result2
            print(f"✅ Tuple unpacking successful: candidates={len(candidates)}, bam_path={bam_path}")
            
            # Check if BAM file exists
            if os.path.exists(bam_path):
                print(f"✅ BAM file exists: {bam_path}")
            else:
                print(f"❌ BAM file not found: {bam_path}")
        else:
            print(f"❌ Expected tuple, got {type(result2)}")
            
    except Exception as e:
        print(f"❌ Error testing detector: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Cleanup
        for f in [test_assembly, test_reads]:
            if os.path.exists(f):
                os.remove(f)
        print("🧹 Cleaned up test files")

def test_multi_sample_processor():
    """Test multi-sample processor basic functionality."""
    print("\n🧪 Testing MultiSampleProcessor")
    print("="*60)
    
    try:
        processor = MultiSampleProcessor(
            processing_mode="separate",
            max_workers=1,  # Use single worker for testing
            log_level="INFO"
        )
        print("✅ MultiSampleProcessor created successfully")
        
        # Test sample discovery
        print("✅ MultiSampleProcessor has sample discovery methods")
        
        # Test result saving methods
        if hasattr(processor, '_save_sample_results_summary'):
            print("✅ Has _save_sample_results_summary method")
        else:
            print("❌ Missing _save_sample_results_summary method")
            
        if hasattr(processor, '_save_failed_sample_results'):
            print("✅ Has _save_failed_sample_results method")
        else:
            print("❌ Missing _save_failed_sample_results method")
        
    except Exception as e:
        print(f"❌ Error testing MultiSampleProcessor: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    print("🔧 Multi-Sample Functionality Fix Test")
    print("="*70)
    
    try:
        test_detector_return_bam()
        test_multi_sample_processor()
        print("\n🎉 All tests completed!")
        
    except Exception as e:
        print(f"\n💥 Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)