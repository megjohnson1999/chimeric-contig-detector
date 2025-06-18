#!/usr/bin/env python3
"""
Debug script to reproduce the exact range() error and find its source.
"""

import sys
import traceback
import tempfile
from pathlib import Path

# Add src to path so we can import
sys.path.insert(0, str(Path(__file__).parent / "src"))

from chimeric_detective.detector import ChimeraDetector

def create_test_data():
    """Create minimal test data to reproduce the error."""
    # Create a simple test assembly
    assembly_content = """>contig1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
"""
    
    # Create test reads (simple paired-end)
    reads1_content = """@read1/1
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
@read2/1  
GCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIII
"""
    
    reads2_content = """@read1/2
CGATCGATCGATCGATCGAT
+
IIIIIIIIIIIIIIIIIIII
@read2/2
AGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIII
"""
    
    # Write to temp files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(assembly_content)
        assembly_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='_R1.fastq', delete=False) as f:
        f.write(reads1_content)
        reads1_file = f.name
        
    with tempfile.NamedTemporaryFile(mode='w', suffix='_R2.fastq', delete=False) as f:
        f.write(reads2_content)
        reads2_file = f.name
    
    return assembly_file, reads1_file, reads2_file

def test_detector_step_by_step():
    """Test each detector method step by step to find where the error occurs."""
    print("Creating test data...")
    assembly_file, reads1_file, reads2_file = create_test_data()
    
    print(f"Assembly: {assembly_file}")
    print(f"Reads1: {reads1_file}")  
    print(f"Reads2: {reads2_file}")
    
    try:
        print("\n1. Creating ChimeraDetector...")
        detector = ChimeraDetector(
            min_contig_length=50,  # Very small for testing
            min_coverage=1.0,
            window_size=20,        # Small window
            step_size=10,
            log_level="DEBUG"
        )
        print("✓ ChimeraDetector created")
        
        print("\n2. Testing detect_chimeras...")
        candidates = detector.detect_chimeras(
            assembly_file=assembly_file,
            reads1=reads1_file, 
            reads2=reads2_file
        )
        print(f"✓ Detection completed. Found {len(candidates)} candidates")
        
    except Exception as e:
        print(f"\n❌ ERROR occurred: {e}")
        print(f"Error type: {type(e).__name__}")
        print("\nFull traceback:")
        traceback.print_exc()
        
        # Try to extract the exact line that failed
        tb = traceback.extract_tb(e.__traceback__)
        for frame in tb:
            if 'range(' in frame.line:
                print(f"\n🎯 FOUND range() error at:")
                print(f"File: {frame.filename}")
                print(f"Line {frame.lineno}: {frame.line}")
                print(f"Function: {frame.name}")
                break
    
    finally:
        # Cleanup
        import os
        for f in [assembly_file, reads1_file, reads2_file]:
            try:
                os.unlink(f)
            except:
                pass

if __name__ == "__main__":
    test_detector_step_by_step()