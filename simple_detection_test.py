#!/usr/bin/env python3
"""
Simple test to verify chimera detection without full pipeline.
"""

import sys
sys.path.append('chimeric_detective/src')

from chimeric_detective.detector import ChimeraDetector

def test_detection():
    """Test chimera detection on our synthetic dataset."""
    print("Testing chimera detection...")
    
    detector = ChimeraDetector(
        min_contig_length=1000,
        min_coverage=5.0,
        coverage_fold_change=2.0,
        gc_content_threshold=0.1,
        kmer_distance_threshold=0.3
    )
    
    assembly_file = "test_data_large/large_test_assembly.fasta"
    reads1 = "test_data_large/reads/sample_001_R1.fastq.gz"
    reads2 = "test_data_large/reads/sample_001_R2.fastq.gz"
    
    try:
        candidates = detector.detect_chimeras(
            assembly_file=assembly_file,
            reads1=reads1,
            reads2=reads2
        )
        
        print(f"Detection successful! Found {len(candidates)} chimera candidates")
        
        for i, candidate in enumerate(candidates[:5]):  # Show first 5
            print(f"\nCandidate {i+1}:")
            print(f"  Contig: {candidate.contig_id}")
            print(f"  Breakpoint: {candidate.breakpoint}")
            print(f"  Confidence: {candidate.confidence_score:.3f}")
            print(f"  Evidence: {', '.join(candidate.evidence_types)}")
            print(f"  Coverage L/R: {candidate.coverage_left:.1f}/{candidate.coverage_right:.1f}")
        
        # Check if we detected our known chimeric contigs
        detected_contigs = {c.contig_id for c in candidates}
        expected_chimeric = ['chimeric_001', 'chimeric_002', 'chimeric_003', 'chimeric_004', 'chimeric_005',
                           'cov_chimeric_001', 'cov_chimeric_002', 'cov_chimeric_003']
        
        found_expected = detected_contigs.intersection(expected_chimeric)
        print(f"\nExpected chimeric contigs found: {len(found_expected)}/{len(expected_chimeric)}")
        print(f"Found: {sorted(found_expected)}")
        
        return True
        
    except Exception as e:
        print(f"Detection failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_detection()