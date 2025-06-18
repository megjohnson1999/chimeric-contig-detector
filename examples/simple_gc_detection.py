#!/usr/bin/env python3
"""
Example script demonstrating the simplified GC-only chimera detection.
"""

import sys
from pathlib import Path

# Add parent directory to path if running as script
sys.path.insert(0, str(Path(__file__).parent.parent))

from chimeric_detective import (
    SimpleChimeraDetector,
    SimpleChimeraAnalyzer,
    SimpleChimeraResolver,
    SimpleChimeraVisualizer
)


def run_simple_detection(assembly_file, output_dir="simple_gc_output"):
    """Run simplified GC-only chimera detection pipeline."""
    
    print("Running Simplified GC-Only Chimera Detection")
    print("=" * 50)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Step 1: Detect chimeras using GC content
    print("\n1. Detecting chimeras based on GC content...")
    detector = SimpleChimeraDetector(
        min_contig_length=1000,
        gc_content_threshold=0.1,
        window_size=1000,
        step_size=500
    )
    
    candidates = detector.detect_chimeras(assembly_file)
    print(f"   Found {len(candidates)} chimera candidates")
    
    # Step 2: Analyze and classify candidates
    print("\n2. Analyzing candidates...")
    analyzer = SimpleChimeraAnalyzer(
        gc_difference_threshold=0.15,
        high_confidence_gc_threshold=0.25
    )
    
    classifications = analyzer.analyze_candidates(candidates)
    high_confidence = analyzer.filter_high_confidence(classifications, min_confidence=0.7)
    
    print(f"   Total candidates: {len(candidates)}")
    print(f"   Likely chimeras: {sum(1 for c in classifications if c.is_likely_chimera)}")
    print(f"   High confidence: {len(high_confidence)}")
    
    # Print some details
    if high_confidence:
        print("\n   High confidence chimeras:")
        for classification in high_confidence[:5]:  # Show first 5
            candidate = classification.candidate
            print(f"   - {candidate.contig_id}: breakpoint at {candidate.breakpoint}, "
                  f"GC diff = {candidate.gc_difference:.3f}, confidence = {classification.confidence:.3f}")
        if len(high_confidence) > 5:
            print(f"   ... and {len(high_confidence) - 5} more")
    
    # Step 3: Resolve chimeras (split contigs)
    if high_confidence:
        print("\n3. Resolving chimeras...")
        resolver = SimpleChimeraResolver(
            min_fragment_length=500,
            output_dir=output_dir
        )
        
        statistics = resolver.resolve_chimeras(
            assembly_file,
            classifications,
            output_prefix="resolved"
        )
        
        print(f"   Split {statistics['contigs_split']} contigs")
        print(f"   Created {statistics['fragments_created']} fragments")
        print(f"   Output: {statistics['output_file']}")
    
    # Step 4: Generate visualization
    print("\n4. Generating report...")
    visualizer = SimpleChimeraVisualizer(output_dir=output_dir)
    
    report_file = visualizer.generate_report(
        assembly_file,
        candidates,
        classifications,
        output_prefix="simple_gc_report"
    )
    
    print(f"   Report: {report_file}")
    
    print("\nDone!")
    return candidates, classifications


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python simple_gc_detection.py <assembly.fasta>")
        sys.exit(1)
    
    assembly_file = sys.argv[1]
    if not Path(assembly_file).exists():
        print(f"Error: Assembly file '{assembly_file}' not found")
        sys.exit(1)
    
    run_simple_detection(assembly_file)