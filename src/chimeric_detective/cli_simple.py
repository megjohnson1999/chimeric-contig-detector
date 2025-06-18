"""
Simplified command-line interface for chimeric detective using only GC-based detection.
"""

import argparse
import sys
from pathlib import Path

from .detector_simple import SimpleChimeraDetector
from .analyzer_simple import SimpleChimeraAnalyzer
from .resolver_simple import SimpleChimeraResolver
from .visualizer_simple import SimpleChimeraVisualizer


def main():
    """Main entry point for simplified chimeric detective."""
    parser = argparse.ArgumentParser(
        description="Simplified Chimeric Detective - GC-based chimera detection for viral assemblies",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "assembly",
        help="Path to assembly FASTA file"
    )
    
    # Optional arguments
    parser.add_argument(
        "-o", "--output-dir",
        default="chimeric_detective_output",
        help="Output directory (default: chimeric_detective_output)"
    )
    
    parser.add_argument(
        "--min-contig-length",
        type=int,
        default=1000,
        help="Minimum contig length to analyze (default: 1000)"
    )
    
    parser.add_argument(
        "--gc-threshold",
        type=float,
        default=0.1,
        help="Minimum GC content difference to consider (default: 0.1)"
    )
    
    parser.add_argument(
        "--window-size",
        type=int,
        default=1000,
        help="Window size for GC analysis (default: 1000)"
    )
    
    parser.add_argument(
        "--step-size",
        type=int,
        default=500,
        help="Step size for sliding window (default: 500)"
    )
    
    parser.add_argument(
        "--min-confidence",
        type=float,
        default=0.7,
        help="Minimum confidence for splitting (default: 0.7)"
    )
    
    parser.add_argument(
        "--min-fragment-length",
        type=int,
        default=500,
        help="Minimum fragment length after splitting (default: 500)"
    )
    
    parser.add_argument(
        "--no-split",
        action="store_true",
        help="Only detect chimeras, don't split contigs"
    )
    
    parser.add_argument(
        "--no-visualize",
        action="store_true",
        help="Skip visualization and report generation"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Set up output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Set logging level
    log_level = "DEBUG" if args.verbose else "INFO"
    
    try:
        # Step 1: Detect chimeras
        print("Step 1: Detecting chimeras using GC content analysis...")
        detector = SimpleChimeraDetector(
            min_contig_length=args.min_contig_length,
            gc_content_threshold=args.gc_threshold,
            window_size=args.window_size,
            step_size=args.step_size,
            log_level=log_level
        )
        
        candidates = detector.detect_chimeras(args.assembly)
        print(f"Found {len(candidates)} potential chimeras")
        
        # Step 2: Analyze and classify candidates
        print("\nStep 2: Analyzing chimera candidates...")
        analyzer = SimpleChimeraAnalyzer(
            gc_difference_threshold=args.gc_threshold,
            log_level=log_level
        )
        
        classifications = analyzer.analyze_candidates(candidates)
        high_confidence = analyzer.filter_high_confidence(classifications, args.min_confidence)
        print(f"Identified {len(high_confidence)} high-confidence chimeras")
        
        # Step 3: Resolve chimeras (if requested)
        if not args.no_split and high_confidence:
            print("\nStep 3: Resolving chimeras...")
            resolver = SimpleChimeraResolver(
                min_fragment_length=args.min_fragment_length,
                output_dir=output_dir,
                log_level=log_level
            )
            
            statistics = resolver.resolve_chimeras(
                args.assembly,
                classifications,
                output_prefix="resolved"
            )
            
            print(f"Split {statistics['contigs_split']} contigs into {statistics['fragments_created']} fragments")
            print(f"Resolved assembly written to: {statistics['output_file']}")
        
        # Step 4: Generate visualization (if requested)
        if not args.no_visualize and candidates:
            print("\nStep 4: Generating visualization...")
            visualizer = SimpleChimeraVisualizer(
                output_dir=output_dir,
                log_level=log_level
            )
            
            report_file = visualizer.generate_report(
                args.assembly,
                candidates,
                classifications,
                output_prefix="chimera_report"
            )
            
            print(f"Report generated: {report_file}")
        
        print("\nAnalysis complete!")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()