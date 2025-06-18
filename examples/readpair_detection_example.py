#!/usr/bin/env python3
"""
Example script demonstrating read-pair based chimera detection.

This example shows how to use the new read-pair focused API for detecting
chimeric contigs in viral metagenomic assemblies.
"""

import sys
from pathlib import Path

# Add parent directory to path if running as script
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chimeric_detective.readpair_based import (
    ReadPairChimeraDetector,
    DetectorConfig,
    DetectorFactory
)


def basic_example(bam_file: str, output_dir: str = "readpair_output"):
    """Basic usage example with default configuration."""
    print("=== Basic Read-Pair Chimera Detection ===")
    
    # Create detector with default settings
    detector = DetectorFactory.create_default()
    
    # Run detection
    print(f"Analyzing BAM file: {bam_file}")
    candidates = detector.detect_chimeras(bam_file)
    
    # Write results
    output_files = detector.write_results(output_dir)
    
    # Print summary
    print(f"\nResults:")
    print(f"  Total candidates: {len(candidates)}")
    print(f"  High confidence: {len(detector.get_high_confidence_candidates())}")
    
    print(f"\nOutput files:")
    for fmt, path in output_files.items():
        print(f"  {fmt}: {path}")
    
    return candidates


def sensitive_example(bam_file: str, output_dir: str = "readpair_sensitive_output"):
    """Example using sensitive detection parameters."""
    print("\n=== Sensitive Detection Mode ===")
    
    # Create detector optimized for sensitivity
    detector = DetectorFactory.create_sensitive()
    
    # Run detection
    candidates = detector.detect_chimeras(bam_file)
    
    # Write results
    output_files = detector.write_results(output_dir)
    
    print(f"Sensitive mode found {len(candidates)} candidates")
    return candidates


def custom_config_example(bam_file: str, output_dir: str = "readpair_custom_output"):
    """Example with custom configuration."""
    print("\n=== Custom Configuration Example ===")
    
    # Create custom configuration
    config = DetectorConfig()
    
    # Customize read quality filters
    config.read_quality.min_mapping_quality = 30
    config.read_quality.require_proper_pairs = False  # Include discordant pairs
    
    # Customize window analysis
    config.sliding_window.window_size = 2000
    config.sliding_window.step_size = 1000
    config.sliding_window.min_pairs_per_window = 20
    
    # Customize detection thresholds
    config.detection.proper_pair_drop_threshold = 0.4
    config.detection.insert_size_shift_threshold = 2.5
    config.detection.discordant_pair_threshold = 0.25
    
    # Customize output
    config.output.formats = ["json", "tsv", "bed"]
    config.output.include_debug_info = True
    
    # Enable debug logging
    config.logging.level = "DEBUG"
    
    # Create detector with custom config
    detector = ReadPairChimeraDetector(config)
    
    # Run detection
    candidates = detector.detect_chimeras(bam_file)
    
    # Write results
    output_files = detector.write_results(output_dir)
    
    print(f"Custom configuration found {len(candidates)} candidates")
    return candidates


def analysis_example(bam_file: str):
    """Example showing detailed analysis of results."""
    print("\n=== Detailed Analysis Example ===")
    
    # Use default detector
    detector = DetectorFactory.create_default()
    candidates = detector.detect_chimeras(bam_file)
    
    if not candidates:
        print("No candidates found for analysis")
        return
    
    # Analyze by contig
    by_contig = detector.get_candidates_by_contig()
    print(f"\nCandidates by contig:")
    for contig, contig_candidates in by_contig.items():
        print(f"  {contig}: {len(contig_candidates)} candidates")
    
    # Show top candidates
    print(f"\nTop 3 candidates:")
    for i, candidate in enumerate(sorted(candidates, key=lambda x: x.confidence, reverse=True)[:3]):
        print(f"\n  {i+1}. {candidate.contig}:{candidate.position}")
        print(f"     Confidence: {candidate.confidence:.3f}")
        print(f"     Anomaly types: {', '.join(set(a.anomaly_type for a in candidate.supporting_anomalies))}")
        print(f"     Left proper pair ratio: {candidate.left_metrics.proper_pair_ratio:.3f}")
        print(f"     Right proper pair ratio: {candidate.right_metrics.proper_pair_ratio:.3f}")
        print(f"     Left insert size median: {candidate.left_metrics.insert_size_median:.1f}")
        print(f"     Right insert size median: {candidate.right_metrics.insert_size_median:.1f}")


def config_file_example(config_file: str, bam_file: str, output_dir: str = "readpair_config_output"):
    """Example using configuration file."""
    print(f"\n=== Configuration File Example ===")
    
    # Create detector from config file
    detector = DetectorFactory.create_from_file(config_file)
    
    # Run detection
    candidates = detector.detect_chimeras(bam_file)
    
    # Write results
    output_files = detector.write_results(output_dir)
    
    print(f"Config file mode found {len(candidates)} candidates")
    return candidates


def create_example_config():
    """Create an example configuration file."""
    print("\n=== Creating Example Configuration ===")
    
    config = DetectorConfig()
    
    # Save as YAML
    config.save("example_config.yaml", "yaml")
    print("Example configuration saved to: example_config.yaml")
    
    # Save as JSON
    config.save("example_config.json", "json")
    print("Example configuration saved to: example_config.json")


def main():
    """Main example function."""
    if len(sys.argv) < 2:
        print("Usage: python readpair_detection_example.py <bam_file> [assembly_file]")
        print("\nThis script demonstrates various ways to use the read-pair based chimera detector.")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    assembly_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not Path(bam_file).exists():
        print(f"Error: BAM file not found: {bam_file}")
        sys.exit(1)
    
    print(f"Read-Pair Based Chimera Detection Examples")
    print(f"BAM file: {bam_file}")
    if assembly_file:
        print(f"Assembly file: {assembly_file}")
    print("=" * 50)
    
    try:
        # Run different examples
        candidates1 = basic_example(bam_file)
        candidates2 = sensitive_example(bam_file)
        candidates3 = custom_config_example(bam_file)
        
        # Detailed analysis
        analysis_example(bam_file)
        
        # Create example configs
        create_example_config()
        
        # Compare results
        print(f"\n=== Comparison ===")
        print(f"Default mode: {len(candidates1)} candidates")
        print(f"Sensitive mode: {len(candidates2)} candidates")
        print(f"Custom mode: {len(candidates3)} candidates")
        
        print(f"\nExample complete! Check the output directories for results.")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()