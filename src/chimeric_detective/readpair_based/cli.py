"""
Command-line interface for read-pair based chimera detection.
"""

import sys
import argparse
from pathlib import Path
import logging

from .config import DetectorConfig, create_argument_parser
from .core import ReadPairChimeraDetector, DetectorFactory


def main():
    """Main CLI entry point."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    try:
        # Load base configuration
        if args.config:
            config = DetectorConfig.from_file(args.config)
        else:
            config = DetectorConfig()
        
        # Apply command-line overrides
        _apply_cli_overrides(config, args)
        
        # Validate configuration
        issues = config.validate()
        if issues:
            print("Configuration issues found:", file=sys.stderr)
            for issue in issues:
                print(f"  - {issue}", file=sys.stderr)
            sys.exit(1)
        
        # Set up logging
        if args.debug:
            config.logging.level = "DEBUG"
        elif args.verbose:
            config.logging.level = "INFO"
        
        # Create detector
        detector = ReadPairChimeraDetector(config)
        
        # Run detection
        results = detector.detect_chimeras(args.bam, args.assembly)
        
        # Write results
        output_files = detector.write_results(args.output_dir)
        
        # Print summary
        print(f"\nDetection Summary:")
        print(f"  Total candidates: {len(results)}")
        print(f"  High confidence (≥0.8): {len(detector.get_high_confidence_candidates(0.8))}")
        print(f"  Medium confidence (0.5-0.8): {len([r for r in results if 0.5 <= r.confidence < 0.8])}")
        print(f"  Low confidence (<0.5): {len([r for r in results if r.confidence < 0.5])}")
        
        print(f"\nOutput files:")
        for format_name, filepath in output_files.items():
            print(f"  {format_name}: {filepath}")
        
        if not results:
            print("\nNo chimeric breakpoints detected.")
            return 0
        
        # Show top candidates
        print(f"\nTop 5 candidates:")
        for i, candidate in enumerate(sorted(results, key=lambda x: x.confidence, reverse=True)[:5]):
            print(f"  {i+1}. {candidate.contig}:{candidate.position} (confidence: {candidate.confidence:.3f})")
        
        return 0
        
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if hasattr(args, 'debug') and args.debug:
            import traceback
            traceback.print_exc()
        return 1


def _apply_cli_overrides(config: DetectorConfig, args: argparse.Namespace):
    """Apply command-line argument overrides to config."""
    # Read quality overrides
    if args.min_mapping_quality is not None:
        config.read_quality.min_mapping_quality = args.min_mapping_quality
    
    # Window size overrides
    if args.window_size is not None:
        config.sliding_window.window_size = args.window_size
    if args.step_size is not None:
        config.sliding_window.step_size = args.step_size
    
    # Insert size method override
    if args.outlier_method is not None:
        config.insert_size.outlier_method = args.outlier_method
    
    # Threading override
    if args.threads is not None:
        config.threads = args.threads
    
    # Output format overrides
    if args.output_formats is not None:
        config.output.formats = args.output_formats


def validate_command():
    """Validation subcommand for checking inputs."""
    parser = argparse.ArgumentParser(
        description="Validate BAM file for chimera detection"
    )
    parser.add_argument("bam", help="BAM file to validate")
    parser.add_argument("-v", "--verbose", action="store_true")
    
    args = parser.parse_args(sys.argv[2:])  # Skip 'validate' subcommand
    
    from .bamparser import BamValidator
    
    print(f"Validating BAM file: {args.bam}")
    issues = BamValidator.validate_bam_file(args.bam)
    
    if not issues:
        print("✓ BAM file validation passed")
        return 0
    else:
        print("⚠ BAM file validation issues found:")
        for issue in issues:
            print(f"  - {issue}")
        return 1


def config_command():
    """Configuration subcommand for creating config templates."""
    parser = argparse.ArgumentParser(
        description="Generate configuration file templates"
    )
    parser.add_argument("template", choices=["default", "sensitive", "specific"],
                       help="Configuration template type")
    parser.add_argument("-o", "--output", default="chimera_config.yaml",
                       help="Output configuration file")
    parser.add_argument("-f", "--format", choices=["yaml", "json"], default="yaml",
                       help="Configuration file format")
    
    args = parser.parse_args(sys.argv[2:])  # Skip 'config' subcommand
    
    # Create appropriate configuration
    if args.template == "default":
        config = DetectorConfig()
    elif args.template == "sensitive":
        detector = DetectorFactory.create_sensitive()
        config = detector.config
    elif args.template == "specific":
        detector = DetectorFactory.create_specific()
        config = detector.config
    
    # Save configuration
    config.save(args.output, args.format)
    print(f"Configuration template saved to {args.output}")
    
    return 0


def main_with_subcommands():
    """Main entry point with subcommands."""
    if len(sys.argv) < 2:
        main()
        return
    
    subcommand = sys.argv[1]
    
    if subcommand == "validate":
        return validate_command()
    elif subcommand == "config":
        return config_command()
    else:
        # Treat as main command
        main()


if __name__ == "__main__":
    sys.exit(main())