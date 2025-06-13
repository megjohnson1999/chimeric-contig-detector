#!/usr/bin/env python3
"""
Example usage of Chimeric Detective for multi-sample processing.
"""

import os
import tempfile
from pathlib import Path

from chimeric_detective.multi_sample import MultiSampleProcessor, process_multi_sample_directory


def create_multi_sample_example_data():
    """Create example data with multiple samples."""
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="chimeric_detective_multi_sample_"))
    
    # Create example assembly
    assembly_file = temp_dir / "viral_assembly.fasta"
    with open(assembly_file, 'w') as f:
        # Normal contig
        f.write(">contig_1\n")
        f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 20 + "\n")
        
        # Chimeric contig - clear composition change
        f.write(">contig_2_chimeric\n")
        f.write("ATATATATATATATATATATATATATATATATATATATATAT" * 15)
        f.write("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" * 15 + "\n")
        
        # Another normal contig
        f.write(">contig_3\n")
        f.write("CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG" * 25 + "\n")
    
    # Create reads directory with multiple samples
    reads_dir = temp_dir / "reads"
    reads_dir.mkdir()
    
    # Sample 1
    with open(reads_dir / "sample1_R1.fastq.gz", 'w') as f:
        f.write("@read1\nATCGATCGATCGATCGATCGATCG\n+\n########################\n")
        f.write("@read2\nGCTAGCTAGCTAGCTAGCTAGCT\n+\n########################\n")
    
    with open(reads_dir / "sample1_R2.fastq.gz", 'w') as f:
        f.write("@read1\nCGATCGATCGATCGATCGATCGA\n+\n########################\n")
        f.write("@read2\nAGCTAGCTAGCTAGCTAGCTAG\n+\n########################\n")
    
    # Sample 2
    with open(reads_dir / "sample2_R1.fastq.gz", 'w') as f:
        f.write("@read1\nATATATATATATATATATATAT\n+\n########################\n")
        f.write("@read2\nGCGCGCGCGCGCGCGCGCGCGC\n+\n########################\n")
    
    with open(reads_dir / "sample2_R2.fastq.gz", 'w') as f:
        f.write("@read1\nATATATATATATATATATATAT\n+\n########################\n")
        f.write("@read2\nGCGCGCGCGCGCGCGCGCGCGC\n+\n########################\n")
    
    # Sample 3
    with open(reads_dir / "sample3_R1.fastq.gz", 'w') as f:
        f.write("@read1\nCGTAGCTAGCTAGCTAGCTAGC\n+\n########################\n")
        f.write("@read2\nTAGCTAGCTAGCTAGCTAGCTA\n+\n########################\n")
    
    with open(reads_dir / "sample3_R2.fastq.gz", 'w') as f:
        f.write("@read1\nGCTAGCTAGCTAGCTAGCTAGC\n+\n########################\n")
        f.write("@read2\nAGCTAGCTAGCTAGCTAGCTAG\n+\n########################\n")
    
    return str(assembly_file), str(reads_dir), str(temp_dir)


def run_separate_analysis_example():
    """Example: Analyze each sample separately."""
    
    print("🔬 Multi-Sample Example: Separate Analysis")
    print("=" * 60)
    
    # Create example data
    assembly_file, reads_dir, temp_dir = create_multi_sample_example_data()
    output_dir = Path(temp_dir) / "separate_analysis"
    
    print(f"📁 Assembly: {assembly_file}")
    print(f"📁 Reads directory: {reads_dir}")
    print(f"📁 Output directory: {output_dir}")
    
    # Process samples separately
    results = process_multi_sample_directory(
        assembly_file=assembly_file,
        reads_dir=reads_dir,
        reads_pattern="*_R{1,2}.fastq.gz",
        output_dir=str(output_dir),
        processing_mode="separate",
        max_workers=2,
        parallel=True,
        min_contig_length=500,
        min_coverage=1.0,
        log_level="INFO"
    )
    
    print(f"\n✅ Separate analysis complete!")
    print(f"📊 Processed {len(results)} samples:")
    for sample_name, sample_output in results.items():
        print(f"   - {sample_name}: {sample_output}")
    
    return str(output_dir)


def run_merged_analysis_example():
    """Example: Merge all samples and analyze together."""
    
    print("\n🔬 Multi-Sample Example: Merged Analysis")
    print("=" * 60)
    
    # Create example data
    assembly_file, reads_dir, temp_dir = create_multi_sample_example_data()
    output_dir = Path(temp_dir) / "merged_analysis"
    
    # Process samples as merged dataset
    results = process_multi_sample_directory(
        assembly_file=assembly_file,
        reads_dir=reads_dir,
        reads_pattern="*_R{1,2}.fastq.gz",
        output_dir=str(output_dir),
        processing_mode="merged",
        min_contig_length=500,
        min_coverage=1.0,
        log_level="INFO"
    )
    
    print(f"\n✅ Merged analysis complete!")
    print(f"📊 Merged results: {results}")
    
    return str(output_dir)


def run_batch_analysis_example():
    """Example: Process samples in batches."""
    
    print("\n🔬 Multi-Sample Example: Batch Processing")
    print("=" * 60)
    
    # Create example data
    assembly_file, reads_dir, temp_dir = create_multi_sample_example_data()
    output_dir = Path(temp_dir) / "batch_analysis"
    
    # Process samples in batches
    processor = MultiSampleProcessor(
        processing_mode="batch",
        max_workers=2,
        log_level="INFO"
    )
    
    results = processor.process_samples_directory(
        assembly_file=assembly_file,
        reads_dir=reads_dir,
        reads_pattern="*_R{1,2}.fastq.gz",
        output_dir=str(output_dir),
        batch_size=2,  # Process 2 samples per batch
        min_contig_length=500,
        min_coverage=1.0
    )
    
    print(f"\n✅ Batch processing complete!")
    print(f"📊 Processed {len(results)} samples in batches")
    
    return str(output_dir)


def demonstrate_custom_patterns():
    """Demonstrate different file naming patterns."""
    
    print("\n🔬 Multi-Sample Example: Custom Patterns")
    print("=" * 60)
    
    # Create temp directory with different naming patterns
    temp_dir = Path(tempfile.mkdtemp(prefix="pattern_test_"))
    reads_dir = temp_dir / "reads"
    reads_dir.mkdir()
    
    # Different naming patterns
    patterns_and_files = [
        ("*_R{1,2}.fastq.gz", ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]),
        ("*_{1,2}.fq.gz", ["sample2_1.fq.gz", "sample2_2.fq.gz"]),
        ("*.R{1,2}.fastq", ["sample3.R1.fastq", "sample3.R2.fastq"]),
    ]
    
    print("📋 Supported file naming patterns:")
    for pattern, example_files in patterns_and_files:
        print(f"   Pattern: {pattern}")
        print(f"   Example: {', '.join(example_files)}")
        
        # Create example files
        for filename in example_files:
            filepath = reads_dir / filename
            with open(filepath, 'w') as f:
                f.write("@read1\nATCG\n+\n####\n")
    
    print(f"\n📁 Created example files in: {reads_dir}")
    
    # Test pattern discovery
    from chimeric_detective.utils import parse_reads_pattern
    
    for pattern, _ in patterns_and_files:
        try:
            read_pairs = parse_reads_pattern(str(reads_dir), pattern)
            print(f"\n✅ Pattern '{pattern}' found {len(read_pairs)} sample(s):")
            for reads1, reads2 in read_pairs:
                print(f"   - R1: {Path(reads1).name}, R2: {Path(reads2).name if reads2 else 'None'}")
        except Exception as e:
            print(f"❌ Pattern '{pattern}' failed: {e}")


def main():
    """Run all multi-sample examples."""
    
    print("🧬 Chimeric Detective Multi-Sample Processing Examples")
    print("=" * 80)
    print()
    print("This script demonstrates different ways to process multiple samples:")
    print("1. Separate analysis - Each sample analyzed independently")
    print("2. Merged analysis - All reads combined for higher coverage")
    print("3. Batch processing - Memory-efficient processing of many samples")
    print("4. Custom file patterns - Different naming conventions")
    print()
    
    try:
        # Run examples
        separate_dir = run_separate_analysis_example()
        merged_dir = run_merged_analysis_example()
        batch_dir = run_batch_analysis_example()
        demonstrate_custom_patterns()
        
        print("\n" + "=" * 80)
        print("📋 SUMMARY")
        print("=" * 80)
        print(f"✅ All examples completed successfully!")
        print()
        print("📁 Example outputs:")
        print(f"   - Separate analysis: {separate_dir}")
        print(f"   - Merged analysis: {merged_dir}")
        print(f"   - Batch processing: {batch_dir}")
        print()
        print("💡 Command-line equivalents:")
        print()
        print("# Separate analysis (default)")
        print("chimeric_detective -a assembly.fasta --reads-dir reads/ -o results/")
        print()
        print("# Merged analysis")
        print("chimeric_detective -a assembly.fasta --reads-dir reads/ \\")
        print("                  --multi-sample-mode merged -o results/")
        print()
        print("# Parallel processing with custom workers")
        print("chimeric_detective -a assembly.fasta --reads-dir reads/ \\")
        print("                  --max-workers 8 --parallel -o results/")
        print()
        print("# Custom pattern")
        print("chimeric_detective -a assembly.fasta --reads-dir reads/ \\")
        print("                  --reads-pattern '*_{1,2}.fq.gz' -o results/")
        print()
        print("🎯 Check the output directories for detailed results!")
        
    except Exception as e:
        print(f"❌ Example failed: {e}")
        raise


if __name__ == "__main__":
    main()