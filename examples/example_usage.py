#!/usr/bin/env python3
"""
Example usage of Chimeric Detective programmatically.
"""

import os
import tempfile
from pathlib import Path

from chimeric_detective.detector import ChimeraDetector
from chimeric_detective.analyzer import ChimeraAnalyzer
from chimeric_detective.resolver import ChimeraResolver
from chimeric_detective.visualizer import ChimeraVisualizer


def create_example_data():
    """Create example assembly and read files for testing."""
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="chimeric_detective_example_"))
    
    # Create example assembly with a chimeric contig
    assembly_file = temp_dir / "example_assembly.fasta"
    with open(assembly_file, 'w') as f:
        # Normal contig
        f.write(">contig_1\n")
        f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 20 + "\n")
        
        # Chimeric contig - clear composition change in the middle
        f.write(">contig_2_chimeric\n")
        # First part: AT-rich
        f.write("ATATATATATATATATATATATATATATATATATATATATAT" * 15)
        # Second part: GC-rich  
        f.write("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" * 15 + "\n")
        
        # Another normal contig
        f.write(">contig_3\n")
        f.write("CGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG" * 25 + "\n")
    
    # Create example reads (simplified - would normally be FASTQ)
    reads_file = temp_dir / "example_reads.fasta"
    with open(reads_file, 'w') as f:
        # Reads covering the normal regions
        for i in range(100):
            f.write(f">read_{i}\n")
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
    
    return str(assembly_file), str(reads_file), str(temp_dir)


def run_example():
    """Run complete example analysis."""
    
    print("🔬 Chimeric Detective Example Analysis")
    print("=" * 50)
    
    # Create example data
    print("📁 Creating example data...")
    assembly_file, reads_file, temp_dir = create_example_data()
    
    output_dir = Path(temp_dir) / "results"
    output_dir.mkdir(exist_ok=True)
    
    print(f"Assembly file: {assembly_file}")
    print(f"Reads file: {reads_file}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Step 1: Chimera Detection
        print("\n🔍 Step 1: Detecting chimeric contigs...")
        
        detector = ChimeraDetector(
            min_contig_length=500,
            min_coverage=1.0,  # Low threshold for example
            coverage_fold_change=1.5,  # Sensitive settings
            gc_content_threshold=0.05,
            kmer_distance_threshold=0.2,
            log_level="INFO"
        )
        
        # For this example, we'll simulate without actual BAM alignment
        # In real usage, you would provide actual reads
        candidates = []
        
        # Manually create a candidate for demonstration
        from chimeric_detective.detector import ChimeraCandidate
        
        example_candidate = ChimeraCandidate(
            contig_id="contig_2_chimeric",
            breakpoint=630,  # Middle of the sequence
            confidence_score=0.85,
            evidence_types=["gc_content_shift", "kmer_composition_change"],
            coverage_left=25.0,
            coverage_right=28.0,
            gc_content_left=0.0,  # AT-rich
            gc_content_right=1.0,  # GC-rich
            kmer_distance=0.8,
            spanning_reads=3,
            read_orientation_score=0.1
        )
        
        candidates = [example_candidate]
        
        print(f"✅ Detected {len(candidates)} chimera candidates")
        
        # Step 2: Chimera Analysis
        print("\n🧬 Step 2: Analyzing and classifying chimeras...")
        
        analyzer = ChimeraAnalyzer(
            reference_db=None,  # No reference database for this example
            log_level="INFO"
        )
        
        # For demonstration, manually create analysis
        from chimeric_detective.analyzer import ChimeraAnalysis
        
        example_analysis = ChimeraAnalysis(
            candidate=example_candidate,
            chimera_type="technical_artifact",
            classification_confidence=0.82,
            explanation=(
                "Contig contig_2_chimeric shows clear evidence of being a technical chimera "
                "resulting from misassembly. At position 630, there is a dramatic shift in "
                "sequence composition from AT-rich (0% GC) to GC-rich (100% GC). The k-mer "
                "distance of 0.8 indicates completely different sequence composition. This "
                "extreme composition change suggests incorrect joining of two unrelated "
                "sequences during assembly. Confidence: 0.82."
            ),
            recommendation="SPLIT - High confidence technical artifact, recommend splitting",
            taxonomic_left="Unknown",
            taxonomic_right="Unknown",
            biological_evidence={
                'taxonomic_similarity': 0.0,
                'junction_support': {'spanning_pairs': 3, 'split_reads': 0},
                'coverage_gradient': {'gradient_type': 'moderate', 'fold_change': 1.12}
            },
            technical_evidence={
                'coverage_ratio': 1.12,
                'spanning_reads_ratio': 0.11,
                'kmer_distance': 0.8,
                'gc_content_diff': 1.0,
                'orientation_score': 0.1
            }
        )
        
        analyses = [example_analysis]
        
        print(f"✅ Analyzed {len(analyses)} chimeras")
        for analysis in analyses:
            print(f"   - {analysis.candidate.contig_id}: {analysis.chimera_type} "
                  f"(confidence: {analysis.classification_confidence:.3f})")
        
        # Step 3: Chimera Resolution
        print("\n✂️ Step 3: Resolving chimeras and creating cleaned assembly...")
        
        resolver = ChimeraResolver(
            split_technical=True,
            split_pcr=True,
            preserve_biological=True,
            min_split_length=200,
            confidence_threshold=0.5,
            log_level="INFO"
        )
        
        output_files = resolver.resolve_chimeras(
            analyses=analyses,
            assembly_file=assembly_file,
            output_dir=str(output_dir)
        )
        
        print(f"✅ Resolution complete")
        print(f"   - Split contigs: {len([d for d in resolver.splitting_decisions if d.action == 'split'])}")
        print(f"   - Preserved contigs: {len([d for d in resolver.splitting_decisions if d.action == 'preserve'])}")
        
        # Step 4: Generate Report
        print("\n📊 Step 4: Generating interactive report...")
        
        visualizer = ChimeraVisualizer(log_level="INFO")
        
        report_path = visualizer.create_report(
            analyses=analyses,
            decisions=resolver.splitting_decisions,
            output_dir=str(output_dir),
            assembly_file=assembly_file
        )
        
        print(f"✅ Interactive report generated: {report_path}")
        
        # Print summary
        print("\n" + "=" * 50)
        print("📋 ANALYSIS SUMMARY")
        print("=" * 50)
        
        print(f"📁 Output files:")
        for file_type, file_path in output_files.items():
            print(f"   - {file_type}: {file_path}")
        
        print(f"\n🎯 Results:")
        print(f"   - Total contigs analyzed: {len(analyses)}")
        print(f"   - Chimeric contigs found: {len([a for a in analyses])}")
        print(f"   - Technical artifacts: {len([a for a in analyses if a.chimera_type == 'technical_artifact'])}")
        print(f"   - Biological recombination: {len([a for a in analyses if a.chimera_type == 'biological_recombination'])}")
        
        print(f"\n📊 Open the interactive report to explore results:")
        print(f"   file://{report_path}")
        
        print(f"\n🗂️ All files are in: {output_dir}")
        
    except Exception as e:
        print(f"❌ Error during analysis: {str(e)}")
        raise
    
    return str(output_dir)


if __name__ == "__main__":
    output_dir = run_example()
    print(f"\n✅ Example completed successfully!")
    print(f"📁 Check the results in: {output_dir}")