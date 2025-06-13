"""
Multi-sample processing module for handling directories with multiple samples.
"""

import os
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import json
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from .detector import ChimeraDetector, ChimeraCandidate
from .analyzer import ChimeraAnalyzer, ChimeraAnalysis
from .resolver import ChimeraResolver, SplittingDecision
from .visualizer import ChimeraVisualizer
from .utils import parse_reads_pattern, setup_logging, create_output_directory


class MultiSampleProcessor:
    """Processor for handling multiple samples in a directory."""
    
    def __init__(self,
                 processing_mode: str = "separate",  # "separate", "merged", "batch"
                 max_workers: int = 4,
                 log_level: str = "INFO"):
        """
        Initialize MultiSampleProcessor.
        
        Args:
            processing_mode: How to process multiple samples
                - "separate": Process each sample independently
                - "merged": Merge all reads and process as one sample
                - "batch": Process in batches for memory efficiency
            max_workers: Maximum number of parallel workers
            log_level: Logging level
        """
        self.processing_mode = processing_mode
        self.max_workers = max_workers
        self.logger = setup_logging(log_level)
        
    def process_samples_directory(self,
                                 assembly_file: str,
                                 reads_dir: str,
                                 reads_pattern: str,
                                 output_dir: str,
                                 **kwargs) -> Dict[str, str]:
        """
        Process all samples in a directory.
        
        Args:
            assembly_file: Path to assembly FASTA file
            reads_dir: Directory containing read files
            reads_pattern: Pattern for finding read files
            output_dir: Output directory for results
            **kwargs: Additional parameters for detection/analysis
            
        Returns:
            Dictionary mapping sample names to their output directories
        """
        self.logger.info(f"Processing samples from directory: {reads_dir}")
        self.logger.info(f"Using pattern: {reads_pattern}")
        self.logger.info(f"Processing mode: {self.processing_mode}")
        
        # Find all sample read files
        sample_files = self._discover_samples(reads_dir, reads_pattern)
        
        if not sample_files:
            raise ValueError(f"No samples found matching pattern: {reads_pattern}")
        
        self.logger.info(f"Found {len(sample_files)} samples")
        
        # Create output directory structure
        create_output_directory(output_dir)
        
        if self.processing_mode == "separate":
            return self._process_samples_separately(
                assembly_file, sample_files, output_dir, **kwargs
            )
        elif self.processing_mode == "merged":
            return self._process_samples_merged(
                assembly_file, sample_files, output_dir, **kwargs
            )
        elif self.processing_mode == "batch":
            return self._process_samples_batch(
                assembly_file, sample_files, output_dir, **kwargs
            )
        else:
            raise ValueError(f"Unknown processing mode: {self.processing_mode}")
    
    def _discover_samples(self, reads_dir: str, reads_pattern: str) -> Dict[str, Tuple[str, Optional[str]]]:
        """Discover samples in the reads directory."""
        
        read_pairs = parse_reads_pattern(reads_dir, reads_pattern)
        
        samples = {}
        for reads1, reads2 in read_pairs:
            # Extract sample name from filename
            sample_name = self._extract_sample_name(reads1, reads_pattern)
            samples[sample_name] = (reads1, reads2)
        
        return samples
    
    def _extract_sample_name(self, reads1_path: str, pattern: str) -> str:
        """Extract sample name from read file path."""
        filename = Path(reads1_path).name
        
        # Remove common read file suffixes
        for suffix in ['_R1.fastq.gz', '_R1.fastq', '_R1.fq.gz', '_R1.fq',
                      '_1.fastq.gz', '_1.fastq', '_1.fq.gz', '_1.fq',
                      '.fastq.gz', '.fastq', '.fq.gz', '.fq']:
            if filename.endswith(suffix):
                return filename[:-len(suffix)]
        
        # If no standard suffix found, use the whole filename without extension
        return Path(filename).stem
    
    def _process_samples_separately(self,
                                  assembly_file: str,
                                  sample_files: Dict[str, Tuple[str, Optional[str]]],
                                  output_dir: str,
                                  **kwargs) -> Dict[str, str]:
        """Process each sample separately."""
        
        results = {}
        sample_outputs = []
        
        if kwargs.get('parallel', True) and len(sample_files) > 1:
            # Parallel processing
            self.logger.info(f"Processing {len(sample_files)} samples in parallel")
            
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {}
                
                for sample_name, (reads1, reads2) in sample_files.items():
                    sample_output_dir = Path(output_dir) / f"sample_{sample_name}"
                    
                    future = executor.submit(
                        self._process_single_sample,
                        assembly_file, reads1, reads2, sample_name, 
                        str(sample_output_dir), **kwargs
                    )
                    futures[future] = sample_name
                
                # Collect results with progress bar
                for future in tqdm(as_completed(futures), total=len(futures), 
                                 desc="Processing samples"):
                    sample_name = futures[future]
                    try:
                        sample_result = future.result()
                        results[sample_name] = sample_result
                        sample_outputs.append(sample_result)
                    except Exception as e:
                        self.logger.error(f"Failed to process sample {sample_name}: {e}")
        
        else:
            # Sequential processing
            self.logger.info(f"Processing {len(sample_files)} samples sequentially")
            
            for sample_name, (reads1, reads2) in tqdm(sample_files.items(), 
                                                    desc="Processing samples"):
                sample_output_dir = Path(output_dir) / f"sample_{sample_name}"
                
                try:
                    sample_result = self._process_single_sample(
                        assembly_file, reads1, reads2, sample_name,
                        str(sample_output_dir), **kwargs
                    )
                    results[sample_name] = sample_result
                    sample_outputs.append(sample_result)
                except Exception as e:
                    self.logger.error(f"Failed to process sample {sample_name}: {e}")
        
        # Generate combined report
        self._generate_multi_sample_report(sample_outputs, output_dir)
        
        return results
    
    def _process_single_sample(self,
                             assembly_file: str,
                             reads1: str,
                             reads2: Optional[str],
                             sample_name: str,
                             sample_output_dir: str,
                             **kwargs) -> str:
        """Process a single sample."""
        
        create_output_directory(sample_output_dir)
        
        # Initialize components
        detector = ChimeraDetector(
            min_contig_length=kwargs.get('min_contig_length', 1000),
            min_coverage=kwargs.get('min_coverage', 5.0),
            coverage_fold_change=kwargs.get('coverage_fold_change', 2.0),
            gc_content_threshold=kwargs.get('gc_content_threshold', 0.1),
            kmer_distance_threshold=kwargs.get('kmer_distance_threshold', 0.3),
            log_level=kwargs.get('log_level', 'INFO')
        )
        
        analyzer = ChimeraAnalyzer(
            reference_db=kwargs.get('reference'),
            log_level=kwargs.get('log_level', 'INFO')
        )
        
        resolver = ChimeraResolver(
            split_technical=kwargs.get('split_technical', True),
            split_pcr=kwargs.get('split_pcr', True),
            preserve_biological=kwargs.get('preserve_biological', True),
            min_split_length=kwargs.get('min_split_length', 500),
            confidence_threshold=kwargs.get('confidence_threshold', 0.5),
            log_level=kwargs.get('log_level', 'INFO')
        )
        
        # Process sample
        candidates = detector.detect_chimeras(
            assembly_file=assembly_file,
            reads1=reads1,
            reads2=reads2,
            temp_dir=kwargs.get('temp_dir')
        )
        
        if candidates:
            # Create a temporary BAM file path for analysis
            # In real implementation, this would come from the detector
            with tempfile.NamedTemporaryFile(suffix='.bam', delete=False) as tmp_bam:
                temp_bam_path = tmp_bam.name
            
            analyses = analyzer.analyze_chimeras(
                candidates=candidates,
                assembly_file=assembly_file,
                bam_file=temp_bam_path
            )
            
            output_files = resolver.resolve_chimeras(
                analyses=analyses,
                assembly_file=assembly_file,
                output_dir=sample_output_dir
            )
            
            # Generate sample report
            if kwargs.get('generate_report', True):
                visualizer = ChimeraVisualizer(log_level=kwargs.get('log_level', 'INFO'))
                report_path = visualizer.create_report(
                    analyses=analyses,
                    decisions=resolver.splitting_decisions,
                    output_dir=sample_output_dir,
                    assembly_file=assembly_file
                )
        
        return sample_output_dir
    
    def _process_samples_merged(self,
                              assembly_file: str,
                              sample_files: Dict[str, Tuple[str, Optional[str]]],
                              output_dir: str,
                              **kwargs) -> Dict[str, str]:
        """Process all samples as one merged dataset."""
        
        self.logger.info("Merging all samples for combined analysis")
        
        # Merge all read files
        merged_reads1, merged_reads2 = self._merge_read_files(sample_files, output_dir)
        
        # Process as single sample
        result = self._process_single_sample(
            assembly_file, merged_reads1, merged_reads2, "merged_all",
            output_dir, **kwargs
        )
        
        return {"merged_all": result}
    
    def _process_samples_batch(self,
                             assembly_file: str,
                             sample_files: Dict[str, Tuple[str, Optional[str]]],
                             output_dir: str,
                             batch_size: int = 5,
                             **kwargs) -> Dict[str, str]:
        """Process samples in batches for memory efficiency."""
        
        self.logger.info(f"Processing {len(sample_files)} samples in batches of {batch_size}")
        
        results = {}
        sample_names = list(sample_files.keys())
        
        for i in range(0, len(sample_names), batch_size):
            batch_samples = {name: sample_files[name] 
                           for name in sample_names[i:i+batch_size]}
            
            self.logger.info(f"Processing batch {i//batch_size + 1} "
                           f"({len(batch_samples)} samples)")
            
            batch_results = self._process_samples_separately(
                assembly_file, batch_samples, output_dir, **kwargs
            )
            results.update(batch_results)
        
        return results
    
    def _merge_read_files(self, 
                         sample_files: Dict[str, Tuple[str, Optional[str]]],
                         output_dir: str) -> Tuple[str, Optional[str]]:
        """Merge read files from multiple samples."""
        
        output_dir = Path(output_dir)
        merged_reads1 = output_dir / "merged_reads_R1.fastq.gz"
        merged_reads2 = output_dir / "merged_reads_R2.fastq.gz" if any(r2 for _, r2 in sample_files.values()) else None
        
        # Merge R1 files
        with open(merged_reads1, 'wb') as outfile:
            for sample_name, (reads1, reads2) in sample_files.items():
                with open(reads1, 'rb') as infile:
                    outfile.write(infile.read())
        
        # Merge R2 files if they exist
        if merged_reads2:
            with open(merged_reads2, 'wb') as outfile:
                for sample_name, (reads1, reads2) in sample_files.items():
                    if reads2:
                        with open(reads2, 'rb') as infile:
                            outfile.write(infile.read())
        
        return str(merged_reads1), str(merged_reads2) if merged_reads2 else None
    
    def _generate_multi_sample_report(self, sample_outputs: List[str], output_dir: str):
        """Generate combined report across all samples."""
        
        self.logger.info("Generating multi-sample summary report")
        
        # Collect results from all samples
        all_results = []
        sample_summaries = []
        
        for sample_output in sample_outputs:
            sample_name = Path(sample_output).name
            results_file = Path(sample_output) / "chimeric_detective_results.json"
            
            if results_file.exists():
                with open(results_file, 'r') as f:
                    sample_data = json.load(f)
                    sample_data['sample_name'] = sample_name
                    all_results.append(sample_data)
                    
                    # Extract summary stats
                    summary = sample_data.get('summary', {})
                    summary['sample_name'] = sample_name
                    sample_summaries.append(summary)
        
        # Create summary DataFrame
        if sample_summaries:
            summary_df = pd.DataFrame(sample_summaries)
            summary_path = Path(output_dir) / "multi_sample_summary.tsv"
            summary_df.to_csv(summary_path, sep='\t', index=False)
            
            # Generate summary statistics
            total_analyses = summary_df['total_analyses'].sum()
            total_split = summary_df['contigs_split'].sum()
            total_preserved = summary_df['contigs_preserved'].sum()
            
            self.logger.info(f"Multi-sample summary:")
            self.logger.info(f"  - Total samples: {len(sample_summaries)}")
            self.logger.info(f"  - Total chimera analyses: {total_analyses}")
            self.logger.info(f"  - Total contigs split: {total_split}")
            self.logger.info(f"  - Total contigs preserved: {total_preserved}")
        
        # Save combined results
        combined_results = {
            'processing_mode': self.processing_mode,
            'total_samples': len(sample_outputs),
            'sample_results': all_results,
            'summary_statistics': {
                'total_samples': len(sample_outputs),
                'total_analyses': sum(r.get('summary', {}).get('total_analyses', 0) for r in all_results),
                'total_split': sum(r.get('summary', {}).get('contigs_split', 0) for r in all_results),
                'total_preserved': sum(r.get('summary', {}).get('contigs_preserved', 0) for r in all_results)
            }
        }
        
        combined_results_path = Path(output_dir) / "combined_results.json"
        with open(combined_results_path, 'w') as f:
            json.dump(combined_results, f, indent=2)
        
        # Generate HTML summary report
        self._generate_multi_sample_html_report(combined_results, output_dir)
    
    def _generate_multi_sample_html_report(self, combined_results: Dict, output_dir: str):
        """Generate HTML report summarizing all samples."""
        
        from jinja2 import Template
        
        html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>Multi-Sample Chimeric Detective Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { text-align: center; color: #333; }
        .summary { background: #f5f5f5; padding: 20px; border-radius: 10px; margin: 20px 0; }
        .sample-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }
        .sample-card { background: white; border: 1px solid #ddd; padding: 15px; border-radius: 8px; }
        .sample-card h3 { margin-top: 0; color: #4CAF50; }
        .stats { display: flex; justify-content: space-around; text-align: center; }
        .stat { margin: 10px; }
        .stat-number { font-size: 2em; font-weight: bold; color: #4CAF50; }
        .stat-label { color: #666; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; }
        .chimera-high { color: #dc3545; font-weight: bold; }
        .chimera-medium { color: #ffc107; font-weight: bold; }
        .chimera-low { color: #28a745; font-weight: bold; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 Multi-Sample Chimeric Detective Report</h1>
        <p>Analysis Summary Across {{ combined_results.total_samples }} Samples</p>
    </div>
    
    <div class="summary">
        <h2>📊 Overall Summary</h2>
        <div class="stats">
            <div class="stat">
                <div class="stat-number">{{ combined_results.total_samples }}</div>
                <div class="stat-label">Total Samples</div>
            </div>
            <div class="stat">
                <div class="stat-number">{{ combined_results.summary_statistics.total_analyses }}</div>
                <div class="stat-label">Total Chimeras</div>
            </div>
            <div class="stat">
                <div class="stat-number">{{ combined_results.summary_statistics.total_split }}</div>
                <div class="stat-label">Contigs Split</div>
            </div>
            <div class="stat">
                <div class="stat-number">{{ combined_results.summary_statistics.total_preserved }}</div>
                <div class="stat-label">Contigs Preserved</div>
            </div>
        </div>
    </div>
    
    <h2>📋 Per-Sample Results</h2>
    <table>
        <thead>
            <tr>
                <th>Sample Name</th>
                <th>Total Analyses</th>
                <th>Contigs Split</th>
                <th>Contigs Preserved</th>
                <th>Contigs Flagged</th>
                <th>Report Link</th>
            </tr>
        </thead>
        <tbody>
            {% for result in combined_results.sample_results %}
            <tr>
                <td>{{ result.sample_name }}</td>
                <td>{{ result.summary.total_analyses }}</td>
                <td class="chimera-high">{{ result.summary.contigs_split }}</td>
                <td class="chimera-low">{{ result.summary.contigs_preserved }}</td>
                <td class="chimera-medium">{{ result.summary.contigs_flagged }}</td>
                <td><a href="{{ result.sample_name }}/chimeric_detective_report.html">View Report</a></td>
            </tr>
            {% endfor %}
        </tbody>
    </table>
    
    <div class="sample-grid">
        {% for result in combined_results.sample_results %}
        <div class="sample-card">
            <h3>{{ result.sample_name }}</h3>
            <p><strong>Processing Mode:</strong> {{ combined_results.processing_mode }}</p>
            <p><strong>Chimeras Detected:</strong> {{ result.summary.total_analyses }}</p>
            <p><strong>High Confidence:</strong> {{ result.summary.high_confidence }}</p>
            <p><strong>Medium Confidence:</strong> {{ result.summary.medium_confidence }}</p>
            <p><strong>Low Confidence:</strong> {{ result.summary.low_confidence }}</p>
            <p><a href="{{ result.sample_name }}/chimeric_detective_report.html">📊 View Detailed Report</a></p>
        </div>
        {% endfor %}
    </div>
    
    <footer style="text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666;">
        <p>Generated by Chimeric Detective Multi-Sample Processor</p>
    </footer>
</body>
</html>
        """
        
        template = Template(html_template)
        html_content = template.render(combined_results=combined_results)
        
        report_path = Path(output_dir) / "multi_sample_report.html"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"Multi-sample HTML report generated: {report_path}")


def process_multi_sample_directory(assembly_file: str,
                                  reads_dir: str,
                                  reads_pattern: str = "*_R{1,2}.fastq.gz",
                                  output_dir: str = "multi_sample_results",
                                  processing_mode: str = "separate",
                                  max_workers: int = 4,
                                  **kwargs) -> Dict[str, str]:
    """
    Convenience function to process multiple samples from a directory.
    
    Args:
        assembly_file: Path to assembly FASTA file
        reads_dir: Directory containing read files
        reads_pattern: Pattern for finding read files
        output_dir: Output directory for results
        processing_mode: How to process samples ("separate", "merged", "batch")
        max_workers: Maximum parallel workers
        **kwargs: Additional parameters for analysis
        
    Returns:
        Dictionary mapping sample names to output directories
    """
    
    processor = MultiSampleProcessor(
        processing_mode=processing_mode,
        max_workers=max_workers,
        log_level=kwargs.get('log_level', 'INFO')
    )
    
    return processor.process_samples_directory(
        assembly_file=assembly_file,
        reads_dir=reads_dir,
        reads_pattern=reads_pattern,
        output_dir=output_dir,
        **kwargs
    )