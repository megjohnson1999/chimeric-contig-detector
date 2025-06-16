"""
Multi-sample processing module for handling directories with multiple samples.
"""

import os
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import json
import numpy as np
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
                 processing_mode: str = "separate",  # "separate", "merged", "batch", "coassembly"
                 max_workers: int = 4,
                 log_level: str = "INFO"):
        """
        Initialize MultiSampleProcessor.
        
        Args:
            processing_mode: How to process multiple samples
                - "separate": Process each sample independently
                - "merged": Merge all reads and process as one sample
                - "batch": Process in batches for memory efficiency
                - "coassembly": Process a co-assembly with multiple samples' reads
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
        elif self.processing_mode == "coassembly":
            return self._process_coassembly(
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
                    
                    # Remove conflicting parameters from kwargs before passing
                    sample_kwargs = kwargs.copy()
                    sample_kwargs.pop('reads1', None)
                    sample_kwargs.pop('reads2', None)
                    sample_kwargs.pop('reads', None)
                    sample_kwargs.pop('bam', None)
                    sample_kwargs.pop('assembly', None)
                    sample_kwargs.pop('reads_dir', None)
                    sample_kwargs.pop('reads_pattern', None)
                    sample_kwargs.pop('out', None)
                    
                    future = executor.submit(
                        self._process_single_sample,
                        assembly_file, reads1, reads2, sample_name, 
                        str(sample_output_dir), **sample_kwargs
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
                        # Create minimal output directory for failed sample
                        failed_output_dir = Path(output_dir) / f"sample_{sample_name}"
                        create_output_directory(str(failed_output_dir))
                        
                        # Save error status
                        self._save_failed_sample_results(sample_name, str(e), str(failed_output_dir))
                        results[sample_name] = str(failed_output_dir)
                        sample_outputs.append(str(failed_output_dir))
        
        else:
            # Sequential processing
            self.logger.info(f"Processing {len(sample_files)} samples sequentially")
            
            for sample_name, (reads1, reads2) in tqdm(sample_files.items(), 
                                                    desc="Processing samples"):
                sample_output_dir = Path(output_dir) / f"sample_{sample_name}"
                
                try:
                    # Remove conflicting parameters from kwargs before passing
                    sample_kwargs = kwargs.copy()
                    sample_kwargs.pop('reads1', None)
                    sample_kwargs.pop('reads2', None)
                    sample_kwargs.pop('reads', None)
                    sample_kwargs.pop('bam', None)
                    sample_kwargs.pop('assembly', None)
                    sample_kwargs.pop('reads_dir', None)
                    sample_kwargs.pop('reads_pattern', None)
                    sample_kwargs.pop('out', None)
                    
                    sample_result = self._process_single_sample(
                        assembly_file, reads1, reads2, sample_name,
                        str(sample_output_dir), **sample_kwargs
                    )
                    results[sample_name] = sample_result
                    sample_outputs.append(sample_result)
                except Exception as e:
                    self.logger.error(f"Failed to process sample {sample_name}: {e}")
                    # Create minimal output directory for failed sample
                    failed_output_dir = Path(output_dir) / f"sample_{sample_name}"
                    create_output_directory(str(failed_output_dir))
                    
                    # Save error status
                    self._save_failed_sample_results(sample_name, str(e), str(failed_output_dir))
                    results[sample_name] = str(failed_output_dir)
                    sample_outputs.append(str(failed_output_dir))
        
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
        candidates, bam_file = detector.detect_chimeras(
            assembly_file=assembly_file,
            reads1=reads1,
            reads2=reads2,
            temp_dir=kwargs.get('temp_dir'),
            return_bam_path=True
        )
        
        if candidates:
            if not bam_file:
                raise ValueError(f"No BAM file generated for sample {sample_name}")
            
            analyses = analyzer.analyze_chimeras(
                candidates=candidates,
                assembly_file=assembly_file,
                bam_file=bam_file
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
            
            # Save sample results summary
            self._save_sample_results_summary(
                sample_name, analyses, resolver.splitting_decisions, sample_output_dir
            )
        else:
            # No chimeras detected - still save empty results
            self._save_sample_results_summary(
                sample_name, [], [], sample_output_dir
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
        
        # Remove conflicting parameters from kwargs before passing
        sample_kwargs = kwargs.copy()
        sample_kwargs.pop('reads1', None)
        sample_kwargs.pop('reads2', None)
        sample_kwargs.pop('reads', None)
        sample_kwargs.pop('bam', None)
        sample_kwargs.pop('assembly', None)
        sample_kwargs.pop('reads_dir', None)
        sample_kwargs.pop('reads_pattern', None)
        sample_kwargs.pop('out', None)
        
        # Process as single sample
        result = self._process_single_sample(
            assembly_file, merged_reads1, merged_reads2, "merged_all",
            output_dir, **sample_kwargs
        )
        
        return {"merged_all": result}
    
    def _process_samples_batch(self,
                             assembly_file: str,
                             sample_files: Dict[str, Tuple[str, Optional[str]]],
                             output_dir: str,
                             **kwargs) -> Dict[str, str]:
        """Process samples in batches for memory efficiency."""
        
        batch_size = kwargs.get('batch_size', 5)
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
    
    def _process_coassembly(self,
                           assembly_file: str,
                           sample_files: Dict[str, Tuple[str, Optional[str]]],
                           output_dir: str,
                           **kwargs) -> Dict[str, str]:
        """Process a co-assembly with multiple samples' reads.
        
        This mode is for when you have a co-assembly created from multiple samples,
        and you want to use all samples' reads together to evaluate chimeras.
        Unlike 'merged' mode, this keeps track of per-sample coverage information.
        """
        
        self.logger.info(f"Processing co-assembly with {len(sample_files)} samples")
        
        # Initialize components with filtered kwargs
        detector_kwargs = {k: v for k, v in kwargs.items() if k in [
            'min_contig_length', 'min_coverage', 'coverage_fold_change', 
            'gc_content_threshold', 'kmer_distance_threshold', 'window_size',
            'step_size', 'min_spanning_reads', 'log_level'
        ]}
        
        analyzer_kwargs = {k: v for k, v in kwargs.items() if k in [
            'reference', 'log_level'
        ]}
        
        resolver_kwargs = {k: v for k, v in kwargs.items() if k in [
            'split_technical', 'split_pcr', 'preserve_biological',
            'min_split_length', 'confidence_threshold', 'log_level'
        ]}
        
        detector = ChimeraDetector(**detector_kwargs)
        analyzer = ChimeraAnalyzer(**analyzer_kwargs)
        resolver = ChimeraResolver(**resolver_kwargs)
        visualizer = ChimeraVisualizer()
        
        # Step 1: Map all samples' reads to the assembly and collect coverage data
        self.logger.info("Mapping all samples to co-assembly")
        sample_bam_files = {}
        sample_coverages = {}
        
        # Step 1: Run chimera detection for each sample separately
        self.logger.info("Running chimera detection for each sample")
        sample_candidates = {}
        sample_bam_files = {}
        
        with tempfile.TemporaryDirectory() as temp_dir:
            for sample_name, (reads1, reads2) in tqdm(sample_files.items(), 
                                                     desc="Processing samples"):
                self.logger.debug(f"Processing {sample_name}")
                
                # Run detection for this sample
                candidates, bam_file = detector.detect_chimeras(
                    assembly_file=assembly_file,
                    reads1=reads1,
                    reads2=reads2,
                    temp_dir=temp_dir,
                    return_bam_path=True
                )
                
                sample_candidates[sample_name] = candidates
                sample_bam_files[sample_name] = bam_file
            
            # Step 2: Aggregate candidates across samples
            self.logger.info("Aggregating chimera candidates across samples")
            candidates = self._aggregate_chimera_candidates(sample_candidates)
            
            self.logger.info(f"Detected {len(candidates)} aggregated chimera candidates")
            
            # Step 3: Analyze chimeras with multi-sample context
            self.logger.info(f"Analyzing {len(candidates)} chimera candidates")
            analyses = []
            
            for candidate in tqdm(candidates, desc="Analyzing candidates"):
                # Analyze with awareness of multiple samples
                analysis = analyzer.analyze_chimera(candidate)
                analysis.multi_sample_support = self._calculate_multi_sample_support(
                    candidate, sample_candidates
                )
                analyses.append(analysis)
            
            # Step 4: Resolve chimeras
            self.logger.info(f"Resolving {len(analyses)} chimera analyses")
            contigs = detector._load_assembly(assembly_file)  # Load contigs for resolver
            decisions = resolver.resolve_chimeras(analyses, contigs)
            
            # Step 5: Generate outputs
            resolver.write_outputs(decisions, contigs, output_dir)
            
            # Generate visualization report
            if kwargs.get('generate_report', True):
                self.logger.info("Creating interactive HTML report")
                report_path = visualizer.create_html_report(
                    analyses, decisions, output_dir
                )
            
            # Generate multi-sample summary
            summary = self._generate_coassembly_summary(
                candidates, analyses, decisions, sample_files
            )
            
            summary_path = Path(output_dir) / "coassembly_summary.json"
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2)
            
            self.logger.info("Co-assembly analysis completed successfully!")
            
            return {"coassembly": str(output_dir)}
    
    def _aggregate_chimera_candidates(self, sample_candidates: Dict[str, List[ChimeraCandidate]]) -> List[ChimeraCandidate]:
        """Aggregate chimera candidates across samples."""
        # Create a mapping of (contig_id, breakpoint) -> list of candidates
        candidate_groups = {}
        
        for sample_name, candidates in sample_candidates.items():
            for candidate in candidates:
                # Group candidates by contig and approximate breakpoint (within 100bp)
                key = self._get_candidate_key(candidate)
                if key not in candidate_groups:
                    candidate_groups[key] = []
                candidate_groups[key].append((sample_name, candidate))
        
        # Create consensus candidates
        aggregated_candidates = []
        for group_key, sample_candidate_pairs in candidate_groups.items():
            consensus_candidate = self._create_consensus_candidate(sample_candidate_pairs)
            aggregated_candidates.append(consensus_candidate)
        
        return aggregated_candidates
    
    def _get_candidate_key(self, candidate: ChimeraCandidate) -> Tuple[str, int]:
        """Get a grouping key for candidate aggregation."""
        # Round breakpoint to nearest 100bp for grouping
        rounded_breakpoint = (candidate.breakpoint // 100) * 100
        return (candidate.contig_id, rounded_breakpoint)
    
    def _create_consensus_candidate(self, sample_candidate_pairs: List[Tuple[str, ChimeraCandidate]]) -> ChimeraCandidate:
        """Create a consensus candidate from multiple samples."""
        sample_names = [pair[0] for pair in sample_candidate_pairs]
        candidates = [pair[1] for pair in sample_candidate_pairs]
        
        # Use the first candidate as base and aggregate evidence
        base_candidate = candidates[0]
        
        # Calculate consensus values
        consensus_breakpoint = int(np.mean([c.breakpoint for c in candidates]))
        consensus_confidence = np.mean([c.confidence_score for c in candidates])
        
        # Aggregate evidence types
        all_evidence = set()
        for candidate in candidates:
            all_evidence.update(candidate.evidence_types)
        
        # Aggregate coverage (average across samples)
        avg_coverage_left = np.mean([c.coverage_left for c in candidates])
        avg_coverage_right = np.mean([c.coverage_right for c in candidates])
        
        # Create consensus candidate
        consensus = ChimeraCandidate(
            contig_id=base_candidate.contig_id,
            breakpoint=consensus_breakpoint,
            confidence_score=consensus_confidence,
            evidence_types=list(all_evidence),
            coverage_left=avg_coverage_left,
            coverage_right=avg_coverage_right,
            gc_content_left=np.mean([c.gc_content_left for c in candidates]),
            gc_content_right=np.mean([c.gc_content_right for c in candidates]),
            kmer_distance=np.mean([c.kmer_distance for c in candidates]),
            spanning_reads=int(np.mean([c.spanning_reads for c in candidates])),
            read_orientation_score=np.mean([c.read_orientation_score for c in candidates])
        )
        
        # Add multi-sample metadata
        consensus.supporting_samples = sample_names
        consensus.sample_count = len(sample_names)
        
        return consensus
    
    def _aggregate_sample_coverages_old(self, sample_coverages: Dict[str, Dict]) -> Dict:
        """Aggregate coverage information across samples."""
        aggregated = {}
        
        # Get all contigs
        all_contigs = set()
        for sample_cov in sample_coverages.values():
            all_contigs.update(sample_cov.keys())
        
        # Aggregate coverage for each contig
        for contig_id in all_contigs:
            contig_coverages = []
            
            for sample_name, coverages in sample_coverages.items():
                if contig_id in coverages:
                    contig_coverages.append(coverages[contig_id])
            
            if contig_coverages:
                # Calculate aggregated statistics
                aggregated[contig_id] = {
                    'mean_coverage': sum(c.get('mean_coverage', 0) for c in contig_coverages) / len(contig_coverages),
                    'coverage_array': self._merge_coverage_arrays([c.get('coverage', []) for c in contig_coverages]),
                    'num_samples': len(contig_coverages),
                    'sample_presence': len(contig_coverages) / len(sample_coverages)
                }
        
        return aggregated
    
    def _merge_coverage_arrays(self, coverage_arrays: List[List[float]]) -> List[float]:
        """Merge coverage arrays from multiple samples."""
        if not coverage_arrays:
            return []
        
        # Ensure all arrays have the same length
        max_len = max(len(arr) for arr in coverage_arrays)
        
        # Pad arrays to same length if needed
        padded_arrays = []
        for arr in coverage_arrays:
            if len(arr) < max_len:
                padded = arr + [0] * (max_len - len(arr))
                padded_arrays.append(padded)
            else:
                padded_arrays.append(arr)
        
        # Calculate mean coverage at each position
        merged = []
        for i in range(max_len):
            values = [arr[i] for arr in padded_arrays if i < len(arr)]
            merged.append(sum(values) / len(values) if values else 0)
        
        return merged
    
    def _calculate_multi_sample_support(self, candidate: ChimeraCandidate, 
                                      sample_candidates: Dict) -> Dict:
        """Calculate how many samples support this chimera candidate."""
        # If candidate has multi-sample metadata, use it
        if hasattr(candidate, 'supporting_samples'):
            support = {
                'supporting_samples': candidate.supporting_samples,
                'total_samples': len(sample_candidates),
                'sample_count': candidate.sample_count,
                'sample_support_ratio': candidate.sample_count / len(sample_candidates),
                'consistent_across_samples': candidate.sample_count > 1
            }
        else:
            # Fallback: calculate based on similar candidates in other samples
            contig_id = candidate.contig_id
            breakpoint = candidate.breakpoint
            supporting_samples = []
            
            for sample_name, candidates in sample_candidates.items():
                for sample_candidate in candidates:
                    if (sample_candidate.contig_id == contig_id and 
                        abs(sample_candidate.breakpoint - breakpoint) <= 100):
                        supporting_samples.append(sample_name)
                        break
            
            support = {
                'supporting_samples': supporting_samples,
                'total_samples': len(sample_candidates),
                'sample_count': len(supporting_samples),
                'sample_support_ratio': len(supporting_samples) / len(sample_candidates),
                'consistent_across_samples': len(supporting_samples) > 1
            }
        
        return support
    
    def _generate_coassembly_summary(self, candidates: List[ChimeraCandidate],
                                   analyses: List[ChimeraAnalysis],
                                   decisions: List[SplittingDecision],
                                   sample_files: Dict) -> Dict:
        """Generate summary for co-assembly analysis."""
        
        decision_dict = {d.contig_id: d for d in decisions}
        
        summary = {
            'mode': 'coassembly',
            'num_samples': len(sample_files),
            'sample_names': list(sample_files.keys()),
            'total_candidates': len(candidates),
            'total_analyses': len(analyses),
            'decisions': {
                'split': sum(1 for d in decisions if d.action == 'split'),
                'preserve': sum(1 for d in decisions if d.action == 'preserve'),
                'flag': sum(1 for d in decisions if d.action == 'flag')
            },
            'multi_sample_statistics': {
                'candidates_in_all_samples': 0,
                'candidates_in_multiple_samples': 0,
                'candidates_in_single_sample': 0
            }
        }
        
        # Calculate multi-sample statistics
        for analysis in analyses:
            if hasattr(analysis, 'multi_sample_support'):
                support = analysis.multi_sample_support
                num_samples = support.get('samples_with_contig', 0)
                
                if num_samples == len(sample_files):
                    summary['multi_sample_statistics']['candidates_in_all_samples'] += 1
                elif num_samples > 1:
                    summary['multi_sample_statistics']['candidates_in_multiple_samples'] += 1
                else:
                    summary['multi_sample_statistics']['candidates_in_single_sample'] += 1
        
        return summary
    
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
    
    def _save_sample_results_summary(self, 
                                   sample_name: str, 
                                   analyses: List, 
                                   decisions: List, 
                                   output_dir: str):
        """Save a summary of sample results to JSON file."""
        import json
        
        # Create summary data
        summary = {
            'sample_name': sample_name,
            'total_analyses': len(analyses),
            'contigs_split': len([d for d in decisions if d.action == 'split']),
            'contigs_preserved': len([d for d in decisions if d.action == 'preserve']),
            'contigs_flagged': len([d for d in decisions if d.action == 'flag']),
            'high_confidence': len([a for a in analyses if hasattr(a, 'confidence_score') and a.confidence_score > 0.8]),
            'medium_confidence': len([a for a in analyses if hasattr(a, 'confidence_score') and 0.5 < a.confidence_score <= 0.8]),
            'low_confidence': len([a for a in analyses if hasattr(a, 'confidence_score') and a.confidence_score <= 0.5]),
        }
        
        # Save to file
        results_data = {
            'sample_name': sample_name,
            'summary': summary,
            'timestamp': pd.Timestamp.now().isoformat(),
            'processing_status': 'completed'
        }
        
        results_file = Path(output_dir) / "chimeric_detective_results.json"
        with open(results_file, 'w') as f:
            json.dump(results_data, f, indent=2)
    
    def _save_failed_sample_results(self, sample_name: str, error_message: str, output_dir: str):
        """Save error information for failed sample."""
        import json
        
        # Create error summary
        summary = {
            'sample_name': sample_name,
            'total_analyses': 0,
            'contigs_split': 0,
            'contigs_preserved': 0,
            'contigs_flagged': 0,
            'high_confidence': 0,
            'medium_confidence': 0,
            'low_confidence': 0,
        }
        
        # Save error info
        results_data = {
            'sample_name': sample_name,
            'summary': summary,
            'timestamp': pd.Timestamp.now().isoformat(),
            'processing_status': 'failed',
            'error_message': error_message
        }
        
        results_file = Path(output_dir) / "chimeric_detective_results.json"
        with open(results_file, 'w') as f:
            json.dump(results_data, f, indent=2)


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
    
    # Remove conflicting parameters from kwargs
    clean_kwargs = kwargs.copy()
    clean_kwargs.pop('assembly_file', None)
    clean_kwargs.pop('reads_dir', None)
    clean_kwargs.pop('reads_pattern', None)
    clean_kwargs.pop('output_dir', None)
    clean_kwargs.pop('processing_mode', None)
    clean_kwargs.pop('max_workers', None)
    
    return processor.process_samples_directory(
        assembly_file=assembly_file,
        reads_dir=reads_dir,
        reads_pattern=reads_pattern,
        output_dir=output_dir,
        **clean_kwargs
    )