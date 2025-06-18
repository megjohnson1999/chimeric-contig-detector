"""
Core read-pair based chimera detector.

Main class that orchestrates the entire detection workflow.
"""

import logging
from typing import List, Dict, Any, Optional
from pathlib import Path
import time

from .config import DetectorConfig
from .bamparser import BamParser, BamValidator
from .analyzer import ReadPairAnalyzer, BreakpointCandidate
from .output import OutputManager


class ReadPairChimeraDetector:
    """Main detector class for read-pair based chimera detection."""
    
    def __init__(self, config: DetectorConfig):
        """
        Initialize detector with configuration.
        
        Args:
            config: Detector configuration
        """
        self.config = config
        self._setup_logging()
        self.logger = logging.getLogger(__name__)
        
        # Initialize components
        self.analyzer = ReadPairAnalyzer(config)
        self.output_manager = OutputManager(config.output)
        
        # Runtime state
        self.results: List[BreakpointCandidate] = []
        self.metadata: Dict[str, Any] = {}
    
    def _setup_logging(self):
        """Configure logging based on config."""
        logger = logging.getLogger()
        logger.setLevel(getattr(logging, self.config.logging.level.upper()))
        
        # Clear existing handlers
        logger.handlers.clear()
        
        # Console handler
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter(self.config.logging.format)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler if specified
        if self.config.logging.file:
            file_handler = logging.FileHandler(self.config.logging.file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    
    def detect_chimeras(self, bam_path: str, assembly_path: Optional[str] = None,
                       contigs: Optional[List[str]] = None) -> List[BreakpointCandidate]:
        """
        Detect chimeric breakpoints in BAM file.
        
        Args:
            bam_path: Path to indexed BAM file
            assembly_path: Optional assembly FASTA file for validation
            contigs: Optional list of specific contigs to analyze
            
        Returns:
            List of breakpoint candidates
        """
        start_time = time.time()
        self.logger.info(f"Starting chimera detection on {bam_path}")
        
        # Validate inputs
        self._validate_inputs(bam_path, assembly_path)
        
        # Initialize BAM parser
        with BamParser(bam_path, self.config.read_quality) as parser:
            # Get contig list
            target_contigs = self._get_target_contigs(parser, contigs)
            self.logger.info(f"Analyzing {len(target_contigs)} contigs")
            
            # Analyze each contig
            all_results = []
            for i, contig in enumerate(target_contigs):
                self.logger.info(f"Processing contig {i+1}/{len(target_contigs)}: {contig}")
                
                try:
                    contig_results = self.analyzer.analyze_contig(parser, contig)
                    
                    # Set contig name for results
                    for result in contig_results:
                        result.contig = contig
                    
                    all_results.extend(contig_results)
                    
                    self.logger.info(f"Found {len(contig_results)} candidates in {contig}")
                    
                except Exception as e:
                    self.logger.error(f"Error analyzing contig {contig}: {e}")
                    if self.config.logging.level == "DEBUG":
                        import traceback
                        traceback.print_exc()
        
        # Filter and rank results
        self.results = self._post_process_results(all_results)
        
        # Store metadata
        end_time = time.time()
        self.metadata = {
            "input_bam": bam_path,
            "input_assembly": assembly_path,
            "analysis_time_seconds": round(end_time - start_time, 2),
            "contigs_analyzed": len(target_contigs),
            "total_candidates": len(self.results),
            "config": self.config.to_dict(),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")
        }
        
        self.logger.info(f"Detection complete. Found {len(self.results)} candidates in {end_time - start_time:.1f}s")
        return self.results
    
    def write_results(self, output_dir: str) -> Dict[str, str]:
        """
        Write detection results to specified output directory.
        
        Args:
            output_dir: Output directory path
            
        Returns:
            Dictionary mapping format to output file path
        """
        if not self.results:
            self.logger.warning("No results to write")
            return {}
        
        self.logger.info(f"Writing results to {output_dir}")
        
        # Write formatted outputs
        output_files = self.output_manager.write_results(
            self.results, self.metadata, output_dir
        )
        
        # Create summary report
        summary_path = self.output_manager.create_summary_report(
            self.results, self.metadata, output_dir
        )
        output_files["summary"] = summary_path
        
        return output_files
    
    def get_high_confidence_candidates(self, threshold: float = 0.8) -> List[BreakpointCandidate]:
        """Get candidates above confidence threshold."""
        return [r for r in self.results if r.confidence >= threshold]
    
    def get_candidates_by_contig(self) -> Dict[str, List[BreakpointCandidate]]:
        """Group candidates by contig."""
        by_contig = {}
        for candidate in self.results:
            if candidate.contig not in by_contig:
                by_contig[candidate.contig] = []
            by_contig[candidate.contig].append(candidate)
        return by_contig
    
    def _validate_inputs(self, bam_path: str, assembly_path: Optional[str]):
        """Validate input files."""
        # Check BAM file
        bam_issues = BamValidator.validate_bam_file(bam_path)
        if bam_issues:
            for issue in bam_issues:
                self.logger.warning(f"BAM validation: {issue}")
        
        # Check assembly file if provided
        if assembly_path:
            if not Path(assembly_path).exists():
                raise FileNotFoundError(f"Assembly file not found: {assembly_path}")
    
    def _get_target_contigs(self, parser: BamParser, 
                           contigs: Optional[List[str]]) -> List[str]:
        """Get list of contigs to analyze."""
        if contigs:
            # Validate specified contigs exist
            missing = [c for c in contigs if c not in parser.contigs]
            if missing:
                self.logger.warning(f"Contigs not found in BAM: {missing}")
            return [c for c in contigs if c in parser.contigs]
        else:
            # Analyze all contigs with sufficient length
            min_length = self.config.sliding_window.window_size * 2
            return [contig for contig, length in parser.contigs.items() 
                   if length >= min_length]
    
    def _post_process_results(self, results: List[BreakpointCandidate]) -> List[BreakpointCandidate]:
        """Post-process and filter results."""
        if not results:
            return []
        
        # Sort by confidence
        results.sort(key=lambda x: x.confidence, reverse=True)
        
        # Apply additional filtering if needed
        filtered_results = []
        for result in results:
            # Filter by minimum anomaly length
            anomaly_length = max(a.end - a.start for a in result.supporting_anomalies)
            if anomaly_length >= self.config.detection.min_anomaly_length:
                filtered_results.append(result)
        
        self.logger.info(f"Post-processing: {len(results)} -> {len(filtered_results)} candidates")
        return filtered_results


class DetectorFactory:
    """Factory for creating detector instances with different configurations."""
    
    @staticmethod
    def create_default() -> ReadPairChimeraDetector:
        """Create detector with default configuration."""
        config = DetectorConfig()
        return ReadPairChimeraDetector(config)
    
    @staticmethod
    def create_from_file(config_path: str) -> ReadPairChimeraDetector:
        """Create detector from configuration file."""
        config = DetectorConfig.from_file(config_path)
        
        # Validate configuration
        issues = config.validate()
        if issues:
            logger = logging.getLogger(__name__)
            for issue in issues:
                logger.warning(f"Config validation: {issue}")
        
        return ReadPairChimeraDetector(config)
    
    @staticmethod
    def create_sensitive() -> ReadPairChimeraDetector:
        """Create detector optimized for high sensitivity."""
        config = DetectorConfig()
        
        # Adjust thresholds for higher sensitivity
        config.detection.proper_pair_drop_threshold = 0.3
        config.detection.insert_size_shift_threshold = 1.5
        config.detection.discordant_pair_threshold = 0.2
        
        # Smaller windows for better resolution
        config.sliding_window.window_size = 500
        config.sliding_window.step_size = 250
        
        return ReadPairChimeraDetector(config)
    
    @staticmethod
    def create_specific() -> ReadPairChimeraDetector:
        """Create detector optimized for high specificity."""
        config = DetectorConfig()
        
        # Adjust thresholds for higher specificity
        config.detection.proper_pair_drop_threshold = 0.7
        config.detection.insert_size_shift_threshold = 3.0
        config.detection.discordant_pair_threshold = 0.5
        
        # Require more evidence
        config.sliding_window.min_pairs_per_window = 100
        config.insert_size.min_pairs_for_baseline = 2000
        
        return ReadPairChimeraDetector(config)