"""
Chimera detection module for identifying potential chimeric contigs.
"""

import logging
import os
import tempfile
from typing import List, Dict, Tuple, Optional, Set, Union
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
import pysam
from scipy import stats
from scipy.signal import find_peaks
from sklearn.cluster import DBSCAN

from .utils import (
    calculate_gc_content, calculate_kmer_frequencies, calculate_kmer_distance,
    run_command, merge_overlapping_intervals, setup_logging
)


@dataclass
class ChimeraCandidate:
    """Data class for storing chimera candidate information."""
    contig_id: str
    breakpoint: int
    confidence_score: float
    evidence_types: List[str]
    coverage_left: float
    coverage_right: float
    gc_content_left: float
    gc_content_right: float
    kmer_distance: float
    spanning_reads: int
    read_orientation_score: float


class ChimeraDetector:
    """Main class for detecting chimeric contigs in viral assemblies."""
    
    def __init__(self, 
                 min_contig_length: int = 1000,
                 min_coverage: float = 5.0,
                 coverage_fold_change: float = 2.0,
                 gc_content_threshold: float = 0.1,
                 kmer_distance_threshold: float = 0.3,
                 window_size: int = 1000,
                 step_size: int = 500,
                 min_spanning_reads: int = 5,
                 log_level: str = "INFO"):
        """
        Initialize ChimeraDetector.
        
        Args:
            min_contig_length: Minimum contig length to analyze
            min_coverage: Minimum coverage threshold
            coverage_fold_change: Minimum fold change in coverage to consider
            gc_content_threshold: Minimum GC content difference to consider
            kmer_distance_threshold: Minimum k-mer distance to consider
            window_size: Window size for sliding window analysis
            step_size: Step size for sliding window
            min_spanning_reads: Minimum spanning reads required
            log_level: Logging level
        """
        self.min_contig_length = min_contig_length
        self.min_coverage = min_coverage
        self.coverage_fold_change = coverage_fold_change
        self.gc_content_threshold = gc_content_threshold
        self.kmer_distance_threshold = kmer_distance_threshold
        self.window_size = window_size
        self.step_size = step_size
        self.min_spanning_reads = min_spanning_reads
        
        self.logger = setup_logging(log_level)
        self.chimera_candidates: List[ChimeraCandidate] = []
        
    def detect_chimeras(self, 
                       assembly_file: str,
                       bam_file: Optional[str] = None,
                       reads1: Optional[str] = None,
                       reads2: Optional[str] = None,
                       temp_dir: Optional[str] = None,
                       return_bam_path: bool = False) -> Union[List[ChimeraCandidate], Tuple[List[ChimeraCandidate], str]]:
        """
        Main method to detect chimeric contigs.
        
        Args:
            assembly_file: Path to assembly FASTA file
            bam_file: Path to BAM file with aligned reads
            reads1: Path to forward reads (if BAM not provided)
            reads2: Path to reverse reads (if BAM not provided)
            temp_dir: Temporary directory for intermediate files
            return_bam_path: If True, return tuple of (candidates, bam_path)
            
        Returns:
            List of ChimeraCandidate objects, or tuple of (candidates, bam_path) if return_bam_path=True
        """
        self.logger.info("Starting chimera detection")
        
        # Load assembly
        contigs = self._load_assembly(assembly_file)
        self.logger.info(f"Loaded {len(contigs)} contigs from assembly")
        
        # Prepare BAM file
        bam_file_to_use = bam_file
        if bam_file_to_use is None:
            if reads1 is None:
                raise ValueError("Either bam_file or reads1 must be provided")
            bam_file_to_use = self._align_reads(assembly_file, reads1, reads2, temp_dir)
        
        # Analyze each contig
        candidates = []
        for contig_id, sequence in contigs.items():
            if len(sequence) < self.min_contig_length:
                continue
                
            self.logger.debug(f"Analyzing contig {contig_id}")
            contig_candidates = self._analyze_contig(contig_id, sequence, bam_file_to_use)
            candidates.extend(contig_candidates)
        
        self.chimera_candidates = candidates
        self.logger.info(f"Detected {len(candidates)} chimera candidates")
        
        if return_bam_path:
            return candidates, bam_file_to_use
        else:
            return candidates
    
    def _load_assembly(self, assembly_file: str) -> Dict[str, str]:
        """Load assembly sequences from FASTA file."""
        contigs = {}
        for record in SeqIO.parse(assembly_file, "fasta"):
            contigs[record.id] = str(record.seq).upper()
        return contigs
    
    def _align_reads(self, 
                    assembly_file: str,
                    reads1: str,
                    reads2: Optional[str] = None,
                    temp_dir: Optional[str] = None) -> str:
        """Align reads to assembly and return BAM file path."""
        if temp_dir is None:
            temp_dir = tempfile.mkdtemp()
        
        temp_dir = Path(temp_dir)
        sam_file = temp_dir / "aligned.sam"
        bam_file = temp_dir / "aligned.bam"
        sorted_bam = temp_dir / "aligned_sorted.bam"
        
        self.logger.info("Aligning reads to assembly")
        
        # Check assembly size for memory considerations
        assembly_size = os.path.getsize(assembly_file) / (1024 * 1024)  # Size in MB
        if assembly_size > 500:  # Large assembly (>500MB)
            self.logger.warning(f"Large assembly detected ({assembly_size:.1f}MB). This may require significant memory for alignment.")
        
        # Check available tools first
        from .utils import check_external_tools
        available_tools = check_external_tools()
        
        aligner_success = False
        
        # Try minimap2 first if available
        if available_tools.get('minimap2', False):
            try:
                if reads2:
                    # Paired-end with memory-efficient options for large assemblies
                    cmd = ["minimap2", "-ax", "sr", "-t", "4", "--split-prefix", str(temp_dir / "tmp"),
                           "-I", "8G", assembly_file, reads1, reads2]
                else:
                    # Single-end with memory-efficient options
                    cmd = ["minimap2", "-ax", "sr", "-t", "4", "--split-prefix", str(temp_dir / "tmp"),
                           "-I", "8G", assembly_file, reads1]
                
                with open(sam_file, 'w') as f:
                    result = run_command(cmd, capture_output=True)
                    f.write(result.stdout)
                aligner_success = True
                self.logger.info("Successfully aligned reads using minimap2")
                
            except Exception as e:
                self.logger.warning(f"minimap2 failed: {e}")
                # If minimap2 fails due to memory, suggest BWA as it's more memory efficient
                if "SIGKILL" in str(e) or "killed" in str(e).lower():
                    self.logger.warning("minimap2 was killed (likely due to memory), trying BWA which is more memory efficient")
        else:
            self.logger.warning("minimap2 not available, trying BWA")
        
        # Try BWA if minimap2 failed or not available
        if not aligner_success and available_tools.get('bwa', False):
            try:
                # Index assembly
                run_command(["bwa", "index", assembly_file])
                
                # Align reads
                if reads2:
                    cmd = ["bwa", "mem", assembly_file, reads1, reads2]
                else:
                    cmd = ["bwa", "mem", assembly_file, reads1]
                
                with open(sam_file, 'w') as f:
                    result = run_command(cmd, capture_output=True)
                    f.write(result.stdout)
                aligner_success = True
                self.logger.info("Successfully aligned reads using BWA")
                
            except Exception as bwa_error:
                self.logger.error(f"BWA failed: {bwa_error}")
        
        if not aligner_success:
            error_msg = "No alignment tool succeeded. "
            if not available_tools.get('minimap2', False) and not available_tools.get('bwa', False):
                error_msg += "Neither minimap2 nor BWA are available. Please install one of these aligners."
            else:
                error_msg += f"Available tools: {[k for k, v in available_tools.items() if v]}. "
                if assembly_size > 500:
                    error_msg += f"Large assembly ({assembly_size:.1f}MB) may be causing memory issues. "
                    error_msg += "Consider: 1) Running on a machine with more RAM, 2) Filtering assembly to remove small contigs, "
                    error_msg += "3) Using --max-workers 1 to reduce parallel processing, or 4) Splitting analysis by contig."
                else:
                    error_msg += "Check tool installation and file paths."
            raise RuntimeError(error_msg)
        
        # Check if samtools is available for BAM processing
        if not available_tools.get('samtools', False):
            raise RuntimeError("samtools is required but not available. Please install samtools.")
        
        # Convert to BAM and sort
        run_command(["samtools", "view", "-bS", str(sam_file), "-o", str(bam_file)])
        run_command(["samtools", "sort", str(bam_file), "-o", str(sorted_bam)])
        run_command(["samtools", "index", str(sorted_bam)])
        
        return str(sorted_bam)
    
    def _analyze_contig(self, contig_id: str, sequence: str, bam_file: str) -> List[ChimeraCandidate]:
        """Analyze a single contig for chimeric signatures with adaptive window sizing."""
        candidates = []
        
        # Critical Fix #2: Adaptive window sizing based on contig length
        contig_length = len(sequence)
        adaptive_window = max(200, min(2000, contig_length // 20))
        adaptive_step = max(50, adaptive_window // 8)  # Finer step size for better resolution
        
        self.logger.debug(f"Contig {contig_id}: length={contig_length}, window={adaptive_window}, step={adaptive_step}")
        
        # Calculate coverage profile
        coverage = self._calculate_coverage(contig_id, len(sequence), bam_file)
        
        # Calculate GC content profile with adaptive window
        gc_profile = self._calculate_adaptive_gc_profile(sequence, adaptive_window, adaptive_step)
        
        # Calculate k-mer profiles with adaptive window
        kmer_profile = self._calculate_adaptive_kmer_profile(sequence, adaptive_window, adaptive_step)
        
        # Detect breakpoints using multiple methods with adaptive parameters
        coverage_breakpoints = self._detect_coverage_breakpoints_adaptive(coverage, adaptive_window)
        gc_breakpoints = self._detect_gc_breakpoints_adaptive(gc_profile, sequence, adaptive_step)
        kmer_breakpoints = self._detect_kmer_breakpoints_adaptive(kmer_profile, sequence, adaptive_step)
        
        # Critical Fix #4: Sub-grid detection - scan between detected breakpoints
        sub_grid_breakpoints = self._detect_sub_grid_breakpoints(
            sequence, coverage_breakpoints + gc_breakpoints + kmer_breakpoints, adaptive_window
        )
        
        # Combine all breakpoints
        all_breakpoints = set(coverage_breakpoints + gc_breakpoints + kmer_breakpoints + sub_grid_breakpoints)
        
        for initial_breakpoint in all_breakpoints:
            # Critical Fix #1: True breakpoint refinement at nucleotide resolution
            refined_breakpoint = self._refine_breakpoint(sequence, initial_breakpoint, adaptive_window // 4)
            
            candidate = self._evaluate_breakpoint_adaptive(
                contig_id, sequence, refined_breakpoint, bam_file,
                coverage, gc_profile, kmer_profile, adaptive_window
            )
            if candidate:
                candidates.append(candidate)
        
        return candidates
    
    def _calculate_coverage(self, contig_id: str, contig_length: int, bam_file: str) -> np.ndarray:
        """Calculate coverage profile for a contig."""
        coverage = np.zeros(contig_length)
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch(contig_id):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                
                start = read.reference_start
                end = read.reference_end
                
                if start is not None and end is not None:
                    coverage[start:end] += 1
        
        return coverage
    
    def _calculate_kmer_profile(self, sequence: str) -> List[Dict[str, int]]:
        """Calculate k-mer frequency profiles in sliding windows."""
        kmer_profile = []
        seq_len = len(sequence)
        
        for i in range(0, seq_len - self.window_size + 1, self.step_size):
            window_seq = sequence[i:i + self.window_size]
            kmers = calculate_kmer_frequencies(window_seq, k=4)
            kmer_profile.append(kmers)
        
        return kmer_profile
    
    def _calculate_adaptive_gc_profile(self, sequence: str, window_size: int, step_size: int) -> List[float]:
        """Calculate GC content profile with adaptive window sizing."""
        gc_profile = []
        seq_len = len(sequence)
        
        for i in range(0, seq_len - window_size + 1, step_size):
            window_seq = sequence[i:i + window_size]
            gc_content = (window_seq.count('G') + window_seq.count('C')) / len(window_seq)
            gc_profile.append(gc_content)
        
        return gc_profile
    
    def _calculate_adaptive_kmer_profile(self, sequence: str, window_size: int, step_size: int) -> List[Dict[str, int]]:
        """Calculate k-mer frequency profiles with adaptive window sizing."""
        kmer_profile = []
        seq_len = len(sequence)
        
        for i in range(0, seq_len - window_size + 1, step_size):
            window_seq = sequence[i:i + window_size]
            kmers = calculate_kmer_frequencies(window_seq, k=4)
            kmer_profile.append(kmers)
        
        return kmer_profile
    
    def _refine_breakpoint(self, sequence: str, initial_pos: int, window: int = 100) -> int:
        """Critical Fix #1: Actually refine breakpoint position at nucleotide resolution."""
        # Ensure we have sufficient sequence on both sides for analysis
        min_margin = 50
        if initial_pos < min_margin or initial_pos > len(sequence) - min_margin:
            return initial_pos
        
        # Limit refinement window to avoid out-of-bounds issues
        safe_window = min(window, initial_pos - min_margin, len(sequence) - initial_pos - min_margin)
        if safe_window < 10:  # Too small to be meaningful
            return initial_pos
        
        # Analyze GC/k-mer patterns at 1bp resolution around initial_pos
        max_signal = 0.0
        best_position = initial_pos
        
        # Scan in fine resolution around the initial position
        for pos in range(initial_pos - safe_window, initial_pos + safe_window + 1):
            if pos < min_margin or pos > len(sequence) - min_margin:
                continue
            
            # Calculate GC content discontinuity at this position with safe bounds
            analysis_size = min(50, pos, len(sequence) - pos)
            if analysis_size < 20:  # Need minimum sequence for meaningful analysis
                continue
                
            left_seq = sequence[pos-analysis_size:pos]
            right_seq = sequence[pos:pos+analysis_size]
            
            if len(left_seq) >= 20 and len(right_seq) >= 20:
                left_gc = (left_seq.count('G') + left_seq.count('C')) / len(left_seq)
                right_gc = (right_seq.count('G') + right_seq.count('C')) / len(right_seq)
                gc_signal = abs(left_gc - right_gc)
                
                # Calculate k-mer composition discontinuity
                left_kmers = calculate_kmer_frequencies(left_seq, k=4)
                right_kmers = calculate_kmer_frequencies(right_seq, k=4)
                kmer_signal = calculate_kmer_distance(left_kmers, right_kmers)
                
                # Combined signal strength
                total_signal = gc_signal + kmer_signal
                
                if total_signal > max_signal:
                    max_signal = total_signal
                    best_position = pos
        
        return best_position
    
    def _detect_sub_grid_breakpoints(self, sequence: str, detected_breakpoints: List[int], window_size: int) -> List[int]:
        """Critical Fix #4: Add fine-scale analysis between grid points to catch breakpoints at arbitrary positions."""
        sub_grid_breakpoints = []
        
        if len(detected_breakpoints) < 2:
            return sub_grid_breakpoints
        
        # Sort breakpoints
        sorted_breakpoints = sorted(detected_breakpoints)
        
        # Scan between each pair of detected breakpoints
        for i in range(len(sorted_breakpoints) - 1):
            start_pos = sorted_breakpoints[i]
            end_pos = sorted_breakpoints[i + 1]
            
            # If gap is large enough, scan for additional breakpoints
            if end_pos - start_pos > window_size * 2:
                scan_start = start_pos + window_size // 2
                scan_end = end_pos - window_size // 2
                scan_step = max(25, window_size // 20)  # Fine scanning resolution
                
                for pos in range(scan_start, scan_end, scan_step):
                    # Quick signal check at this position
                    if self._quick_signal_check(sequence, pos):
                        # Refine this position
                        refined_pos = self._refine_breakpoint(sequence, pos, window_size // 8)
                        
                        # Only add if it's not too close to existing breakpoints
                        if all(abs(refined_pos - bp) > window_size // 4 for bp in detected_breakpoints + sub_grid_breakpoints):
                            sub_grid_breakpoints.append(refined_pos)
        
        return sub_grid_breakpoints
    
    def _quick_signal_check(self, sequence: str, pos: int, window: int = 100) -> bool:
        """Quick check for potential breakpoint signal at a position."""
        min_margin = 25  # Minimum sequence needed for meaningful GC analysis
        if pos < min_margin or pos > len(sequence) - min_margin:
            return False
        
        # Use safe window size
        safe_window = min(window, pos, len(sequence) - pos)
        if safe_window < min_margin:
            return False
        
        left_seq = sequence[pos-safe_window:pos]
        right_seq = sequence[pos:pos+safe_window]
        
        # Quick GC content check
        if len(left_seq) < min_margin or len(right_seq) < min_margin:
            return False
            
        left_gc = (left_seq.count('G') + left_seq.count('C')) / len(left_seq)
        right_gc = (right_seq.count('G') + right_seq.count('C')) / len(right_seq)
        
        # Return True if GC difference exceeds threshold
        return abs(left_gc - right_gc) > self.gc_content_threshold * 0.5
    
    def _detect_coverage_breakpoints_adaptive(self, coverage: np.ndarray, window_size: int) -> List[int]:
        """Detect coverage breakpoints with adaptive window sizing."""
        return self._detect_coverage_breakpoints(coverage, window_size)
    
    def _detect_coverage_breakpoints(self, coverage: np.ndarray, window_size: int = None) -> List[int]:
        """Detect potential breakpoints based on coverage discontinuities."""
        breakpoints = []
        
        if window_size is None:
            window_size = self.window_size
        
        # Smooth coverage to reduce noise
        window = np.ones(window_size) / window_size
        smoothed_coverage = np.convolve(coverage, window, mode='same')
        
        # Find significant changes in coverage
        for i in range(self.window_size, len(smoothed_coverage) - self.window_size):
            left_cov = np.mean(smoothed_coverage[i-self.window_size:i])
            right_cov = np.mean(smoothed_coverage[i:i+self.window_size])
            
            if left_cov > self.min_coverage and right_cov > self.min_coverage:
                fold_change = max(left_cov, right_cov) / min(left_cov, right_cov)
                if fold_change >= self.coverage_fold_change:
                    breakpoints.append(i)
        
        # Merge nearby breakpoints
        if breakpoints:
            merged_breakpoints = merge_overlapping_intervals(
                [(bp - 50, bp + 50) for bp in breakpoints], min_gap=200
            )
            breakpoints = [int((start + end) / 2) for start, end in merged_breakpoints]
        
        return breakpoints
    
    def _detect_gc_breakpoints(self, gc_profile: List[float]) -> List[int]:
        """Detect potential breakpoints based on GC content changes using signal analysis."""
        breakpoints = []
        
        if len(gc_profile) < 4:
            return breakpoints
        
        # Calculate GC differences between adjacent windows
        gc_diffs = []
        for i in range(len(gc_profile) - 1):
            diff = abs(gc_profile[i+1] - gc_profile[i])
            gc_diffs.append(diff)
        
        # Find peaks in GC differences that exceed threshold
        step = self.window_size // 2  # Position increment between windows
        
        for i, diff in enumerate(gc_diffs):
            if diff >= self.gc_content_threshold:
                # Check if this is a local maximum (actual discontinuity)
                is_peak = True
                
                # Compare with neighbors (if they exist)
                if i > 0 and gc_diffs[i-1] >= diff:
                    is_peak = False
                if i < len(gc_diffs) - 1 and gc_diffs[i+1] >= diff:
                    is_peak = False
                
                if is_peak:
                    # Calculate position as the boundary between windows
                    position = (i + 1) * step
                    
                    # Refine breakpoint position by analyzing at higher resolution
                    refined_position = self._refine_gc_breakpoint(gc_profile, i, step, position)
                    breakpoints.append(refined_position)
        
        return breakpoints
    
    def _detect_kmer_breakpoints(self, kmer_profile: List[Dict[str, int]]) -> List[int]:
        """Detect potential breakpoints based on k-mer composition changes using signal analysis."""
        breakpoints = []
        
        if len(kmer_profile) < 2:
            return breakpoints
        
        # Calculate k-mer distances between adjacent windows
        kmer_dists = []
        for i in range(len(kmer_profile) - 1):
            dist = calculate_kmer_distance(kmer_profile[i], kmer_profile[i+1])
            kmer_dists.append(dist)
        
        # Find peaks in k-mer distances that exceed threshold
        step = self.step_size  # Position increment between windows
        
        for i, dist in enumerate(kmer_dists):
            if dist >= self.kmer_distance_threshold:
                # Check if this is a local maximum (actual discontinuity)
                is_peak = True
                
                # Compare with neighbors (if they exist) 
                if i > 0 and kmer_dists[i-1] >= dist:
                    is_peak = False
                if i < len(kmer_dists) - 1 and kmer_dists[i+1] >= dist:
                    is_peak = False
                
                if is_peak:
                    # Calculate position as the boundary between windows
                    position = (i + 1) * step
                    
                    # Refine breakpoint position by analyzing at higher resolution
                    refined_position = self._refine_kmer_breakpoint(kmer_profile, i, step, position)
                    breakpoints.append(refined_position)
        
        return breakpoints
    
    def _detect_gc_breakpoints_adaptive(self, gc_profile: List[float], sequence: str, step_size: int) -> List[int]:
        """Detect GC breakpoints with adaptive parameters."""
        breakpoints = []
        
        if len(gc_profile) < 4:
            return breakpoints
        
        # Calculate GC differences between adjacent windows
        gc_diffs = []
        for i in range(len(gc_profile) - 1):
            diff = abs(gc_profile[i+1] - gc_profile[i])
            gc_diffs.append(diff)
        
        # Find peaks in GC differences that exceed threshold
        for i, diff in enumerate(gc_diffs):
            if diff >= self.gc_content_threshold:
                # Check if this is a local maximum (actual discontinuity)
                is_peak = True
                
                # Compare with neighbors (if they exist)
                if i > 0 and gc_diffs[i-1] >= diff:
                    is_peak = False
                if i < len(gc_diffs) - 1 and gc_diffs[i+1] >= diff:
                    is_peak = False
                
                if is_peak:
                    # Calculate position as the boundary between windows
                    position = (i + 1) * step_size
                    
                    # Ensure position is valid
                    if 0 < position < len(sequence):
                        breakpoints.append(position)
        
        return breakpoints
    
    def _detect_kmer_breakpoints_adaptive(self, kmer_profile: List[Dict[str, int]], sequence: str, step_size: int) -> List[int]:
        """Detect k-mer breakpoints with adaptive parameters."""
        breakpoints = []
        
        if len(kmer_profile) < 2:
            return breakpoints
        
        # Calculate k-mer distances between adjacent windows
        kmer_dists = []
        for i in range(len(kmer_profile) - 1):
            dist = calculate_kmer_distance(kmer_profile[i], kmer_profile[i+1])
            kmer_dists.append(dist)
        
        # Find peaks in k-mer distances that exceed threshold
        for i, dist in enumerate(kmer_dists):
            if dist >= self.kmer_distance_threshold:
                # Check if this is a local maximum (actual discontinuity)
                is_peak = True
                
                # Compare with neighbors (if they exist) 
                if i > 0 and kmer_dists[i-1] >= dist:
                    is_peak = False
                if i < len(kmer_dists) - 1 and kmer_dists[i+1] >= dist:
                    is_peak = False
                
                if is_peak:
                    # Calculate position as the boundary between windows
                    position = (i + 1) * step_size
                    
                    # Ensure position is valid
                    if 0 < position < len(sequence):
                        breakpoints.append(position)
        
        return breakpoints
    
    def _refine_gc_breakpoint(self, gc_profile: List[float], peak_idx: int, step: int, initial_position: int) -> int:
        """Refine GC breakpoint position by analyzing at higher resolution."""
        # For now, return the window boundary position
        # Future enhancement: analyze GC content at single nucleotide resolution around this position
        return initial_position
    
    def _refine_kmer_breakpoint(self, kmer_profile: List[Dict[str, int]], peak_idx: int, step: int, initial_position: int) -> int:
        """Refine k-mer breakpoint position by analyzing at higher resolution."""
        # For now, return the window boundary position  
        # Future enhancement: analyze k-mer composition at higher resolution around this position
        return initial_position
    
    def _evaluate_breakpoint(self, 
                           contig_id: str,
                           sequence: str,
                           breakpoint: int,
                           bam_file: str,
                           coverage: np.ndarray,
                           gc_profile: List[float],
                           kmer_profile: List[Dict[str, int]]) -> Optional[ChimeraCandidate]:
        """Evaluate a potential breakpoint and create a candidate."""
        
        # Calculate evidence scores
        evidence_types = []
        
        # Coverage evidence
        left_cov = np.mean(coverage[max(0, breakpoint-self.window_size):breakpoint])
        right_cov = np.mean(coverage[breakpoint:min(len(coverage), breakpoint+self.window_size)])
        
        if left_cov > 0 and right_cov > 0:
            cov_fold_change = max(left_cov, right_cov) / min(left_cov, right_cov)
            if cov_fold_change >= self.coverage_fold_change:
                evidence_types.append("coverage_discontinuity")
        
        # GC content evidence
        gc_idx = breakpoint // (self.window_size // 2)
        if 0 < gc_idx < len(gc_profile) - 1:
            left_gc = gc_profile[max(0, gc_idx-1)]
            right_gc = gc_profile[min(len(gc_profile)-1, gc_idx+1)]
            gc_diff = abs(left_gc - right_gc)
            
            if gc_diff >= self.gc_content_threshold:
                evidence_types.append("gc_content_shift")
        
        # K-mer evidence
        kmer_idx = breakpoint // self.step_size
        if 0 < kmer_idx < len(kmer_profile) - 1:
            left_kmers = kmer_profile[max(0, kmer_idx-1)]
            right_kmers = kmer_profile[min(len(kmer_profile)-1, kmer_idx+1)]
            kmer_dist = calculate_kmer_distance(left_kmers, right_kmers)
            
            if kmer_dist >= self.kmer_distance_threshold:
                evidence_types.append("kmer_composition_change")
        
        # Spanning reads evidence
        spanning_reads = self._count_spanning_reads(contig_id, breakpoint, bam_file)
        
        # Read orientation evidence
        orientation_score = self._calculate_orientation_score(contig_id, breakpoint, bam_file)
        
        # Calculate confidence score
        confidence_score = self._calculate_confidence_score(
            evidence_types, cov_fold_change if 'coverage_discontinuity' in evidence_types else 1.0,
            gc_diff if 'gc_content_shift' in evidence_types else 0.0,
            kmer_dist if 'kmer_composition_change' in evidence_types else 0.0,
            spanning_reads, orientation_score
        )
        
        # Only create candidate if we have sufficient evidence
        if len(evidence_types) >= 2 or confidence_score > 0.7:
            return ChimeraCandidate(
                contig_id=contig_id,
                breakpoint=breakpoint,
                confidence_score=confidence_score,
                evidence_types=evidence_types,
                coverage_left=left_cov,
                coverage_right=right_cov,
                gc_content_left=left_gc if 'gc_content_shift' in evidence_types else 0.0,
                gc_content_right=right_gc if 'gc_content_shift' in evidence_types else 0.0,
                kmer_distance=kmer_dist if 'kmer_composition_change' in evidence_types else 0.0,
                spanning_reads=spanning_reads,
                read_orientation_score=orientation_score
            )
        
        return None
    
    def _evaluate_breakpoint_adaptive(self, 
                                    contig_id: str,
                                    sequence: str,
                                    breakpoint: int,
                                    bam_file: str,
                                    coverage: np.ndarray,
                                    gc_profile: List[float],
                                    kmer_profile: List[Dict[str, int]],
                                    window_size: int) -> Optional[ChimeraCandidate]:
        """Evaluate a potential breakpoint with adaptive window sizing."""
        
        # Calculate evidence scores with adaptive window
        evidence_types = []
        
        # Ensure breakpoint is within valid bounds for analysis
        if breakpoint < 50 or breakpoint > len(sequence) - 50:
            return None
        
        # Use smaller analysis window to avoid out-of-bounds errors  
        analysis_window = min(window_size, breakpoint, len(sequence) - breakpoint) // 2
        analysis_window = max(25, analysis_window)  # Minimum 25bp window
        
        # Coverage evidence with safe bounds
        left_start = max(0, breakpoint - analysis_window)
        right_end = min(len(coverage), breakpoint + analysis_window)
        
        left_cov = np.mean(coverage[left_start:breakpoint]) if left_start < breakpoint else 0.0
        right_cov = np.mean(coverage[breakpoint:right_end]) if breakpoint < right_end else 0.0
        
        cov_fold_change = 1.0
        if left_cov > 0 and right_cov > 0:
            cov_fold_change = max(left_cov, right_cov) / min(left_cov, right_cov)
            if cov_fold_change >= self.coverage_fold_change:
                evidence_types.append("coverage_discontinuity")
        
        # GC content evidence using high-resolution analysis with safe bounds
        seq_window = min(50, analysis_window)
        left_seq_start = max(0, breakpoint - seq_window)
        right_seq_end = min(len(sequence), breakpoint + seq_window)
        
        left_seq = sequence[left_seq_start:breakpoint]
        right_seq = sequence[breakpoint:right_seq_end]
        
        left_gc = 0.0
        right_gc = 0.0
        gc_diff = 0.0
        
        if len(left_seq) >= 20 and len(right_seq) >= 20:
            left_gc = (left_seq.count('G') + left_seq.count('C')) / len(left_seq)
            right_gc = (right_seq.count('G') + right_seq.count('C')) / len(right_seq)
            gc_diff = abs(left_gc - right_gc)
            
            if gc_diff >= self.gc_content_threshold:
                evidence_types.append("gc_content_shift")
        
        # K-mer evidence using high-resolution analysis
        kmer_dist = 0.0
        if len(left_seq) >= 20 and len(right_seq) >= 20:
            left_kmers = calculate_kmer_frequencies(left_seq, k=4)
            right_kmers = calculate_kmer_frequencies(right_seq, k=4)
            kmer_dist = calculate_kmer_distance(left_kmers, right_kmers)
            
            if kmer_dist >= self.kmer_distance_threshold:
                evidence_types.append("kmer_composition_change")
        
        # Spanning reads evidence
        spanning_reads = self._count_spanning_reads(contig_id, breakpoint, bam_file)
        
        # Read orientation evidence
        orientation_score = self._calculate_orientation_score(contig_id, breakpoint, bam_file)
        
        # Calculate confidence score
        confidence_score = self._calculate_confidence_score(
            evidence_types, cov_fold_change,
            gc_diff, kmer_dist, spanning_reads, orientation_score
        )
        
        # Only create candidate if we have sufficient evidence
        if len(evidence_types) >= 2 or confidence_score > 0.7:
            return ChimeraCandidate(
                contig_id=contig_id,
                breakpoint=breakpoint,
                confidence_score=confidence_score,
                evidence_types=evidence_types,
                coverage_left=left_cov,
                coverage_right=right_cov,
                gc_content_left=left_gc,
                gc_content_right=right_gc,
                kmer_distance=kmer_dist,
                spanning_reads=spanning_reads,
                read_orientation_score=orientation_score
            )
        
        return None
    
    def _count_spanning_reads(self, contig_id: str, breakpoint: int, bam_file: str) -> int:
        """Count reads that span the breakpoint."""
        spanning_count = 0
        window = 100  # Window around breakpoint
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch(contig_id, 
                                breakpoint - window, 
                                breakpoint + window):
                if (read.reference_start <= breakpoint - 50 and 
                    read.reference_end >= breakpoint + 50):
                    spanning_count += 1
        
        return spanning_count
    
    def _calculate_orientation_score(self, contig_id: str, breakpoint: int, bam_file: str) -> float:
        """Calculate read pair orientation score around breakpoint."""
        proper_pairs = 0
        improper_pairs = 0
        window = 500
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam.fetch(contig_id, 
                                breakpoint - window, 
                                breakpoint + window):
                if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
                    if read.is_proper_pair:
                        proper_pairs += 1
                    else:
                        improper_pairs += 1
        
        total_pairs = proper_pairs + improper_pairs
        if total_pairs == 0:
            return 0.0
        
        return improper_pairs / total_pairs
    
    def _calculate_confidence_score(self, 
                                  evidence_types: List[str],
                                  cov_fold_change: float,
                                  gc_diff: float,
                                  kmer_dist: float,
                                  spanning_reads: int,
                                  orientation_score: float) -> float:
        """Calculate overall confidence score for a chimera candidate."""
        
        # Base score from number of evidence types
        base_score = len(evidence_types) * 0.2
        
        # Coverage evidence weight
        if "coverage_discontinuity" in evidence_types:
            cov_score = min(1.0, (cov_fold_change - 2.0) / 8.0)  # Normalize to 0-1
            base_score += cov_score * 0.3
        
        # GC content evidence weight
        if "gc_content_shift" in evidence_types:
            gc_score = min(1.0, gc_diff / 0.3)  # Normalize to 0-1
            base_score += gc_score * 0.2
        
        # K-mer evidence weight
        if "kmer_composition_change" in evidence_types:
            kmer_score = min(1.0, kmer_dist / 0.5)  # Normalize to 0-1
            base_score += kmer_score * 0.2
        
        # Spanning reads penalty (fewer spanning reads = higher confidence for technical chimera)
        spanning_penalty = max(0, (spanning_reads - self.min_spanning_reads) / 20.0)
        base_score -= spanning_penalty * 0.1
        
        # Orientation score bonus
        base_score += orientation_score * 0.1
        
        return max(0.0, min(1.0, base_score))
    
    def get_summary_stats(self) -> Dict:
        """Get summary statistics of detected chimeras."""
        if not self.chimera_candidates:
            return {}
        
        confidence_scores = [c.confidence_score for c in self.chimera_candidates]
        evidence_counts = [len(c.evidence_types) for c in self.chimera_candidates]
        
        stats = {
            "total_candidates": len(self.chimera_candidates),
            "high_confidence": len([c for c in self.chimera_candidates if c.confidence_score > 0.8]),
            "medium_confidence": len([c for c in self.chimera_candidates if 0.5 < c.confidence_score <= 0.8]),
            "low_confidence": len([c for c in self.chimera_candidates if c.confidence_score <= 0.5]),
            "mean_confidence": np.mean(confidence_scores),
            "median_confidence": np.median(confidence_scores),
            "mean_evidence_types": np.mean(evidence_counts),
            "evidence_type_counts": {}
        }
        
        # Count evidence types
        all_evidence_types = []
        for candidate in self.chimera_candidates:
            all_evidence_types.extend(candidate.evidence_types)
        
        for evidence_type in set(all_evidence_types):
            stats["evidence_type_counts"][evidence_type] = all_evidence_types.count(evidence_type)
        
        return stats