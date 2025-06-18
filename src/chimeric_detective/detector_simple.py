"""
Simplified chimera detection module using only GC content analysis.
"""

import logging
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import numpy as np
from Bio import SeqIO
from scipy import stats
from scipy.signal import find_peaks

from .utils import setup_logging
from .constants import (
    SequenceConstants, WindowConstants, DetectionThresholds,
    get_adaptive_window_size, get_adaptive_step_size, adjust_window_for_sequence_length
)


def calculate_gc_content_simple(sequence: str) -> float:
    """Calculate GC content of a sequence as a single value."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return gc_count / len(sequence)


@dataclass
class ChimeraCandidate:
    """Simplified data class for storing chimera candidate information."""
    contig_id: str
    breakpoint: int
    confidence_score: float
    gc_content_left: float
    gc_content_right: float
    gc_difference: float
    breakpoint_region: Optional[Tuple[int, int]] = None


class SimpleChimeraDetector:
    """Simplified class for detecting chimeric contigs using only GC content."""
    
    def __init__(self, 
                 min_contig_length: int = 1000,
                 gc_content_threshold: float = 0.1,
                 window_size: int = 1000,
                 step_size: int = 500,
                 log_level: str = "INFO"):
        """
        Initialize SimpleChimeraDetector.
        
        Args:
            min_contig_length: Minimum contig length to analyze
            gc_content_threshold: Minimum GC content difference to consider
            window_size: Window size for sliding window analysis
            step_size: Step size for sliding window
            log_level: Logging level
        """
        self.min_contig_length = min_contig_length
        self.gc_content_threshold = gc_content_threshold
        self.window_size = window_size
        self.step_size = step_size
        
        self.logger = setup_logging(log_level)
        self.chimera_candidates: List[ChimeraCandidate] = []
        
    def detect_chimeras(self, assembly_file: str) -> List[ChimeraCandidate]:
        """
        Detect chimeric contigs in an assembly file using GC content analysis.
        
        Args:
            assembly_file: Path to assembly FASTA file
            
        Returns:
            List of chimera candidates
        """
        self.logger.info(f"Starting simplified chimera detection on {assembly_file}")
        
        self.chimera_candidates = []
        
        # Process each contig
        for record in SeqIO.parse(assembly_file, "fasta"):
            contig_id = record.id
            sequence = str(record.seq).upper()
            
            # Analyze contig
            contig_candidates = self._analyze_contig(contig_id, sequence)
            self.chimera_candidates.extend(contig_candidates)
        
        self.logger.info(f"Detected {len(self.chimera_candidates)} chimera candidates")
        return self.chimera_candidates
    
    def _analyze_contig(self, contig_id: str, sequence: str) -> List[ChimeraCandidate]:
        """Analyze a single contig for chimeric signatures using GC content."""
        # Check if contig is suitable for analysis
        if len(sequence) < self.min_contig_length:
            self.logger.debug(f"Skipping contig {contig_id}: too short ({len(sequence)}bp)")
            return []
        
        # Calculate adaptive window parameters
        window_params = self._calculate_adaptive_window_parameters(sequence)
        
        # Calculate GC content profile
        gc_profile = self._calculate_gc_profile(
            sequence, window_params['window'], window_params['step']
        )
        
        # Detect breakpoints
        breakpoints = self._detect_gc_breakpoints(gc_profile, window_params['step'])
        
        # Create candidates from breakpoints
        candidates = []
        for breakpoint in breakpoints:
            candidate = self._evaluate_breakpoint(
                contig_id, sequence, breakpoint, gc_profile, window_params
            )
            if candidate:
                candidates.append(candidate)
        
        return candidates
    
    def _calculate_adaptive_window_parameters(self, sequence: str) -> dict:
        """Calculate adaptive window and step sizes based on sequence length."""
        contig_length = len(sequence)
        
        # Calculate adaptive window size
        adaptive_window = get_adaptive_window_size(contig_length)
        adaptive_step = get_adaptive_step_size(adaptive_window)
        
        # Adjust window if it's too large for the sequence
        adaptive_window = adjust_window_for_sequence_length(adaptive_window, contig_length)
        
        # Recalculate step size if window was adjusted
        if adaptive_window < get_adaptive_window_size(contig_length):
            adaptive_step = max(50, adaptive_window // 4)
        
        return {
            'length': contig_length,
            'window': adaptive_window,
            'step': adaptive_step
        }
    
    def _calculate_gc_profile(self, sequence: str, window_size: int, step_size: int) -> List[float]:
        """Calculate GC content profile using sliding window."""
        gc_profile = []
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window_seq = sequence[i:i + window_size]
            gc = calculate_gc_content_simple(window_seq)
            gc_profile.append(gc)
        
        return gc_profile
    
    def _detect_gc_breakpoints(self, gc_profile: List[float], step_size: int) -> List[int]:
        """Detect potential breakpoints from GC content changes."""
        if len(gc_profile) < 3:
            return []
        
        # Convert to numpy array and calculate differences between adjacent windows
        gc_array = np.array(gc_profile)
        gc_diffs = np.abs(np.diff(gc_array))
        
        # Debug logging
        self.logger.debug(f"GC profile length: {len(gc_profile)}, type: {type(gc_profile)}")
        self.logger.debug(f"GC array shape: {gc_array.shape}, dtype: {gc_array.dtype}")
        self.logger.debug(f"GC diffs shape: {gc_diffs.shape}, dtype: {gc_diffs.dtype}")
        
        # Ensure gc_diffs is 1D and float
        gc_diffs = gc_diffs.flatten().astype(float)
        
        # Find peaks in GC content changes
        peaks, properties = find_peaks(
            gc_diffs,
            height=self.gc_content_threshold,
            distance=2  # Minimum distance between peaks
        )
        
        # Convert peak indices to genomic positions
        breakpoints = []
        for peak_idx in peaks:
            # Calculate the actual position in the sequence
            # Peak in diff array corresponds to change between windows peak_idx and peak_idx+1
            position = (peak_idx + 1) * step_size
            breakpoints.append(position)
        
        return breakpoints
    
    def _evaluate_breakpoint(self, contig_id: str, sequence: str, breakpoint: int,
                           gc_profile: List[float], window_params: dict) -> Optional[ChimeraCandidate]:
        """Evaluate a potential breakpoint and create a candidate if significant."""
        window_size = window_params['window']
        step_size = window_params['step']
        
        # Calculate flanking region boundaries
        left_start = max(0, breakpoint - window_size)
        left_end = breakpoint
        right_start = breakpoint
        right_end = min(len(sequence), breakpoint + window_size)
        
        # Calculate GC content for flanking regions
        left_seq = sequence[left_start:left_end]
        right_seq = sequence[right_start:right_end]
        
        if len(left_seq) < 100 or len(right_seq) < 100:
            return None
        
        gc_left = calculate_gc_content_simple(left_seq)
        gc_right = calculate_gc_content_simple(right_seq)
        gc_diff = abs(gc_left - gc_right)
        
        # Check if difference is significant
        if gc_diff < self.gc_content_threshold:
            return None
        
        # Calculate confidence score based on GC difference
        confidence = min(1.0, gc_diff / 0.3)  # Normalize to 0-1 range
        
        # Define breakpoint region
        region_start = max(0, breakpoint - step_size)
        region_end = min(len(sequence), breakpoint + step_size)
        
        return ChimeraCandidate(
            contig_id=contig_id,
            breakpoint=breakpoint,
            confidence_score=confidence,
            gc_content_left=gc_left,
            gc_content_right=gc_right,
            gc_difference=gc_diff,
            breakpoint_region=(region_start, region_end)
        )