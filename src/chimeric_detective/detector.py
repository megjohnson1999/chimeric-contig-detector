"""
Chimera detection module for identifying potential chimeric contigs.
"""

import logging
import os
import tempfile
from typing import List, Dict, Tuple, Optional, Set
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
                       temp_dir: Optional[str] = None) -> List[ChimeraCandidate]:
        """
        Main method to detect chimeric contigs.
        
        Args:
            assembly_file: Path to assembly FASTA file
            bam_file: Path to BAM file with aligned reads
            reads1: Path to forward reads (if BAM not provided)
            reads2: Path to reverse reads (if BAM not provided)
            temp_dir: Temporary directory for intermediate files
            
        Returns:
            List of ChimeraCandidate objects
        """
        self.logger.info("Starting chimera detection")
        
        # Load assembly
        contigs = self._load_assembly(assembly_file)
        self.logger.info(f"Loaded {len(contigs)} contigs from assembly")
        
        # Prepare BAM file
        if bam_file is None:
            if reads1 is None:
                raise ValueError("Either bam_file or reads1 must be provided")
            bam_file = self._align_reads(assembly_file, reads1, reads2, temp_dir)
        
        # Analyze each contig
        candidates = []
        for contig_id, sequence in contigs.items():
            if len(sequence) < self.min_contig_length:
                continue
                
            self.logger.debug(f"Analyzing contig {contig_id}")
            contig_candidates = self._analyze_contig(contig_id, sequence, bam_file)
            candidates.extend(contig_candidates)
        
        self.chimera_candidates = candidates
        self.logger.info(f"Detected {len(candidates)} chimera candidates")
        
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
        
        # Try minimap2 first, then BWA
        try:
            if reads2:
                # Paired-end
                cmd = ["minimap2", "-ax", "sr", assembly_file, reads1, reads2]
            else:
                # Single-end
                cmd = ["minimap2", "-ax", "sr", assembly_file, reads1]
            
            with open(sam_file, 'w') as f:
                result = run_command(cmd, capture_output=True)
                f.write(result.stdout)
            
        except Exception as e:
            self.logger.warning(f"minimap2 failed: {e}, trying BWA")
            
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
            except Exception as bwa_error:
                self.logger.error(f"Both minimap2 and BWA failed. minimap2 error: {e}, BWA error: {bwa_error}")
                raise RuntimeError("No aligner available. Please install minimap2 or BWA.")
        
        # Convert to BAM and sort
        run_command(["samtools", "view", "-bS", str(sam_file), "-o", str(bam_file)])
        run_command(["samtools", "sort", str(bam_file), "-o", str(sorted_bam)])
        run_command(["samtools", "index", str(sorted_bam)])
        
        return str(sorted_bam)
    
    def _analyze_contig(self, contig_id: str, sequence: str, bam_file: str) -> List[ChimeraCandidate]:
        """Analyze a single contig for chimeric signatures."""
        candidates = []
        
        # Calculate coverage profile
        coverage = self._calculate_coverage(contig_id, len(sequence), bam_file)
        
        # Calculate GC content profile
        gc_profile = calculate_gc_content(sequence, self.window_size)
        
        # Calculate k-mer profiles
        kmer_profile = self._calculate_kmer_profile(sequence)
        
        # Detect breakpoints using multiple methods
        coverage_breakpoints = self._detect_coverage_breakpoints(coverage)
        gc_breakpoints = self._detect_gc_breakpoints(gc_profile)
        kmer_breakpoints = self._detect_kmer_breakpoints(kmer_profile)
        
        # Combine and score breakpoints
        all_breakpoints = set(coverage_breakpoints + gc_breakpoints + kmer_breakpoints)
        
        for breakpoint in all_breakpoints:
            candidate = self._evaluate_breakpoint(
                contig_id, sequence, breakpoint, bam_file,
                coverage, gc_profile, kmer_profile
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
    
    def _detect_coverage_breakpoints(self, coverage: np.ndarray) -> List[int]:
        """Detect potential breakpoints based on coverage discontinuities."""
        breakpoints = []
        
        # Smooth coverage to reduce noise
        window = np.ones(self.window_size) / self.window_size
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
        """Detect potential breakpoints based on GC content changes."""
        breakpoints = []
        
        if len(gc_profile) < 4:
            return breakpoints
        
        # Convert to positions in the sequence
        step = self.window_size // 2
        
        for i in range(1, len(gc_profile) - 1):
            left_gc = gc_profile[i-1]
            right_gc = gc_profile[i+1]
            
            gc_diff = abs(left_gc - right_gc)
            if gc_diff >= self.gc_content_threshold:
                position = i * step
                breakpoints.append(position)
        
        return breakpoints
    
    def _detect_kmer_breakpoints(self, kmer_profile: List[Dict[str, int]]) -> List[int]:
        """Detect potential breakpoints based on k-mer composition changes."""
        breakpoints = []
        
        if len(kmer_profile) < 4:
            return breakpoints
        
        step = self.step_size
        
        for i in range(1, len(kmer_profile) - 1):
            left_kmers = kmer_profile[i-1]
            right_kmers = kmer_profile[i+1]
            
            kmer_dist = calculate_kmer_distance(left_kmers, right_kmers)
            if kmer_dist >= self.kmer_distance_threshold:
                position = i * step
                breakpoints.append(position)
        
        return breakpoints
    
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