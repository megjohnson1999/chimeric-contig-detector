#!/usr/bin/env python3
"""
Simple Chimera Detector

A focused tool for detecting chimeric contigs using two reliable signals:
1. GC content shifts
2. Read pair mapping anomalies

Designed for simplicity, reliability, and biological meaningfulness.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import numpy as np
import pysam
from Bio import SeqIO


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ChimeraCandidate:
    """Simple data class for chimera candidates."""
    contig_id: str
    breakpoint: int
    gc_left: float
    gc_right: float
    gc_shift: float
    improper_pair_ratio: float
    confidence: float


class SimpleChimeraDetector:
    """Simple, focused chimera detector using GC shifts and read pair anomalies."""
    
    def __init__(self):
        # Fixed, biologically meaningful parameters
        self.gc_shift_threshold = 0.05      # 5% GC change (detectable shift)
        self.improper_pair_threshold = 0.20  # 20% improper pairs (significant anomaly)
        self.window_size = 300              # 300bp windows (higher resolution)
        self.step_size = 50                 # 50bp steps (more overlap)
        self.edge_exclusion_ratio = 0.2     # Ignore 20% from each end (minimum 500bp)
        self.min_contig_length = 2000       # Only analyze contigs >2kb
        self.consensus_window = 200         # Signals must agree within 200bp
        
        logger.info("Simple Chimera Detector initialized")
        logger.info(f"GC shift threshold: {self.gc_shift_threshold:.1%}")
        logger.info(f"Improper pair threshold: {self.improper_pair_threshold:.1%}")
        logger.info(f"Edge exclusion: {self.edge_exclusion_ratio:.0%} from each end")
    
    def detect_chimeras(self, assembly_file: str, bam_file: str) -> List[ChimeraCandidate]:
        """
        Detect chimeric contigs using GC shifts and read pair anomalies.
        
        Args:
            assembly_file: Path to FASTA assembly
            bam_file: Path to sorted, indexed BAM file
            
        Returns:
            List of chimera candidates with both types of evidence
        """
        logger.info(f"Analyzing assembly: {assembly_file}")
        logger.info(f"Using BAM file: {bam_file}")
        
        candidates = []
        contigs_analyzed = 0
        contigs_skipped = 0
        
        # Load assembly
        sequences = {}
        for record in SeqIO.parse(assembly_file, "fasta"):
            sequences[record.id] = str(record.seq)
        
        logger.info(f"Loaded {len(sequences)} contigs from assembly")
        
        # Analyze each contig
        for contig_id, sequence in sequences.items():
            if len(sequence) < self.min_contig_length:
                contigs_skipped += 1
                continue
                
            logger.debug(f"Analyzing {contig_id} (length: {len(sequence)})")
            
            try:
                contig_candidates = self._analyze_contig(contig_id, sequence, bam_file)
                candidates.extend(contig_candidates)
                contigs_analyzed += 1
                
            except Exception as e:
                logger.error(f"Error analyzing {contig_id}: {e}")
                continue
        
        logger.info(f"Analysis complete: {contigs_analyzed} contigs analyzed, {contigs_skipped} skipped")
        logger.info(f"Found {len(candidates)} chimera candidates")
        
        return candidates
    
    def _analyze_contig(self, contig_id: str, sequence: str, bam_file: str) -> List[ChimeraCandidate]:
        """Analyze a single contig for chimeric signatures."""
        candidates = []
        
        # Find GC content shifts
        gc_breakpoints = self._find_gc_shifts(sequence)
        if not gc_breakpoints:
            return candidates
            
        logger.debug(f"Found {len(gc_breakpoints)} GC shifts in {contig_id}")
        
        # Find read pair anomalies
        read_breakpoints = self._find_read_anomalies(contig_id, sequence, bam_file)
        logger.debug(f"Found {len(read_breakpoints)} read anomalies in {contig_id}")
        if not read_breakpoints:
            return candidates
        
        # Find consensus breakpoints (both signals agree)
        consensus_breakpoints = self._find_consensus_breakpoints(
            gc_breakpoints, read_breakpoints
        )
        
        if not consensus_breakpoints:
            return candidates
            
        logger.debug(f"Found {len(consensus_breakpoints)} consensus breakpoints in {contig_id}")
        
        # Create candidates for consensus breakpoints
        for breakpoint in consensus_breakpoints:
            candidate = self._create_candidate(contig_id, sequence, breakpoint, bam_file)
            if candidate:
                candidates.append(candidate)
        
        return candidates
    
    def _find_gc_shifts(self, sequence: str) -> List[int]:
        """Find positions with significant GC content shifts."""
        breakpoints = []
        seq_len = len(sequence)
        
        # Calculate adaptive edge exclusion
        edge_exclusion = max(500, int(seq_len * self.edge_exclusion_ratio))
        edge_exclusion = min(edge_exclusion, seq_len // 4)  # Never exclude more than 25%
        
        # Calculate GC content in sliding windows
        gc_values = []
        positions = []
        
        logger.debug(f"GC analysis range: {edge_exclusion} to {seq_len - edge_exclusion - self.window_size}")
        
        for i in range(edge_exclusion, seq_len - edge_exclusion - self.window_size, self.step_size):
            window = sequence[i:i + self.window_size]
            gc_content = (window.count('G') + window.count('C')) / len(window)
            gc_values.append(gc_content)
            positions.append(i + self.window_size // 2)  # Center of window
            logger.debug(f"GC window {i}-{i+self.window_size}: {gc_content:.3f}")
        
        if len(gc_values) < 2:
            return breakpoints
        
        # Find significant changes between adjacent windows
        for i in range(len(gc_values) - 1):
            gc_shift = abs(gc_values[i + 1] - gc_values[i])
            
            logger.debug(f"Window {positions[i]} -> {positions[i+1]}: GC {gc_values[i]:.3f} -> {gc_values[i+1]:.3f}, shift={gc_shift:.3f}")
            
            if gc_shift >= self.gc_shift_threshold:
                # Breakpoint is between the two windows
                breakpoint = (positions[i] + positions[i + 1]) // 2
                breakpoints.append(breakpoint)
                logger.debug(f"*** GC shift of {gc_shift:.3f} at position {breakpoint}")
        
        return breakpoints
    
    def _find_read_anomalies(self, contig_id: str, sequence: str, bam_file: str) -> List[int]:
        """Find positions with read pair mapping anomalies."""
        breakpoints = []
        seq_len = len(sequence)
        
        # Calculate adaptive edge exclusion
        edge_exclusion = max(500, int(seq_len * self.edge_exclusion_ratio))
        edge_exclusion = min(edge_exclusion, seq_len // 4)  # Never exclude more than 25%
        
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                if contig_id not in bam.references:
                    logger.debug(f"Contig {contig_id} not found in BAM file")
                    return breakpoints
                
                # Calculate proper pair ratios in sliding windows
                ratios = []
                positions = []
                
                for i in range(edge_exclusion, seq_len - edge_exclusion - self.window_size, self.step_size):
                    window_start = i
                    window_end = i + self.window_size
                    
                    proper_pairs = 0
                    improper_pairs = 0
                    
                    # Count reads in this window
                    for read in bam.fetch(contig_id, window_start, window_end):
                        if (read.is_paired and 
                            not read.is_unmapped and 
                            not read.mate_is_unmapped and
                            read.reference_start >= window_start and
                            read.reference_start < window_end):
                            
                            if read.is_proper_pair:
                                proper_pairs += 1
                            else:
                                improper_pairs += 1
                    
                    total_pairs = proper_pairs + improper_pairs
                    if total_pairs >= 10:  # Minimum coverage for reliable calculation
                        improper_ratio = improper_pairs / total_pairs
                        ratios.append(improper_ratio)
                        positions.append(i + self.window_size // 2)
                
                # Find significant increases in improper pair ratio
                for i in range(len(ratios) - 1):
                    if (ratios[i + 1] >= self.improper_pair_threshold and 
                        ratios[i + 1] > ratios[i] + 0.1):  # At least 10% increase
                        
                        breakpoint = (positions[i] + positions[i + 1]) // 2
                        breakpoints.append(breakpoint)
                        logger.debug(f"Read anomaly: {ratios[i + 1]:.3f} improper pairs at position {breakpoint}")
        
        except Exception as e:
            logger.error(f"Error analyzing read pairs for {contig_id}: {e}")
            return []
        
        return breakpoints
    
    def _find_consensus_breakpoints(self, gc_breakpoints: List[int], read_breakpoints: List[int]) -> List[int]:
        """Find breakpoints where both GC and read signals agree."""
        consensus = []
        
        for gc_bp in gc_breakpoints:
            for read_bp in read_breakpoints:
                if abs(gc_bp - read_bp) <= self.consensus_window:
                    # Use the average position as the consensus breakpoint
                    consensus_bp = (gc_bp + read_bp) // 2
                    consensus.append(consensus_bp)
                    logger.debug(f"Consensus breakpoint at {consensus_bp} (GC: {gc_bp}, reads: {read_bp})")
                    break  # Only use each GC breakpoint once
        
        return consensus
    
    def _create_candidate(self, contig_id: str, sequence: str, breakpoint: int, bam_file: str) -> Optional[ChimeraCandidate]:
        """Create a chimera candidate with detailed information."""
        try:
            # Calculate GC content on left and right sides
            left_start = max(0, breakpoint - 500)
            left_seq = sequence[left_start:breakpoint]
            
            right_end = min(len(sequence), breakpoint + 500)
            right_seq = sequence[breakpoint:right_end]
            
            if len(left_seq) < 100 or len(right_seq) < 100:
                return None
            
            gc_left = (left_seq.count('G') + left_seq.count('C')) / len(left_seq)
            gc_right = (right_seq.count('G') + right_seq.count('C')) / len(right_seq)
            gc_shift = abs(gc_right - gc_left)
            
            # Calculate improper pair ratio around breakpoint
            improper_ratio = self._calculate_improper_pair_ratio(contig_id, breakpoint, bam_file)
            
            # Calculate confidence score
            confidence = self._calculate_confidence(gc_shift, improper_ratio)
            
            return ChimeraCandidate(
                contig_id=contig_id,
                breakpoint=breakpoint,
                gc_left=gc_left,
                gc_right=gc_right,
                gc_shift=gc_shift,
                improper_pair_ratio=improper_ratio,
                confidence=confidence
            )
            
        except Exception as e:
            logger.error(f"Error creating candidate for {contig_id} at {breakpoint}: {e}")
            return None
    
    def _calculate_improper_pair_ratio(self, contig_id: str, breakpoint: int, bam_file: str) -> float:
        """Calculate improper pair ratio around a breakpoint."""
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                proper_pairs = 0
                improper_pairs = 0
                window = 300
                
                start = max(0, breakpoint - window)
                end = breakpoint + window
                
                for read in bam.fetch(contig_id, start, end):
                    if (read.is_paired and 
                        not read.is_unmapped and 
                        not read.mate_is_unmapped):
                        
                        if read.is_proper_pair:
                            proper_pairs += 1
                        else:
                            improper_pairs += 1
                
                total = proper_pairs + improper_pairs
                return improper_pairs / total if total > 0 else 0.0
                
        except Exception:
            return 0.0
    
    def _calculate_confidence(self, gc_shift: float, improper_ratio: float) -> float:
        """Calculate confidence score based on signal strength."""
        # Normalize signals to 0-1 range
        gc_score = min(1.0, gc_shift / 0.3)  # Max at 30% shift
        read_score = min(1.0, improper_ratio / 0.5)  # Max at 50% improper
        
        # Simple average
        return (gc_score + read_score) / 2
    
    def split_contigs(self, assembly_file: str, candidates: List[ChimeraCandidate], output_dir: str):
        """Split contigs at validated breakpoints."""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Load original assembly
        sequences = {}
        for record in SeqIO.parse(assembly_file, "fasta"):
            sequences[record.id] = str(record.seq)
        
        # Group candidates by contig (take highest confidence if multiple)
        contig_breakpoints = {}
        for candidate in candidates:
            if (candidate.contig_id not in contig_breakpoints or 
                candidate.confidence > contig_breakpoints[candidate.contig_id].confidence):
                contig_breakpoints[candidate.contig_id] = candidate
        
        # Write cleaned assembly
        cleaned_file = output_path / "cleaned_assembly.fasta"
        with open(cleaned_file, 'w') as f:
            for contig_id, sequence in sequences.items():
                if contig_id in contig_breakpoints:
                    # Split this contig
                    candidate = contig_breakpoints[contig_id]
                    breakpoint = candidate.breakpoint
                    
                    left_seq = sequence[:breakpoint]
                    right_seq = sequence[breakpoint:]
                    
                    # Only split if both parts are substantial
                    if len(left_seq) >= 500 and len(right_seq) >= 500:
                        f.write(f">{contig_id}_split_A\n{left_seq}\n")
                        f.write(f">{contig_id}_split_B\n{right_seq}\n")
                        logger.info(f"Split {contig_id} at position {breakpoint} "
                                  f"(confidence: {candidate.confidence:.3f})")
                    else:
                        f.write(f">{contig_id}\n{sequence}\n")
                        logger.warning(f"Skipped splitting {contig_id} - resulting segments too short")
                else:
                    # Keep original contig
                    f.write(f">{contig_id}\n{sequence}\n")
        
        logger.info(f"Cleaned assembly written to {cleaned_file}")
    
    def write_report(self, candidates: List[ChimeraCandidate], output_file: str):
        """Write a simple TSV report of chimera candidates."""
        with open(output_file, 'w') as f:
            f.write("contig_id\tbreakpoint\tgc_left\tgc_right\tgc_shift\t"
                   "improper_pair_ratio\tconfidence\n")
            
            for candidate in sorted(candidates, key=lambda x: x.confidence, reverse=True):
                f.write(f"{candidate.contig_id}\t{candidate.breakpoint}\t"
                       f"{candidate.gc_left:.3f}\t{candidate.gc_right:.3f}\t"
                       f"{candidate.gc_shift:.3f}\t{candidate.improper_pair_ratio:.3f}\t"
                       f"{candidate.confidence:.3f}\n")
        
        logger.info(f"Report written to {output_file}")


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Simple chimera detector using GC shifts and read pair anomalies"
    )
    parser.add_argument("-a", "--assembly", required=True,
                       help="Path to assembly FASTA file")
    parser.add_argument("-b", "--bam", required=True,
                       help="Path to sorted, indexed BAM file")
    parser.add_argument("-o", "--output", required=True,
                       help="Output directory")
    parser.add_argument("--debug", action="store_true",
                       help="Enable debug logging")
    
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate input files
    if not Path(args.assembly).exists():
        logger.error(f"Assembly file not found: {args.assembly}")
        sys.exit(1)
    
    if not Path(args.bam).exists():
        logger.error(f"BAM file not found: {args.bam}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    # Run detection
    detector = SimpleChimeraDetector()
    candidates = detector.detect_chimeras(args.assembly, args.bam)
    
    if candidates:
        logger.info(f"Detected {len(candidates)} chimeric contigs")
        
        # Write report
        report_file = output_dir / "chimera_report.tsv"
        detector.write_report(candidates, str(report_file))
        
        # Split contigs
        detector.split_contigs(args.assembly, candidates, str(output_dir))
        
    else:
        logger.info("No chimeric contigs detected")
        
        # Still create cleaned assembly (same as input)
        import shutil
        shutil.copy(args.assembly, output_dir / "cleaned_assembly.fasta")
    
    logger.info("Analysis complete!")


if __name__ == "__main__":
    main()