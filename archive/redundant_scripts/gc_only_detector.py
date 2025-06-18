#!/usr/bin/env python3
"""
GC-Only Chimera Detector

Simplified version that detects chimeric contigs using only GC content shifts.
This version doesn't require read pair analysis and can work with assembly-only data.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
import numpy as np
from Bio import SeqIO


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ChimeraCandidate:
    """Data class for GC-based chimera candidates."""
    contig_id: str
    breakpoint: int
    gc_left: float
    gc_right: float
    gc_shift: float
    confidence: float


class GCOnlyDetector:
    """Simple chimera detector using only GC content shifts."""
    
    def __init__(self):
        # Tuned parameters based on testing
        self.gc_shift_threshold = 0.08      # 8% GC change 
        self.window_size = 300              # 300bp windows
        self.step_size = 50                 # 50bp steps
        self.edge_exclusion_ratio = 0.15    # Ignore 15% from each end
        self.min_contig_length = 2000       # Only analyze contigs >2kb
        
        logger.info("GC-Only Chimera Detector initialized")
        logger.info(f"GC shift threshold: {self.gc_shift_threshold:.1%}")
        logger.info(f"Window size: {self.window_size}bp")
        logger.info(f"Edge exclusion: {self.edge_exclusion_ratio:.0%} from each end")
    
    def detect_chimeras(self, assembly_file: str) -> List[ChimeraCandidate]:
        """Detect chimeric contigs using GC content shifts only."""
        logger.info(f"Analyzing assembly: {assembly_file}")
        
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
                contig_candidates = self._analyze_contig(contig_id, sequence)
                candidates.extend(contig_candidates)
                contigs_analyzed += 1
                
            except Exception as e:
                logger.error(f"Error analyzing {contig_id}: {e}")
                continue
        
        logger.info(f"Analysis complete: {contigs_analyzed} contigs analyzed, {contigs_skipped} skipped")
        logger.info(f"Found {len(candidates)} chimera candidates")
        
        return candidates
    
    def _analyze_contig(self, contig_id: str, sequence: str) -> List[ChimeraCandidate]:
        """Analyze a single contig for GC content shifts."""
        candidates = []
        
        # Find GC content shifts
        gc_breakpoints = self._find_gc_shifts(sequence)
        logger.debug(f"Found {len(gc_breakpoints)} GC shifts in {contig_id}")
        
        if not gc_breakpoints:
            return candidates
        
        # Create candidates for each breakpoint
        for breakpoint in gc_breakpoints:
            candidate = self._create_candidate(contig_id, sequence, breakpoint)
            if candidate:
                candidates.append(candidate)
        
        return candidates
    
    def _find_gc_shifts(self, sequence: str) -> List[int]:
        """Find positions with significant GC content shifts."""
        breakpoints = []
        seq_len = len(sequence)
        
        # Calculate adaptive edge exclusion
        edge_exclusion = max(300, int(seq_len * self.edge_exclusion_ratio))
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
        
        if len(gc_values) < 2:
            return breakpoints
        
        # Find significant changes between adjacent windows
        for i in range(len(gc_values) - 1):
            gc_shift = abs(gc_values[i + 1] - gc_values[i])
            
            if gc_shift >= self.gc_shift_threshold:
                # Breakpoint is between the two windows
                breakpoint = (positions[i] + positions[i + 1]) // 2
                breakpoints.append(breakpoint)
                logger.debug(f"*** GC shift of {gc_shift:.3f} at position {breakpoint}")
        
        return breakpoints
    
    def _create_candidate(self, contig_id: str, sequence: str, breakpoint: int) -> Optional[ChimeraCandidate]:
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
            
            # Simple confidence based on GC shift magnitude
            confidence = min(1.0, gc_shift / 0.3)  # Max confidence at 30% shift
            
            return ChimeraCandidate(
                contig_id=contig_id,
                breakpoint=breakpoint,
                gc_left=gc_left,
                gc_right=gc_right,
                gc_shift=gc_shift,
                confidence=confidence
            )
            
        except Exception as e:
            logger.error(f"Error creating candidate for {contig_id} at {breakpoint}: {e}")
            return None
    
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
                                  f"(GC shift: {candidate.gc_shift:.1%}, confidence: {candidate.confidence:.3f})")
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
            f.write("contig_id\tbreakpoint\tgc_left\tgc_right\tgc_shift\tconfidence\n")
            
            for candidate in sorted(candidates, key=lambda x: x.confidence, reverse=True):
                f.write(f"{candidate.contig_id}\t{candidate.breakpoint}\t"
                       f"{candidate.gc_left:.3f}\t{candidate.gc_right:.3f}\t"
                       f"{candidate.gc_shift:.3f}\t{candidate.confidence:.3f}\n")
        
        logger.info(f"Report written to {output_file}")


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="GC-only chimera detector using GC content shifts"
    )
    parser.add_argument("-a", "--assembly", required=True,
                       help="Path to assembly FASTA file")
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
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    # Run detection
    detector = GCOnlyDetector()
    candidates = detector.detect_chimeras(args.assembly)
    
    if candidates:
        logger.info(f"Detected {len(candidates)} chimeric contigs")
        
        # Write report
        report_file = output_dir / "gc_chimera_report.tsv"
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