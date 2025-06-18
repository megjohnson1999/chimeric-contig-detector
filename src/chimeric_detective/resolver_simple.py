"""
Simplified chimera resolution module for splitting contigs at GC breakpoints.
"""

import logging
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .analyzer_simple import ChimeraClassification
from .utils import setup_logging


class SimpleChimeraResolver:
    """Simplified resolver for splitting chimeric contigs."""
    
    def __init__(self,
                 min_fragment_length: int = 500,
                 output_dir: str = "resolved_contigs",
                 log_level: str = "INFO"):
        """
        Initialize SimpleChimeraResolver.
        
        Args:
            min_fragment_length: Minimum length for split fragments
            output_dir: Directory for output files
            log_level: Logging level
        """
        self.min_fragment_length = min_fragment_length
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = setup_logging(log_level)
    
    def resolve_chimeras(self, 
                        assembly_file: str,
                        classifications: List[ChimeraClassification],
                        output_prefix: str = "resolved") -> Dict[str, any]:
        """
        Resolve chimeric contigs by splitting at breakpoints.
        
        Args:
            assembly_file: Path to original assembly file
            classifications: List of chimera classifications
            output_prefix: Prefix for output files
            
        Returns:
            Dictionary with resolution statistics
        """
        self.logger.info(f"Resolving chimeras in {assembly_file}")
        
        # Filter to only include chimeras that should be split
        chimeras_to_split = [c for c in classifications 
                           if c.is_likely_chimera and c.classification == 'technical_artifact']
        
        self.logger.info(f"Found {len(chimeras_to_split)} chimeras to split")
        
        # Load sequences
        sequences = {record.id: record for record in SeqIO.parse(assembly_file, "fasta")}
        
        # Process chimeras
        split_contigs = []
        kept_contigs = []
        statistics = {
            'total_contigs': len(sequences),
            'chimeras_found': len(chimeras_to_split),
            'contigs_split': 0,
            'fragments_created': 0,
            'fragments_discarded': 0
        }
        
        # Group classifications by contig
        chimeras_by_contig = {}
        for classification in chimeras_to_split:
            contig_id = classification.candidate.contig_id
            if contig_id not in chimeras_by_contig:
                chimeras_by_contig[contig_id] = []
            chimeras_by_contig[contig_id].append(classification)
        
        # Process each contig
        for contig_id, record in sequences.items():
            if contig_id in chimeras_by_contig:
                # Split this contig
                fragments = self._split_contig(record, chimeras_by_contig[contig_id])
                valid_fragments = [f for f in fragments if len(f.seq) >= self.min_fragment_length]
                
                if len(valid_fragments) > 1:
                    split_contigs.extend(valid_fragments)
                    statistics['contigs_split'] += 1
                    statistics['fragments_created'] += len(valid_fragments)
                    statistics['fragments_discarded'] += len(fragments) - len(valid_fragments)
                else:
                    # Keep original if splitting didn't produce valid fragments
                    kept_contigs.append(record)
            else:
                # Keep non-chimeric contigs
                kept_contigs.append(record)
        
        # Write output files
        output_file = self.output_dir / f"{output_prefix}_assembly.fasta"
        all_contigs = kept_contigs + split_contigs
        SeqIO.write(all_contigs, output_file, "fasta")
        
        # Write splitting report
        report_file = self.output_dir / f"{output_prefix}_splitting_report.txt"
        self._write_report(report_file, chimeras_to_split, statistics)
        
        statistics['output_file'] = str(output_file)
        statistics['report_file'] = str(report_file)
        statistics['final_contig_count'] = len(all_contigs)
        
        return statistics
    
    def _split_contig(self, record: SeqRecord, classifications: List[ChimeraClassification]) -> List[SeqRecord]:
        """Split a contig at breakpoint positions."""
        # Sort breakpoints by position
        breakpoints = sorted([c.candidate.breakpoint for c in classifications])
        
        # Add start and end positions
        positions = [0] + breakpoints + [len(record.seq)]
        
        # Create fragments
        fragments = []
        for i in range(len(positions) - 1):
            start = positions[i]
            end = positions[i + 1]
            
            fragment_seq = record.seq[start:end]
            fragment_id = f"{record.id}_fragment_{i+1}"
            fragment_desc = f"Split from {record.id} at positions {breakpoints}"
            
            fragment_record = SeqRecord(
                fragment_seq,
                id=fragment_id,
                description=fragment_desc
            )
            fragments.append(fragment_record)
        
        return fragments
    
    def _write_report(self, report_file: Path, chimeras: List[ChimeraClassification], 
                     statistics: Dict[str, any]):
        """Write a summary report of the resolution process."""
        with open(report_file, 'w') as f:
            f.write("Chimera Resolution Report\n")
            f.write("=" * 50 + "\n\n")
            
            # Write statistics
            f.write("Summary Statistics:\n")
            for key, value in statistics.items():
                if not key.endswith('_file'):
                    f.write(f"  {key}: {value}\n")
            f.write("\n")
            
            # Write details of each split
            f.write("Chimeras Split:\n")
            f.write("-" * 50 + "\n")
            for classification in chimeras:
                candidate = classification.candidate
                f.write(f"\nContig: {candidate.contig_id}\n")
                f.write(f"  Breakpoint: {candidate.breakpoint}\n")
                f.write(f"  GC Left: {candidate.gc_content_left:.3f}\n")
                f.write(f"  GC Right: {candidate.gc_content_right:.3f}\n")
                f.write(f"  GC Difference: {candidate.gc_difference:.3f}\n")
                f.write(f"  Confidence: {classification.confidence:.3f}\n")
                f.write(f"  Reason: {classification.reason}\n")