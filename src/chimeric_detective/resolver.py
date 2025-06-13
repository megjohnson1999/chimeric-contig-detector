"""
Chimera resolution module for splitting technical chimeras and cleaning assemblies.
"""

import logging
import os
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from dataclasses import dataclass
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .analyzer import ChimeraAnalysis
from .utils import setup_logging, calculate_n50


@dataclass
class SplittingDecision:
    """Data class for storing splitting decisions."""
    contig_id: str
    original_length: int
    breakpoint: int
    action: str  # 'split', 'preserve', 'flag'
    reason: str
    new_contig_ids: List[str]
    new_lengths: List[int]
    confidence: float


class ChimeraResolver:
    """Resolver for processing chimera analysis results and creating cleaned assemblies."""
    
    def __init__(self,
                 split_technical: bool = True,
                 split_pcr: bool = True,
                 preserve_biological: bool = True,
                 min_split_length: int = 500,
                 confidence_threshold: float = 0.5,
                 log_level: str = "INFO"):
        """
        Initialize ChimeraResolver.
        
        Args:
            split_technical: Whether to split technical artifacts
            split_pcr: Whether to split PCR chimeras
            preserve_biological: Whether to preserve biological recombination
            min_split_length: Minimum length for split contigs
            confidence_threshold: Minimum confidence for splitting
            log_level: Logging level
        """
        self.split_technical = split_technical
        self.split_pcr = split_pcr
        self.preserve_biological = preserve_biological
        self.min_split_length = min_split_length
        self.confidence_threshold = confidence_threshold
        
        self.logger = setup_logging(log_level)
        self.splitting_decisions: List[SplittingDecision] = []
    
    def resolve_chimeras(self,
                        analyses: List[ChimeraAnalysis],
                        assembly_file: str,
                        output_dir: str) -> Dict[str, str]:
        """
        Process chimera analyses and create cleaned assembly.
        
        Args:
            analyses: List of ChimeraAnalysis objects
            assembly_file: Path to original assembly file
            output_dir: Output directory for results
            
        Returns:
            Dictionary with paths to output files
        """
        self.logger.info(f"Resolving {len(analyses)} chimera analyses")
        
        # Load original assembly
        original_contigs = self._load_assembly(assembly_file)
        
        # Make splitting decisions
        decisions = self._make_splitting_decisions(analyses)
        self.splitting_decisions = decisions
        
        # Create cleaned assembly
        cleaned_contigs = self._create_cleaned_assembly(original_contigs, decisions)
        
        # Write output files
        output_files = self._write_output_files(
            cleaned_contigs, decisions, analyses, output_dir
        )
        
        # Generate summary statistics
        summary_stats = self._generate_summary_stats(
            original_contigs, cleaned_contigs, decisions
        )
        
        self.logger.info(f"Resolution complete. {summary_stats['contigs_split']} contigs split, "
                        f"{summary_stats['contigs_preserved']} preserved")
        
        return output_files
    
    def _load_assembly(self, assembly_file: str) -> Dict[str, str]:
        """Load assembly sequences from FASTA file."""
        contigs = {}
        for record in SeqIO.parse(assembly_file, "fasta"):
            contigs[record.id] = str(record.seq)
        return contigs
    
    def _make_splitting_decisions(self, analyses: List[ChimeraAnalysis]) -> List[SplittingDecision]:
        """Make decisions about which chimeras to split."""
        decisions = []
        
        # Group analyses by contig (a contig might have multiple chimera candidates)
        contig_analyses = {}
        for analysis in analyses:
            contig_id = analysis.candidate.contig_id
            if contig_id not in contig_analyses:
                contig_analyses[contig_id] = []
            contig_analyses[contig_id].append(analysis)
        
        for contig_id, contig_analyses_list in contig_analyses.items():
            # Sort by confidence score (highest first)
            contig_analyses_list.sort(key=lambda x: x.classification_confidence, reverse=True)
            
            # Process the highest confidence analysis for this contig
            best_analysis = contig_analyses_list[0]
            decision = self._make_single_decision(best_analysis)
            decisions.append(decision)
        
        return decisions
    
    def _make_single_decision(self, analysis: ChimeraAnalysis) -> SplittingDecision:
        """Make splitting decision for a single chimera analysis."""
        
        candidate = analysis.candidate
        chimera_type = analysis.chimera_type
        confidence = analysis.classification_confidence
        recommendation = analysis.recommendation
        
        # Default decision
        action = 'preserve'
        reason = 'Default preservation'
        new_contig_ids = [candidate.contig_id]
        new_lengths = [0]  # Will be filled later
        
        # Apply decision rules
        if confidence >= self.confidence_threshold:
            if chimera_type == 'technical_artifact' and self.split_technical:
                if recommendation.startswith('SPLIT'):
                    action = 'split'
                    reason = f'Technical artifact (confidence: {confidence:.2f})'
                    new_contig_ids = [
                        f"{candidate.contig_id}_split_A",
                        f"{candidate.contig_id}_split_B"
                    ]
            
            elif chimera_type == 'pcr_chimera' and self.split_pcr:
                if recommendation.startswith('SPLIT'):
                    action = 'split'
                    reason = f'PCR chimera (confidence: {confidence:.2f})'
                    new_contig_ids = [
                        f"{candidate.contig_id}_split_A",
                        f"{candidate.contig_id}_split_B"
                    ]
            
            elif chimera_type == 'biological_recombination' and self.preserve_biological:
                action = 'flag'
                reason = f'Biological recombination preserved (confidence: {confidence:.2f})'
        
        else:
            action = 'preserve'
            reason = f'Low confidence ({confidence:.2f}), preserved for safety'
        
        return SplittingDecision(
            contig_id=candidate.contig_id,
            original_length=0,  # Will be filled later
            breakpoint=candidate.breakpoint,
            action=action,
            reason=reason,
            new_contig_ids=new_contig_ids,
            new_lengths=new_lengths,
            confidence=confidence
        )
    
    def _create_cleaned_assembly(self,
                               original_contigs: Dict[str, str],
                               decisions: List[SplittingDecision]) -> Dict[str, str]:
        """Create cleaned assembly based on splitting decisions."""
        cleaned_contigs = {}
        
        # Create decision lookup
        decision_lookup = {d.contig_id: d for d in decisions}
        
        for contig_id, sequence in original_contigs.items():
            if contig_id in decision_lookup:
                decision = decision_lookup[contig_id]
                decision.original_length = len(sequence)
                
                if decision.action == 'split':
                    # Split the contig
                    left_seq = sequence[:decision.breakpoint]
                    right_seq = sequence[decision.breakpoint:]
                    
                    # Check minimum length requirements
                    if len(left_seq) >= self.min_split_length and len(right_seq) >= self.min_split_length:
                        cleaned_contigs[decision.new_contig_ids[0]] = left_seq
                        cleaned_contigs[decision.new_contig_ids[1]] = right_seq
                        decision.new_lengths = [len(left_seq), len(right_seq)]
                        
                        self.logger.info(f"Split {contig_id} at position {decision.breakpoint} "
                                       f"into {len(left_seq)} bp and {len(right_seq)} bp contigs")
                    else:
                        # Don't split if resulting contigs would be too short
                        cleaned_contigs[contig_id] = sequence
                        decision.action = 'preserve'
                        decision.reason += ' (split contigs too short)'
                        decision.new_contig_ids = [contig_id]
                        decision.new_lengths = [len(sequence)]
                        
                        self.logger.warning(f"Skipped splitting {contig_id} - resulting contigs too short")
                
                else:
                    # Preserve the contig
                    cleaned_contigs[contig_id] = sequence
                    decision.new_lengths = [len(sequence)]
            
            else:
                # Contig not analyzed, preserve as-is
                cleaned_contigs[contig_id] = sequence
        
        return cleaned_contigs
    
    def _write_output_files(self,
                          cleaned_contigs: Dict[str, str],
                          decisions: List[SplittingDecision],
                          analyses: List[ChimeraAnalysis],
                          output_dir: str) -> Dict[str, str]:
        """Write all output files."""
        
        output_dir = Path(output_dir)
        output_files = {}
        
        # Write cleaned assembly
        cleaned_assembly_path = output_dir / "cleaned_assembly.fasta"
        self._write_fasta(cleaned_contigs, cleaned_assembly_path)
        output_files['cleaned_assembly'] = str(cleaned_assembly_path)
        
        # Write individual chimeric contigs
        chimeric_dir = output_dir / "chimeric_contigs"
        self._write_individual_contigs(decisions, analyses, chimeric_dir)
        output_files['chimeric_contigs_dir'] = str(chimeric_dir)
        
        # Write splitting decisions table
        decisions_path = output_dir / "splitting_decisions.tsv"
        self._write_decisions_table(decisions, decisions_path)
        output_files['splitting_decisions'] = str(decisions_path)
        
        # Write detailed results JSON
        results_path = output_dir / "chimeric_detective_results.json"
        self._write_results_json(analyses, decisions, results_path)
        output_files['results_json'] = str(results_path)
        
        return output_files
    
    def _write_fasta(self, contigs: Dict[str, str], output_path: Path):
        """Write contigs to FASTA file."""
        records = []
        for contig_id, sequence in contigs.items():
            record = SeqRecord(
                Seq(sequence),
                id=contig_id,
                description=f"length={len(sequence)}"
            )
            records.append(record)
        
        with open(output_path, 'w') as f:
            SeqIO.write(records, f, "fasta")
    
    def _write_individual_contigs(self,
                                decisions: List[SplittingDecision],
                                analyses: List[ChimeraAnalysis],
                                output_dir: Path):
        """Write individual files for each chimeric contig."""
        
        output_dir.mkdir(exist_ok=True)
        
        # Create analysis lookup
        analysis_lookup = {a.candidate.contig_id: a for a in analyses}
        
        for decision in decisions:
            if decision.contig_id in analysis_lookup:
                analysis = analysis_lookup[decision.contig_id]
                
                # Write original chimeric contig
                original_path = output_dir / f"{decision.contig_id}_chimeric.fasta"
                with open(original_path, 'w') as f:
                    f.write(f">{decision.contig_id}_chimeric\n")
                    f.write(f"# Chimera type: {analysis.chimera_type}\n")
                    f.write(f"# Confidence: {analysis.classification_confidence:.3f}\n")
                    f.write(f"# Breakpoint: {decision.breakpoint}\n")
                    f.write(f"# Action: {decision.action}\n")
                    f.write(f"# Explanation: {analysis.explanation}\n")
                
                # Write split contigs if applicable
                if decision.action == 'split' and len(decision.new_contig_ids) == 2:
                    # This would require access to the original sequence
                    # For now, just create placeholder files
                    for i, new_id in enumerate(decision.new_contig_ids):
                        split_path = output_dir / f"{new_id}.fasta"
                        with open(split_path, 'w') as f:
                            f.write(f">{new_id}\n")
                            f.write(f"# Split from {decision.contig_id} at position {decision.breakpoint}\n")
                            f.write(f"# Length: {decision.new_lengths[i] if i < len(decision.new_lengths) else 'unknown'}\n")
    
    def _write_decisions_table(self, decisions: List[SplittingDecision], output_path: Path):
        """Write splitting decisions to TSV file."""
        
        rows = []
        for decision in decisions:
            row = {
                'contig_id': decision.contig_id,
                'original_length': decision.original_length,
                'breakpoint': decision.breakpoint,
                'action': decision.action,
                'reason': decision.reason,
                'confidence': decision.confidence,
                'new_contigs': ';'.join(decision.new_contig_ids),
                'new_lengths': ';'.join(map(str, decision.new_lengths))
            }
            rows.append(row)
        
        df = pd.DataFrame(rows)
        df.to_csv(output_path, sep='\t', index=False)
    
    def _write_results_json(self,
                          analyses: List[ChimeraAnalysis],
                          decisions: List[SplittingDecision],
                          output_path: Path):
        """Write comprehensive results to JSON file."""
        import json
        
        results = {
            'summary': {
                'total_analyses': len(analyses),
                'total_decisions': len(decisions),
                'contigs_split': len([d for d in decisions if d.action == 'split']),
                'contigs_preserved': len([d for d in decisions if d.action == 'preserve']),
                'contigs_flagged': len([d for d in decisions if d.action == 'flag'])
            },
            'analyses': [],
            'decisions': []
        }
        
        # Add analyses
        for analysis in analyses:
            analysis_dict = {
                'contig_id': analysis.candidate.contig_id,
                'breakpoint': analysis.candidate.breakpoint,
                'confidence_score': analysis.candidate.confidence_score,
                'evidence_types': analysis.candidate.evidence_types,
                'chimera_type': analysis.chimera_type,
                'classification_confidence': analysis.classification_confidence,
                'explanation': analysis.explanation,
                'recommendation': analysis.recommendation,
                'taxonomic_left': analysis.taxonomic_left,
                'taxonomic_right': analysis.taxonomic_right
            }
            results['analyses'].append(analysis_dict)
        
        # Add decisions
        for decision in decisions:
            decision_dict = {
                'contig_id': decision.contig_id,
                'original_length': decision.original_length,
                'breakpoint': decision.breakpoint,
                'action': decision.action,
                'reason': decision.reason,
                'new_contig_ids': decision.new_contig_ids,
                'new_lengths': decision.new_lengths,
                'confidence': decision.confidence
            }
            results['decisions'].append(decision_dict)
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
    
    def _generate_summary_stats(self,
                              original_contigs: Dict[str, str],
                              cleaned_contigs: Dict[str, str],
                              decisions: List[SplittingDecision]) -> Dict:
        """Generate summary statistics."""
        
        original_lengths = [len(seq) for seq in original_contigs.values()]
        cleaned_lengths = [len(seq) for seq in cleaned_contigs.values()]
        
        stats = {
            'original_contigs': len(original_contigs),
            'cleaned_contigs': len(cleaned_contigs),
            'contigs_split': len([d for d in decisions if d.action == 'split']),
            'contigs_preserved': len([d for d in decisions if d.action == 'preserve']),
            'contigs_flagged': len([d for d in decisions if d.action == 'flag']),
            'original_n50': calculate_n50(original_lengths),
            'cleaned_n50': calculate_n50(cleaned_lengths),
            'original_total_length': sum(original_lengths),
            'cleaned_total_length': sum(cleaned_lengths)
        }
        
        return stats
    
    def get_splitting_summary(self) -> pd.DataFrame:
        """Get summary of splitting decisions as DataFrame."""
        if not self.splitting_decisions:
            return pd.DataFrame()
        
        rows = []
        for decision in self.splitting_decisions:
            row = {
                'contig_id': decision.contig_id,
                'action': decision.action,
                'breakpoint': decision.breakpoint,
                'confidence': decision.confidence,
                'reason': decision.reason
            }
            rows.append(row)
        
        return pd.DataFrame(rows)