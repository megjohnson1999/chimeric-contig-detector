"""
Output formatting module for chimera detection results.

Supports multiple output formats and customizable formatting options.
"""

import json
import csv
from typing import List, Dict, Any, Optional, TextIO
from pathlib import Path
import logging
from datetime import datetime
import gzip
from abc import ABC, abstractmethod

from .analyzer import BreakpointCandidate, AnomalyRegion, WindowMetrics
from .config import OutputConfig


class OutputFormatter(ABC):
    """Abstract base class for output formatters."""
    
    @abstractmethod
    def format(self, results: List[BreakpointCandidate], 
              metadata: Dict[str, Any]) -> str:
        """Format results to string."""
        pass
    
    @abstractmethod
    def write(self, results: List[BreakpointCandidate], 
             metadata: Dict[str, Any], 
             output_file: TextIO) -> None:
        """Write results to file."""
        pass


class JSONFormatter(OutputFormatter):
    """JSON output formatter."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def format(self, results: List[BreakpointCandidate], 
              metadata: Dict[str, Any]) -> str:
        """Format results as JSON string."""
        output = {
            "metadata": metadata,
            "summary": {
                "total_candidates": len(results),
                "high_confidence": len([r for r in results if r.confidence >= 0.8]),
                "medium_confidence": len([r for r in results if 0.5 <= r.confidence < 0.8]),
                "low_confidence": len([r for r in results if r.confidence < 0.5])
            },
            "breakpoints": []
        }
        
        for candidate in results:
            breakpoint_data = {
                "contig": candidate.contig,
                "position": candidate.position,
                "confidence": round(candidate.confidence, self.config.decimal_places),
                "anomaly_count": len(candidate.supporting_anomalies),
                "anomaly_types": list(set(a.anomaly_type for a in candidate.supporting_anomalies))
            }
            
            if self.config.include_debug_info:
                breakpoint_data["debug"] = {
                    "left_metrics": self._metrics_to_dict(candidate.left_metrics),
                    "right_metrics": self._metrics_to_dict(candidate.right_metrics),
                    "anomalies": [self._anomaly_to_dict(a) for a in candidate.supporting_anomalies]
                }
            
            output["breakpoints"].append(breakpoint_data)
        
        # Sort by confidence
        output["breakpoints"].sort(key=lambda x: x["confidence"], reverse=True)
        
        return json.dumps(output, indent=2)
    
    def write(self, results: List[BreakpointCandidate], 
             metadata: Dict[str, Any], 
             output_file: TextIO) -> None:
        """Write JSON results to file."""
        json_str = self.format(results, metadata)
        output_file.write(json_str)
    
    def _metrics_to_dict(self, metrics: WindowMetrics) -> Dict[str, Any]:
        """Convert WindowMetrics to dictionary."""
        return {
            "start": metrics.start,
            "end": metrics.end,
            "total_pairs": metrics.total_pairs,
            "proper_pairs": metrics.proper_pairs,
            "proper_pair_ratio": round(metrics.proper_pair_ratio, self.config.decimal_places),
            "insert_size_median": round(metrics.insert_size_median, 2),
            "insert_size_mad": round(metrics.insert_size_mad, 2),
            "discordant_ratio": round(metrics.discordant_ratio, self.config.decimal_places)
        }
    
    def _anomaly_to_dict(self, anomaly: AnomalyRegion) -> Dict[str, Any]:
        """Convert AnomalyRegion to dictionary."""
        return {
            "start": anomaly.start,
            "end": anomaly.end,
            "type": anomaly.anomaly_type,
            "score": round(anomaly.score, self.config.decimal_places)
        }


class TSVFormatter(OutputFormatter):
    """TSV/CSV output formatter."""
    
    def __init__(self, config: OutputConfig, delimiter: str = "\t"):
        self.config = config
        self.delimiter = delimiter
        self.logger = logging.getLogger(__name__)
    
    def format(self, results: List[BreakpointCandidate], 
              metadata: Dict[str, Any]) -> str:
        """Format results as TSV string."""
        lines = []
        
        # Header
        headers = ["contig", "position", "confidence", "anomaly_types", 
                  "left_proper_ratio", "right_proper_ratio",
                  "left_insert_median", "right_insert_median"]
        
        if self.config.include_debug_info:
            headers.extend(["anomaly_count", "anomaly_scores"])
        
        lines.append(self.delimiter.join(headers))
        
        # Data rows
        for candidate in results:
            anomaly_types = ",".join(set(a.anomaly_type for a in candidate.supporting_anomalies))
            
            row = [
                candidate.contig,
                str(candidate.position),
                f"{candidate.confidence:.{self.config.decimal_places}f}",
                anomaly_types,
                f"{candidate.left_metrics.proper_pair_ratio:.{self.config.decimal_places}f}",
                f"{candidate.right_metrics.proper_pair_ratio:.{self.config.decimal_places}f}",
                f"{candidate.left_metrics.insert_size_median:.2f}",
                f"{candidate.right_metrics.insert_size_median:.2f}"
            ]
            
            if self.config.include_debug_info:
                anomaly_scores = ",".join(f"{a.score:.2f}" for a in candidate.supporting_anomalies)
                row.extend([
                    str(len(candidate.supporting_anomalies)),
                    anomaly_scores
                ])
            
            lines.append(self.delimiter.join(row))
        
        return "\n".join(lines)
    
    def write(self, results: List[BreakpointCandidate], 
             metadata: Dict[str, Any], 
             output_file: TextIO) -> None:
        """Write TSV results to file."""
        tsv_str = self.format(results, metadata)
        output_file.write(tsv_str)


class BEDFormatter(OutputFormatter):
    """BED format output for genome browsers."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def format(self, results: List[BreakpointCandidate], 
              metadata: Dict[str, Any]) -> str:
        """Format results as BED."""
        lines = []
        
        # BED header
        lines.append(f"track name=\"Chimera_Breakpoints\" description=\"{metadata.get('description', 'Chimeric breakpoints')}\"")
        
        for i, candidate in enumerate(results):
            # BED format: chrom start end name score
            # Use a window around the breakpoint
            window = 100
            start = max(0, candidate.position - window)
            end = candidate.position + window
            
            # Convert confidence to BED score (0-1000)
            score = int(candidate.confidence * 1000)
            
            # Color based on confidence
            if candidate.confidence >= 0.8:
                color = "255,0,0"  # Red for high confidence
            elif candidate.confidence >= 0.5:
                color = "255,165,0"  # Orange for medium
            else:
                color = "255,255,0"  # Yellow for low
            
            line = f"{candidate.contig}\t{start}\t{end}\tBreakpoint_{i+1}\t{score}\t+\t{start}\t{end}\t{color}"
            lines.append(line)
        
        return "\n".join(lines)
    
    def write(self, results: List[BreakpointCandidate], 
             metadata: Dict[str, Any], 
             output_file: TextIO) -> None:
        """Write BED results to file."""
        bed_str = self.format(results, metadata)
        output_file.write(bed_str)


class VCFFormatter(OutputFormatter):
    """VCF format output for variant representation."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def format(self, results: List[BreakpointCandidate], 
              metadata: Dict[str, Any]) -> str:
        """Format results as VCF."""
        lines = []
        
        # VCF header
        lines.append("##fileformat=VCFv4.3")
        lines.append(f"##fileDate={datetime.now().strftime('%Y%m%d')}")
        lines.append(f"##source=ChimericDetective_ReadPair_v2.0")
        lines.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
        lines.append("##INFO=<ID=CONF,Number=1,Type=Float,Description=\"Confidence score\">")
        lines.append("##INFO=<ID=ANOM,Number=.,Type=String,Description=\"Supporting anomaly types\">")
        lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        
        for i, candidate in enumerate(results):
            # VCF fields
            chrom = candidate.contig
            pos = candidate.position
            id_field = f"CHIMERA_{i+1}"
            ref = "N"  # Unknown reference
            alt = "<BND>"  # Breakend notation
            qual = int(candidate.confidence * 100)  # Scale to 0-100
            filter_field = "PASS" if candidate.confidence >= 0.5 else "LowConf"
            
            # INFO field
            anomaly_types = ",".join(set(a.anomaly_type for a in candidate.supporting_anomalies))
            info_parts = [
                "SVTYPE=BND",
                f"CONF={candidate.confidence:.{self.config.decimal_places}f}",
                f"ANOM={anomaly_types}"
            ]
            info = ";".join(info_parts)
            
            line = f"{chrom}\t{pos}\t{id_field}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}"
            lines.append(line)
        
        return "\n".join(lines)
    
    def write(self, results: List[BreakpointCandidate], 
             metadata: Dict[str, Any], 
             output_file: TextIO) -> None:
        """Write VCF results to file."""
        vcf_str = self.format(results, metadata)
        output_file.write(vcf_str)


class OutputManager:
    """Manager for handling multiple output formats."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Initialize formatters
        self.formatters = {
            "json": JSONFormatter(config),
            "tsv": TSVFormatter(config),
            "csv": TSVFormatter(config, delimiter=","),
            "bed": BEDFormatter(config),
            "vcf": VCFFormatter(config)
        }
    
    def write_results(self, results: List[BreakpointCandidate],
                     metadata: Dict[str, Any],
                     output_dir: str) -> Dict[str, str]:
        """
        Write results in all configured formats.
        
        Returns:
            Dictionary mapping format to output file path
        """
        output_files = {}
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for format_name in self.config.formats:
            if format_name not in self.formatters:
                self.logger.warning(f"Unknown format: {format_name}")
                continue
            
            formatter = self.formatters[format_name]
            filename = f"{self.config.output_prefix}.{format_name}"
            
            if self.config.compress_output and format_name in ["json", "tsv", "csv"]:
                filename += ".gz"
                filepath = output_path / filename
                
                with gzip.open(filepath, 'wt') as f:
                    formatter.write(results, metadata, f)
            else:
                filepath = output_path / filename
                
                with open(filepath, 'w') as f:
                    formatter.write(results, metadata, f)
            
            output_files[format_name] = str(filepath)
            self.logger.info(f"Wrote {format_name} output to {filepath}")
        
        return output_files
    
    def create_summary_report(self, results: List[BreakpointCandidate],
                            metadata: Dict[str, Any],
                            output_dir: str) -> str:
        """Create a human-readable summary report."""
        output_path = Path(output_dir)
        report_path = output_path / f"{self.config.output_prefix}_summary.txt"
        
        with open(report_path, 'w') as f:
            f.write("CHIMERA DETECTION SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Metadata
            f.write("Analysis Details:\n")
            for key, value in metadata.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
            
            # Summary statistics
            f.write("Results Summary:\n")
            f.write(f"  Total candidates: {len(results)}\n")
            f.write(f"  High confidence (≥0.8): {len([r for r in results if r.confidence >= 0.8])}\n")
            f.write(f"  Medium confidence (0.5-0.8): {len([r for r in results if 0.5 <= r.confidence < 0.8])}\n")
            f.write(f"  Low confidence (<0.5): {len([r for r in results if r.confidence < 0.5])}\n")
            f.write("\n")
            
            # Top candidates
            f.write("Top Candidates:\n")
            f.write("-" * 50 + "\n")
            
            for i, candidate in enumerate(sorted(results, key=lambda x: x.confidence, reverse=True)[:10]):
                f.write(f"\n{i+1}. {candidate.contig}:{candidate.position}\n")
                f.write(f"   Confidence: {candidate.confidence:.3f}\n")
                f.write(f"   Anomalies: {', '.join(set(a.anomaly_type for a in candidate.supporting_anomalies))}\n")
                f.write(f"   Left/Right proper pair ratio: {candidate.left_metrics.proper_pair_ratio:.3f} / {candidate.right_metrics.proper_pair_ratio:.3f}\n")
        
        return str(report_path)