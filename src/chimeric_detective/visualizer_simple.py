"""
Simplified visualization module for GC-based chimera detection results.
"""

import logging
from typing import List, Dict, Optional
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import SeqIO

from .detector_simple import ChimeraCandidate
from .analyzer_simple import ChimeraClassification
from .utils import setup_logging
from .detector_simple import calculate_gc_content_simple


class SimpleChimeraVisualizer:
    """Simplified visualizer for chimera detection results."""
    
    def __init__(self,
                 output_dir: str = "chimeric_detective_output",
                 log_level: str = "INFO"):
        """
        Initialize SimpleChimeraVisualizer.
        
        Args:
            output_dir: Directory for output files
            log_level: Logging level
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = setup_logging(log_level)
    
    def generate_report(self,
                       assembly_file: str,
                       candidates: List[ChimeraCandidate],
                       classifications: List[ChimeraClassification],
                       output_prefix: str = "chimera_report") -> str:
        """
        Generate an HTML report with visualizations.
        
        Args:
            assembly_file: Path to assembly file
            candidates: List of chimera candidates
            classifications: List of classifications
            output_prefix: Prefix for output files
            
        Returns:
            Path to generated HTML report
        """
        self.logger.info("Generating visualization report")
        
        # Load sequences for visualization
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(assembly_file, "fasta")}
        
        # Create HTML report
        html_content = self._create_html_report(sequences, candidates, classifications)
        
        # Write report
        report_file = self.output_dir / f"{output_prefix}.html"
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"Report generated: {report_file}")
        return str(report_file)
    
    def _create_html_report(self,
                           sequences: Dict[str, str],
                           candidates: List[ChimeraCandidate],
                           classifications: List[ChimeraClassification]) -> str:
        """Create HTML report with embedded visualizations."""
        # Group candidates by contig
        candidates_by_contig = {}
        for candidate in candidates:
            if candidate.contig_id not in candidates_by_contig:
                candidates_by_contig[candidate.contig_id] = []
            candidates_by_contig[candidate.contig_id].append(candidate)
        
        # Create classification lookup
        classification_lookup = {
            (c.candidate.contig_id, c.candidate.breakpoint): c 
            for c in classifications
        }
        
        # Generate plots for each contig with candidates
        plots_html = []
        for contig_id, contig_candidates in candidates_by_contig.items():
            if contig_id in sequences:
                plot_html = self._create_contig_plot(
                    contig_id, sequences[contig_id], contig_candidates, classification_lookup
                )
                plots_html.append(plot_html)
        
        # Create summary statistics
        total_candidates = len(candidates)
        high_confidence = sum(1 for c in classifications if c.is_likely_chimera and c.confidence >= 0.7)
        
        # Build HTML
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Simplified Chimera Detection Report</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 20px;
                    background-color: #f5f5f5;
                }}
                .container {{
                    max-width: 1200px;
                    margin: 0 auto;
                    background-color: white;
                    padding: 20px;
                    border-radius: 5px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                }}
                h1, h2 {{
                    color: #333;
                }}
                .summary {{
                    background-color: #f0f0f0;
                    padding: 15px;
                    border-radius: 5px;
                    margin-bottom: 20px;
                }}
                .plot-container {{
                    margin-bottom: 30px;
                    border: 1px solid #ddd;
                    padding: 10px;
                    border-radius: 5px;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>Simplified Chimera Detection Report</h1>
                
                <div class="summary">
                    <h2>Summary</h2>
                    <p><strong>Total candidates detected:</strong> {total_candidates}</p>
                    <p><strong>High confidence chimeras:</strong> {high_confidence}</p>
                    <p><strong>Detection method:</strong> GC content analysis only</p>
                </div>
                
                <h2>Detected Chimeras</h2>
                {"".join(plots_html)}
            </div>
        </body>
        </html>
        """
        
        return html
    
    def _create_contig_plot(self,
                           contig_id: str,
                           sequence: str,
                           candidates: List[ChimeraCandidate],
                           classification_lookup: Dict) -> str:
        """Create a plot for a single contig showing GC content and breakpoints."""
        # Calculate GC content profile
        window_size = 100
        step_size = 50
        positions = []
        gc_values = []
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window_seq = sequence[i:i + window_size]
            gc = calculate_gc_content_simple(window_seq)
            positions.append(i + window_size // 2)
            gc_values.append(gc)
        
        # Create plot
        fig = go.Figure()
        
        # Add GC content trace
        fig.add_trace(go.Scatter(
            x=positions,
            y=gc_values,
            mode='lines',
            name='GC Content',
            line=dict(color='blue', width=2)
        ))
        
        # Add breakpoint markers
        for candidate in candidates:
            classification = classification_lookup.get((contig_id, candidate.breakpoint))
            if classification and classification.is_likely_chimera:
                color = 'red' if classification.confidence >= 0.7 else 'orange'
            else:
                color = 'gray'
            
            fig.add_vline(
                x=candidate.breakpoint,
                line_color=color,
                line_width=2,
                line_dash="dash"
            )
            
            # Add annotation
            fig.add_annotation(
                x=candidate.breakpoint,
                y=max(gc_values),
                text=f"GC diff: {candidate.gc_difference:.2f}",
                showarrow=True,
                arrowhead=2,
                bgcolor=color,
                bordercolor=color,
                font=dict(color="white")
            )
        
        # Update layout
        fig.update_layout(
            title=f"Contig: {contig_id} (Length: {len(sequence):,} bp)",
            xaxis_title="Position (bp)",
            yaxis_title="GC Content",
            height=400,
            showlegend=True,
            hovermode='x unified'
        )
        
        # Convert to HTML div
        plot_html = f"""
        <div class="plot-container">
            <h3>{contig_id}</h3>
            {fig.to_html(div_id=f"plot_{contig_id}", include_plotlyjs=False)}
        </div>
        """
        
        return plot_html