"""
Visualization module for read-pair based chimera detection results.

Creates interactive plots for manual inspection of detected anomalies.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from typing import List, Dict, Any, Optional
from pathlib import Path
import logging

from .analyzer import BreakpointCandidate, WindowMetrics
from .bamparser import BamParser
from .config import DetectorConfig


class ChimeraVisualizer:
    """Visualizer for chimera detection results."""
    
    def __init__(self, config: DetectorConfig):
        """Initialize visualizer with configuration."""
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def create_contig_plot(self, parser: BamParser, contig: str,
                          candidates: List[BreakpointCandidate],
                          window_metrics: List[WindowMetrics]) -> go.Figure:
        """
        Create comprehensive plot for a single contig.
        
        Args:
            parser: BAM parser instance
            contig: Contig name
            candidates: Breakpoint candidates for this contig
            window_metrics: Window metrics along the contig
            
        Returns:
            Plotly figure
        """
        # Create subplots
        fig = make_subplots(
            rows=4, cols=1,
            subplot_titles=["Proper Pair Ratio", "Insert Size Distribution", 
                           "Discordant Pair Ratio", "Coverage Uniformity"],
            vertical_spacing=0.08,
            shared_xaxes=True
        )
        
        # Extract data from window metrics
        positions = [(m.start + m.end) // 2 for m in window_metrics]
        proper_ratios = [m.proper_pair_ratio for m in window_metrics]
        insert_medians = [m.insert_size_median for m in window_metrics]
        discordant_ratios = [m.discordant_ratio for m in window_metrics]
        uniformity = [m.coverage_uniformity for m in window_metrics]
        
        # Plot 1: Proper pair ratio
        fig.add_trace(
            go.Scatter(x=positions, y=proper_ratios, mode='lines+markers',
                      name='Proper Pair Ratio', line=dict(color='blue')),
            row=1, col=1
        )
        
        # Plot 2: Insert size median
        fig.add_trace(
            go.Scatter(x=positions, y=insert_medians, mode='lines+markers',
                      name='Insert Size Median', line=dict(color='green')),
            row=2, col=1
        )
        
        # Plot 3: Discordant pair ratio
        fig.add_trace(
            go.Scatter(x=positions, y=discordant_ratios, mode='lines+markers',
                      name='Discordant Ratio', line=dict(color='red')),
            row=3, col=1
        )
        
        # Plot 4: Coverage uniformity
        fig.add_trace(
            go.Scatter(x=positions, y=uniformity, mode='lines+markers',
                      name='Coverage Uniformity', line=dict(color='purple')),
            row=4, col=1
        )
        
        # Add breakpoint markers
        for candidate in candidates:
            if candidate.contig == contig:
                # Add vertical lines for breakpoints
                for row in [1, 2, 3, 4]:
                    fig.add_vline(
                        x=candidate.position,
                        line_color='orange',
                        line_width=2,
                        line_dash='dash',
                        row=row, col=1
                    )
                
                # Add annotation on first subplot
                fig.add_annotation(
                    x=candidate.position,
                    y=max(proper_ratios) if proper_ratios else 1,
                    text=f"BP: {candidate.confidence:.2f}",
                    showarrow=True,
                    arrowhead=2,
                    bgcolor='orange',
                    bordercolor='orange',
                    font=dict(color='white'),
                    row=1, col=1
                )
        
        # Update layout
        fig.update_layout(
            title=f"Chimera Analysis: {contig}",
            height=800,
            showlegend=False
        )
        
        # Update x-axis for bottom subplot
        fig.update_xaxes(title_text="Position (bp)", row=4, col=1)
        
        # Update y-axes
        fig.update_yaxes(title_text="Ratio", row=1, col=1)
        fig.update_yaxes(title_text="Insert Size (bp)", row=2, col=1)
        fig.update_yaxes(title_text="Ratio", row=3, col=1)
        fig.update_yaxes(title_text="Uniformity", row=4, col=1)
        
        return fig
    
    def create_overview_plot(self, candidates: List[BreakpointCandidate]) -> go.Figure:
        """
        Create overview plot of all candidates.
        
        Args:
            candidates: All breakpoint candidates
            
        Returns:
            Plotly figure
        """
        if not candidates:
            # Empty plot
            fig = go.Figure()
            fig.add_annotation(
                text="No candidates found",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig
        
        # Group by contig
        by_contig = {}
        for candidate in candidates:
            if candidate.contig not in by_contig:
                by_contig[candidate.contig] = []
            by_contig[candidate.contig].append(candidate)
        
        # Create scatter plot
        fig = go.Figure()
        
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        for i, (contig, contig_candidates) in enumerate(by_contig.items()):
            positions = [c.position for c in contig_candidates]
            confidences = [c.confidence for c in contig_candidates]
            
            color = colors[i % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=positions,
                y=confidences,
                mode='markers',
                name=contig,
                marker=dict(
                    size=10,
                    color=color,
                    symbol='circle'
                ),
                hovertemplate=f"<b>{contig}</b><br>" +
                             "Position: %{x}<br>" +
                             "Confidence: %{y:.3f}<extra></extra>"
            ))
        
        fig.update_layout(
            title="Chimera Detection Overview",
            xaxis_title="Position (bp)",
            yaxis_title="Confidence Score",
            height=500
        )
        
        return fig
    
    def create_insert_size_plot(self, parser: BamParser, contig: str,
                               breakpoint_position: int) -> go.Figure:
        """
        Create detailed insert size distribution plot around a breakpoint.
        
        Args:
            parser: BAM parser instance
            contig: Contig name
            breakpoint_position: Position of breakpoint
            
        Returns:
            Plotly figure
        """
        window_size = 2000  # Window around breakpoint
        start = max(0, breakpoint_position - window_size)
        end = breakpoint_position + window_size
        
        # Get read pairs in region
        pairs = list(parser.get_read_pairs(contig, start, end))
        
        if not pairs:
            fig = go.Figure()
            fig.add_annotation(
                text="No read pairs found in region",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig
        
        # Separate into left and right of breakpoint
        left_pairs = [p for p in pairs if (p.read1_start + p.read1_end) // 2 < breakpoint_position]
        right_pairs = [p for p in pairs if (p.read1_start + p.read1_end) // 2 >= breakpoint_position]
        
        # Extract insert sizes
        left_inserts = [p.insert_size for p in left_pairs if 0 < p.insert_size < 10000]
        right_inserts = [p.insert_size for p in right_pairs if 0 < p.insert_size < 10000]
        
        # Create histograms
        fig = go.Figure()
        
        if left_inserts:
            fig.add_trace(go.Histogram(
                x=left_inserts,
                name=f"Left of breakpoint (n={len(left_inserts)})",
                opacity=0.7,
                nbinsx=50
            ))
        
        if right_inserts:
            fig.add_trace(go.Histogram(
                x=right_inserts,
                name=f"Right of breakpoint (n={len(right_inserts)})",
                opacity=0.7,
                nbinsx=50
            ))
        
        fig.update_layout(
            title=f"Insert Size Distribution: {contig}:{breakpoint_position}",
            xaxis_title="Insert Size (bp)",
            yaxis_title="Count",
            barmode='overlay',
            height=400
        )
        
        return fig
    
    def generate_html_report(self, candidates: List[BreakpointCandidate],
                           metadata: Dict[str, Any],
                           output_dir: str,
                           parser: Optional[BamParser] = None) -> str:
        """
        Generate comprehensive HTML report.
        
        Args:
            candidates: Breakpoint candidates
            metadata: Analysis metadata
            output_dir: Output directory
            parser: Optional BAM parser for detailed plots
            
        Returns:
            Path to HTML report
        """
        output_path = Path(output_dir)
        report_path = output_path / "chimera_report.html"
        
        # Create overview plot
        overview_fig = self.create_overview_plot(candidates)
        
        # Start HTML content
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Chimera Detection Report</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .container {{ max-width: 1200px; margin: 0 auto; }}
                .summary {{ background-color: #f0f0f0; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
                .candidate {{ border: 1px solid #ddd; padding: 10px; margin: 10px 0; border-radius: 5px; }}
                .high-conf {{ border-left: 5px solid #ff0000; }}
                .med-conf {{ border-left: 5px solid #ff9900; }}
                .low-conf {{ border-left: 5px solid #ffff00; }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>Chimera Detection Report</h1>
                
                <div class="summary">
                    <h2>Summary</h2>
                    <p><strong>Input BAM:</strong> {metadata.get('input_bam', 'N/A')}</p>
                    <p><strong>Analysis Time:</strong> {metadata.get('analysis_time_seconds', 'N/A')} seconds</p>
                    <p><strong>Contigs Analyzed:</strong> {metadata.get('contigs_analyzed', 'N/A')}</p>
                    <p><strong>Total Candidates:</strong> {len(candidates)}</p>
                    <p><strong>High Confidence (≥0.8):</strong> {len([c for c in candidates if c.confidence >= 0.8])}</p>
                    <p><strong>Medium Confidence (0.5-0.8):</strong> {len([c for c in candidates if 0.5 <= c.confidence < 0.8])}</p>
                    <p><strong>Low Confidence (<0.5):</strong> {len([c for c in candidates if c.confidence < 0.5])}</p>
                </div>
                
                <h2>Overview</h2>
                <div id="overview-plot"></div>
                
                <h2>Detailed Results</h2>
        """
        
        # Add candidate details
        for i, candidate in enumerate(sorted(candidates, key=lambda x: x.confidence, reverse=True)):
            confidence_class = "high-conf" if candidate.confidence >= 0.8 else "med-conf" if candidate.confidence >= 0.5 else "low-conf"
            
            anomaly_types = ", ".join(set(a.anomaly_type for a in candidate.supporting_anomalies))
            
            html_content += f"""
                <div class="candidate {confidence_class}">
                    <h3>Candidate {i+1}: {candidate.contig}:{candidate.position}</h3>
                    <p><strong>Confidence:</strong> {candidate.confidence:.3f}</p>
                    <p><strong>Anomaly Types:</strong> {anomaly_types}</p>
                    <p><strong>Supporting Evidence:</strong> {len(candidate.supporting_anomalies)} anomalies</p>
                    <p><strong>Left Proper Pair Ratio:</strong> {candidate.left_metrics.proper_pair_ratio:.3f}</p>
                    <p><strong>Right Proper Pair Ratio:</strong> {candidate.right_metrics.proper_pair_ratio:.3f}</p>
                </div>
            """
        
        # Close HTML
        html_content += """
                </div>
                
                <script>
                    // Add overview plot
                    var overview_data = """ + overview_fig.to_json() + """;
                    Plotly.newPlot('overview-plot', overview_data.data, overview_data.layout);
                </script>
            </body>
        </html>
        """
        
        # Write HTML file
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report generated: {report_path}")
        return str(report_path)