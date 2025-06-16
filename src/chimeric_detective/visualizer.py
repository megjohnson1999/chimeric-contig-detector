"""
Visualization module for creating interactive HTML reports and static plots.
"""

import logging
import os
import json
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
from jinja2 import Template

from .analyzer import ChimeraAnalysis
from .resolver import SplittingDecision
from .utils import setup_logging


class ChimeraVisualizer:
    """Visualizer for creating comprehensive chimera analysis reports."""
    
    def __init__(self,
                 output_format: str = "html",
                 include_static: bool = True,
                 color_scheme: str = "viridis",
                 log_level: str = "INFO"):
        """
        Initialize ChimeraVisualizer.
        
        Args:
            output_format: Output format ('html', 'png', 'both')
            include_static: Whether to include static PNG exports
            color_scheme: Color scheme for plots
            log_level: Logging level
        """
        self.output_format = output_format
        self.include_static = include_static
        self.color_scheme = color_scheme
        
        self.logger = setup_logging(log_level)
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        sns.set_palette(color_scheme)
    
    def create_report(self,
                     analyses: List[ChimeraAnalysis],
                     decisions: List[SplittingDecision],
                     output_dir: str,
                     assembly_file: str) -> str:
        """
        Create comprehensive HTML report with interactive visualizations.
        
        Args:
            analyses: List of ChimeraAnalysis objects
            decisions: List of SplittingDecision objects
            output_dir: Output directory
            assembly_file: Path to original assembly file
            
        Returns:
            Path to generated HTML report
        """
        self.logger.info("Creating interactive HTML report")
        
        output_dir = Path(output_dir)
        figures_dir = output_dir / "figures"
        figures_dir.mkdir(exist_ok=True)
        
        # Generate all visualizations
        try:
            plots = self._generate_all_plots(analyses, decisions, figures_dir)
            self.logger.info(f"Generated {len(plots)} plot types")
        except Exception as e:
            self.logger.error(f"Failed to generate plots: {e}")
            plots = {}
        
        # Create summary statistics
        summary_stats = self._calculate_summary_statistics(analyses, decisions)
        
        # Generate HTML report
        report_path = self._generate_html_report(
            analyses, decisions, plots, summary_stats, output_dir
        )
        
        self.logger.info(f"Report generated: {report_path}")
        return str(report_path)
    
    def _generate_all_plots(self,
                          analyses: List[ChimeraAnalysis],
                          decisions: List[SplittingDecision],
                          figures_dir: Path) -> Dict[str, str]:
        """Generate all plots and return their paths."""
        plots = {}
        
        # Summary plots
        try:
            self.logger.info("Generating chimera types plot...")
            plots['chimera_types'] = self._plot_chimera_types(analyses, figures_dir)
        except Exception as e:
            self.logger.warning(f"Failed to generate chimera types plot: {e}")
            
        try:
            self.logger.info("Generating confidence distribution plot...")
            plots['confidence_distribution'] = self._plot_confidence_distribution(analyses, figures_dir)
        except Exception as e:
            self.logger.warning(f"Failed to generate confidence distribution plot: {e}")
            
        try:
            self.logger.info("Generating decision summary plot...")
            plots['decision_summary'] = self._plot_decision_summary(decisions, figures_dir)
        except Exception as e:
            self.logger.warning(f"Failed to generate decision summary plot: {e}")
        
        # Individual contig plots
        try:
            self.logger.info("Generating individual contig plots...")
            plots['individual_contigs'] = self._plot_individual_contigs(analyses, figures_dir)
        except Exception as e:
            self.logger.warning(f"Failed to generate individual contig plots: {e}")
        
        # Coverage and evidence plots
        try:
            self.logger.info("Generating evidence overview plot...")
            plots['evidence_overview'] = self._plot_evidence_overview(analyses, figures_dir)
        except Exception as e:
            self.logger.warning(f"Failed to generate evidence overview plot: {e}")
        
        return plots
    
    def _plot_chimera_types(self, analyses: List[ChimeraAnalysis], figures_dir: Path) -> str:
        """Plot distribution of chimera types."""
        
        # Count chimera types
        type_counts = {}
        for analysis in analyses:
            chimera_type = analysis.chimera_type
            type_counts[chimera_type] = type_counts.get(chimera_type, 0) + 1
        
        # Create interactive plot
        fig = go.Figure(data=[
            go.Pie(
                labels=list(type_counts.keys()),
                values=list(type_counts.values()),
                hole=0.3,
                textinfo='label+percent+value',
                textfont_size=12
            )
        ])
        
        fig.update_layout(
            title="Distribution of Chimera Types",
            font=dict(size=14),
            showlegend=True,
            height=500
        )
        
        # Save interactive version
        html_path = figures_dir / "chimera_types.html"
        fig.write_html(str(html_path))
        
        # Save static version if requested
        if self.include_static:
            try:
                png_path = figures_dir / "chimera_types.png"
                fig.write_image(str(png_path), width=800, height=500)
            except Exception as e:
                self.logger.warning(f"Could not save static image: {e}. Install kaleido for PNG export.")
        
        return str(html_path)
    
    def _plot_confidence_distribution(self, analyses: List[ChimeraAnalysis], figures_dir: Path) -> str:
        """Plot confidence score distributions."""
        
        # Prepare data
        data = []
        for analysis in analyses:
            data.append({
                'contig_id': analysis.candidate.contig_id,
                'detection_confidence': analysis.candidate.confidence_score,
                'classification_confidence': analysis.classification_confidence,
                'chimera_type': analysis.chimera_type
            })
        
        df = pd.DataFrame(data)
        
        # Create subplot with two histograms
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('Detection Confidence', 'Classification Confidence'),
            horizontal_spacing=0.1
        )
        
        # Detection confidence histogram
        fig.add_trace(
            go.Histogram(
                x=df['detection_confidence'],
                nbinsx=20,
                name='Detection',
                marker_color='lightblue',
                opacity=0.7
            ),
            row=1, col=1
        )
        
        # Classification confidence histogram
        fig.add_trace(
            go.Histogram(
                x=df['classification_confidence'],
                nbinsx=20,
                name='Classification',
                marker_color='lightcoral',
                opacity=0.7
            ),
            row=1, col=2
        )
        
        fig.update_layout(
            title="Confidence Score Distributions",
            height=400,
            showlegend=False
        )
        
        fig.update_xaxes(title_text="Confidence Score", row=1, col=1)
        fig.update_xaxes(title_text="Confidence Score", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=1, col=1)
        fig.update_yaxes(title_text="Count", row=1, col=2)
        
        # Save files
        html_path = figures_dir / "confidence_distribution.html"
        fig.write_html(str(html_path))
        
        if self.include_static:
            try:
                png_path = figures_dir / "confidence_distribution.png"
                fig.write_image(str(png_path), width=1000, height=400)
            except Exception as e:
                self.logger.warning(f"Could not save static image: {e}. Install kaleido for PNG export.")
        
        return str(html_path)
    
    def _plot_decision_summary(self, decisions: List[SplittingDecision], figures_dir: Path) -> str:
        """Plot summary of splitting decisions."""
        
        # Count decision types
        action_counts = {}
        for decision in decisions:
            action = decision.action
            action_counts[action] = action_counts.get(action, 0) + 1
        
        # Create bar plot
        fig = go.Figure(data=[
            go.Bar(
                x=list(action_counts.keys()),
                y=list(action_counts.values()),
                marker_color=['green' if x == 'preserve' else 'red' if x == 'split' else 'orange' 
                             for x in action_counts.keys()],
                text=list(action_counts.values()),
                textposition='auto'
            )
        ])
        
        fig.update_layout(
            title="Splitting Decisions Summary",
            xaxis_title="Decision Type",
            yaxis_title="Number of Contigs",
            height=400
        )
        
        # Save files
        html_path = figures_dir / "decision_summary.html"
        fig.write_html(str(html_path))
        
        if self.include_static:
            try:
                png_path = figures_dir / "decision_summary.png"
                fig.write_image(str(png_path), width=600, height=400)
            except Exception as e:
                self.logger.warning(f"Could not save static image: {e}. Install kaleido for PNG export.")
        
        return str(html_path)
    
    def _plot_individual_contigs(self, analyses: List[ChimeraAnalysis], figures_dir: Path) -> List[str]:
        """Create detailed plots for individual chimeric contigs."""
        
        individual_plots = []
        
        for analysis in analyses:
            candidate = analysis.candidate
            contig_id = candidate.contig_id
            
            # Create subplot for this contig
            fig = make_subplots(
                rows=3, cols=1,
                subplot_titles=(
                    f'Coverage Profile - {contig_id}',
                    f'GC Content - {contig_id}',
                    f'Evidence Summary - {contig_id}'
                ),
                vertical_spacing=0.1,
                row_heights=[0.4, 0.3, 0.3]
            )
            
            # Mock coverage data (would be real data in production)
            positions = np.arange(0, 5000, 50)
            coverage = np.random.normal(20, 5, len(positions))
            
            # Add breakpoint effect
            breakpoint_idx = candidate.breakpoint // 50
            if breakpoint_idx < len(coverage):
                left_cov = candidate.coverage_left
                right_cov = candidate.coverage_right
                coverage[:breakpoint_idx] = np.random.normal(left_cov, left_cov*0.2, breakpoint_idx)
                coverage[breakpoint_idx:] = np.random.normal(right_cov, right_cov*0.2, 
                                                           len(coverage) - breakpoint_idx)
            
            # Coverage plot
            fig.add_trace(
                go.Scatter(
                    x=positions,
                    y=coverage,
                    mode='lines',
                    name='Coverage',
                    line=dict(color='blue', width=2)
                ),
                row=1, col=1
            )
            
            # Add breakpoint line
            fig.add_vline(
                x=candidate.breakpoint,
                line_dash="dash",
                line_color="red",
                annotation_text="Breakpoint",
                row=1, col=1
            )
            
            # GC content plot (mock data)
            gc_positions = np.arange(0, 5000, 500)
            gc_content = np.random.normal(0.5, 0.1, len(gc_positions))
            
            fig.add_trace(
                go.Scatter(
                    x=gc_positions,
                    y=gc_content,
                    mode='lines+markers',
                    name='GC Content',
                    line=dict(color='green', width=2)
                ),
                row=2, col=1
            )
            
            # Evidence summary (bar plot)
            evidence_types = candidate.evidence_types
            evidence_counts = [1] * len(evidence_types)  # Each evidence type present once
            
            fig.add_trace(
                go.Bar(
                    x=evidence_types,
                    y=evidence_counts,
                    name='Evidence',
                    marker_color='orange'
                ),
                row=3, col=1
            )
            
            # Update layout
            fig.update_layout(
                height=800,
                title=f"Detailed Analysis: {contig_id}",
                showlegend=False
            )
            
            # Update axes
            fig.update_xaxes(title_text="Position (bp)", row=1, col=1)
            fig.update_yaxes(title_text="Coverage", row=1, col=1)
            fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
            fig.update_yaxes(title_text="GC Content", row=2, col=1)
            fig.update_xaxes(title_text="Evidence Type", row=3, col=1)
            fig.update_yaxes(title_text="Present", row=3, col=1)
            
            # Save individual contig plot
            html_path = figures_dir / f"{contig_id}_detailed.html"
            fig.write_html(str(html_path))
            individual_plots.append(str(html_path))
            
            if self.include_static:
                try:
                    png_path = figures_dir / f"{contig_id}_detailed.png"
                    fig.write_image(str(png_path), width=1000, height=800)
                except Exception as e:
                    self.logger.warning(f"Could not save static image: {e}. Install kaleido for PNG export.")
        
        return individual_plots
    
    def _plot_evidence_overview(self, analyses: List[ChimeraAnalysis], figures_dir: Path) -> str:
        """Create overview plot of evidence types."""
        
        # Count evidence types across all analyses
        evidence_counts = {}
        for analysis in analyses:
            for evidence_type in analysis.candidate.evidence_types:
                evidence_counts[evidence_type] = evidence_counts.get(evidence_type, 0) + 1
        
        # Create horizontal bar plot
        fig = go.Figure(data=[
            go.Bar(
                y=list(evidence_counts.keys()),
                x=list(evidence_counts.values()),
                orientation='h',
                marker_color='lightblue',
                text=list(evidence_counts.values()),
                textposition='auto'
            )
        ])
        
        fig.update_layout(
            title="Evidence Types Overview",
            xaxis_title="Number of Occurrences",
            yaxis_title="Evidence Type",
            height=400
        )
        
        # Save files
        html_path = figures_dir / "evidence_overview.html"
        fig.write_html(str(html_path))
        
        if self.include_static:
            try:
                png_path = figures_dir / "evidence_overview.png"
                fig.write_image(str(png_path), width=800, height=400)
            except Exception as e:
                self.logger.warning(f"Could not save static image: {e}. Install kaleido for PNG export.")
        
        return str(html_path)
    
    def _calculate_summary_statistics(self,
                                    analyses: List[ChimeraAnalysis],
                                    decisions: List[SplittingDecision]) -> Dict:
        """Calculate summary statistics for the report."""
        
        stats = {
            'total_contigs_analyzed': len(analyses),
            'chimera_types': {},
            'decision_types': {},
            'mean_detection_confidence': 0.0,
            'mean_classification_confidence': 0.0,
            'high_confidence_analyses': 0,
            'evidence_type_counts': {}
        }
        
        # Chimera type distribution
        for analysis in analyses:
            chimera_type = analysis.chimera_type
            stats['chimera_types'][chimera_type] = stats['chimera_types'].get(chimera_type, 0) + 1
        
        # Decision type distribution
        for decision in decisions:
            action = decision.action
            stats['decision_types'][action] = stats['decision_types'].get(action, 0) + 1
        
        # Confidence statistics
        if analyses:
            detection_confidences = [a.candidate.confidence_score for a in analyses]
            classification_confidences = [a.classification_confidence for a in analyses]
            
            stats['mean_detection_confidence'] = np.mean(detection_confidences)
            stats['mean_classification_confidence'] = np.mean(classification_confidences)
            stats['high_confidence_analyses'] = len([c for c in classification_confidences if c > 0.8])
        
        # Evidence type counts
        all_evidence_types = []
        for analysis in analyses:
            all_evidence_types.extend(analysis.candidate.evidence_types)
        
        for evidence_type in set(all_evidence_types):
            stats['evidence_type_counts'][evidence_type] = all_evidence_types.count(evidence_type)
        
        return stats
    
    def _generate_html_report(self,
                            analyses: List[ChimeraAnalysis],
                            decisions: List[SplittingDecision],
                            plots: Dict[str, str],
                            summary_stats: Dict,
                            output_dir: Path) -> str:
        """Generate the main HTML report."""
        
        # HTML template
        template_str = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Chimeric Detective Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }
        .header {
            text-align: center;
            color: #333;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        .stat-card {
            background-color: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid #4CAF50;
        }
        .stat-number {
            font-size: 2em;
            font-weight: bold;
            color: #4CAF50;
        }
        .stat-label {
            color: #666;
            margin-top: 5px;
        }
        .plot-section {
            margin-bottom: 40px;
        }
        .plot-container {
            width: 100%;
            height: 500px;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
        }
        .details-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }
        .details-table th, .details-table td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }
        .details-table th {
            background-color: #f2f2f2;
            font-weight: bold;
        }
        .details-table tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        .confidence-high { color: #28a745; font-weight: bold; }
        .confidence-medium { color: #ffc107; font-weight: bold; }
        .confidence-low { color: #dc3545; font-weight: bold; }
        .action-split { color: #dc3545; font-weight: bold; }
        .action-preserve { color: #28a745; font-weight: bold; }
        .action-flag { color: #ffc107; font-weight: bold; }
        .explanation-box {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
            border-left: 4px solid #4CAF50;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🔬 Chimeric Detective Report</h1>
            <p>Comprehensive Analysis of Chimeric Contigs in Viral Metagenomic Assembly</p>
        </div>
        
        <div class="summary-stats">
            <div class="stat-card">
                <div class="stat-number">{{ summary_stats.total_contigs_analyzed }}</div>
                <div class="stat-label">Contigs Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{{ summary_stats.decision_types.get('split', 0) }}</div>
                <div class="stat-label">Contigs Split</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{{ summary_stats.decision_types.get('preserve', 0) }}</div>
                <div class="stat-label">Contigs Preserved</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{{ summary_stats.mean_classification_confidence|round(2) }}</div>
                <div class="stat-label">Mean Confidence</div>
            </div>
        </div>
        
        <div class="plot-section">
            <h2>📊 Summary Visualizations</h2>
            
            <h3>Chimera Type Distribution</h3>
            <div class="plot-container">
                {{ chimera_types_plot | safe }}
            </div>
            
            <h3>Confidence Score Distributions</h3>
            <div class="plot-container">
                {{ confidence_distribution_plot | safe }}
            </div>
            
            <h3>Decision Summary</h3>
            <div class="plot-container">
                {{ decision_summary_plot | safe }}
            </div>
            
            <h3>Evidence Types Overview</h3>
            <div class="plot-container">
                {{ evidence_overview_plot | safe }}
            </div>
        </div>
        
        <div class="plot-section">
            <h2>🔍 Detailed Analysis Results</h2>
            
            <table class="details-table">
                <thead>
                    <tr>
                        <th>Contig ID</th>
                        <th>Chimera Type</th>
                        <th>Confidence</th>
                        <th>Decision</th>
                        <th>Breakpoint</th>
                        <th>Evidence Types</th>
                        <th>Explanation</th>
                    </tr>
                </thead>
                <tbody>
                    {% for analysis in analyses %}
                    <tr>
                        <td>{{ analysis.candidate.contig_id }}</td>
                        <td>{{ analysis.chimera_type }}</td>
                        <td class="{% if analysis.classification_confidence > 0.8 %}confidence-high{% elif analysis.classification_confidence > 0.5 %}confidence-medium{% else %}confidence-low{% endif %}">
                            {{ analysis.classification_confidence|round(3) }}
                        </td>
                        <td class="action-{{ decisions_dict[analysis.candidate.contig_id].action }}">
                            {{ decisions_dict[analysis.candidate.contig_id].action.upper() }}
                        </td>
                        <td>{{ analysis.candidate.breakpoint }}</td>
                        <td>{{ analysis.candidate.evidence_types|join(", ") }}</td>
                        <td>
                            <div class="explanation-box">
                                {{ analysis.explanation }}
                            </div>
                        </td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
        
        <div class="plot-section">
            <h2>📈 Individual Contig Details</h2>
            <p>Click on the links below to view detailed analysis for each chimeric contig:</p>
            <ul>
                {% for analysis in analyses %}
                <li>
                    <a href="figures/{{ analysis.candidate.contig_id }}_detailed.html" target="_blank">
                        {{ analysis.candidate.contig_id }} - {{ analysis.chimera_type }}
                    </a>
                </li>
                {% endfor %}
            </ul>
        </div>
        
        <div class="plot-section">
            <h2>ℹ️ Methodology & Interpretation</h2>
            <div class="explanation-box">
                <h3>Detection Methods</h3>
                <p>Chimeric contigs are detected using multiple complementary approaches:</p>
                <ul>
                    <li><strong>Coverage Discontinuities:</strong> Sharp changes in read coverage depth</li>
                    <li><strong>Sequence Composition:</strong> Changes in GC content and k-mer frequencies</li>
                    <li><strong>Taxonomic Classification:</strong> Transitions between different viral/host lineages</li>
                    <li><strong>Read Pair Orientation:</strong> Inconsistent paired-end read orientations</li>
                </ul>
                
                <h3>Classification Categories</h3>
                <ul>
                    <li><strong>Technical Artifacts:</strong> Assembly errors, typically split</li>
                    <li><strong>PCR Chimeras:</strong> Amplification artifacts, typically split</li>
                    <li><strong>Biological Recombination:</strong> Genuine recombination events, preserved</li>
                    <li><strong>Provirus Integration:</strong> Virus-host integration sites, flagged</li>
                </ul>
                
                <h3>Confidence Scores</h3>
                <p>Confidence scores range from 0-1, with higher scores indicating stronger evidence for the classification. 
                Scores above 0.8 are considered high confidence, 0.5-0.8 medium confidence, and below 0.5 low confidence.</p>
            </div>
        </div>
        
        <footer style="text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666;">
            <p>Generated by Chimeric Detective v1.0.0 | 
            <a href="https://github.com/yourusername/chimeric-detective">Documentation</a></p>
        </footer>
    </div>
</body>
</html>
        """
        
        # Create decisions lookup for template
        decisions_dict = {d.contig_id: d for d in decisions}
        
        # Read plot HTML content for direct embedding
        plot_content = {}
        figures_dir = output_dir / "figures"
        
        # Read each plot's HTML content
        plot_files = {
            'chimera_types_plot': 'chimera_types.html',
            'confidence_distribution_plot': 'confidence_distribution.html', 
            'decision_summary_plot': 'decision_summary.html',
            'evidence_overview_plot': 'evidence_overview.html'
        }
        
        for template_var, filename in plot_files.items():
            plot_path = figures_dir / filename
            if plot_path.exists():
                try:
                    with open(plot_path, 'r', encoding='utf-8') as f:
                        # Read the plot HTML and extract just the plot div content
                        plot_html = f.read()
                        # Extract the plotly div content (between <body> tags)
                        import re
                        body_match = re.search(r'<body[^>]*>(.*?)</body>', plot_html, re.DOTALL)
                        if body_match:
                            plot_content[template_var] = body_match.group(1)
                        else:
                            plot_content[template_var] = plot_html
                except Exception as e:
                    self.logger.warning(f"Could not read plot file {filename}: {e}")
                    plot_content[template_var] = f"<p>Plot could not be loaded: {filename}</p>"
            else:
                plot_content[template_var] = f"<p>Plot file not found: {filename}</p>"
        
        # Render template
        template = Template(template_str)
        html_content = template.render(
            analyses=analyses,
            decisions=decisions,
            decisions_dict=decisions_dict,
            summary_stats=summary_stats,
            plots=plots,
            **plot_content
        )
        
        # Write HTML file
        report_path = output_dir / "chimeric_detective_report.html"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        return str(report_path)