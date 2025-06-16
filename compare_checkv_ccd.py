#!/usr/bin/env python3
"""
Compare CheckV and Chimeric Detective results to understand assembly quality.
Identifies overlapping issues and categorizes different types of problems.
"""

import pandas as pd
import numpy as np
import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Set, Tuple, List

def load_checkv_results(checkv_dir: str) -> Dict[str, pd.DataFrame]:
    """Load CheckV result files."""
    checkv_path = Path(checkv_dir)
    results = {}
    
    # Quality summary - main results
    quality_file = checkv_path / "quality_summary.tsv"
    if quality_file.exists():
        results['quality'] = pd.read_csv(quality_file, sep='\t')
        print(f"Loaded CheckV quality summary: {len(results['quality'])} contigs")
    else:
        print(f"Warning: {quality_file} not found")
    
    # Contamination details
    contamination_file = checkv_path / "contamination.tsv"
    if contamination_file.exists():
        results['contamination'] = pd.read_csv(contamination_file, sep='\t')
        print(f"Loaded CheckV contamination: {len(results['contamination'])} contigs")
    else:
        print(f"Warning: {contamination_file} not found")
    
    # Complete genomes
    complete_file = checkv_path / "complete_genomes.tsv"
    if complete_file.exists():
        results['complete'] = pd.read_csv(complete_file, sep='\t')
        print(f"Loaded CheckV complete genomes: {len(results['complete'])} contigs")
    
    # Proviruses
    proviruses_file = checkv_path / "proviruses.tsv"
    if proviruses_file.exists():
        results['proviruses'] = pd.read_csv(proviruses_file, sep='\t')
        print(f"Loaded CheckV proviruses: {len(results['proviruses'])} contigs")
    
    return results

def load_ccd_results(ccd_dir: str) -> Dict[str, pd.DataFrame]:
    """Load Chimeric Detective result files."""
    ccd_path = Path(ccd_dir)
    results = {}
    
    # Splitting decisions - main results
    decisions_file = ccd_path / "splitting_decisions.tsv"
    if decisions_file.exists():
        results['decisions'] = pd.read_csv(decisions_file, sep='\t')
        print(f"Loaded CCD splitting decisions: {len(results['decisions'])} contigs")
    else:
        print(f"Warning: {decisions_file} not found")
    
    # Try to find results JSON for more detailed info
    json_files = list(ccd_path.glob("*results.json"))
    if json_files:
        import json
        with open(json_files[0], 'r') as f:
            results['detailed'] = json.load(f)
        print(f"Loaded CCD detailed results from {json_files[0].name}")
    
    return results

def categorize_checkv_quality(row) -> str:
    """Categorize CheckV quality assessment."""
    if pd.isna(row.get('checkv_quality')):
        return 'Unknown'
    
    quality = row['checkv_quality'].lower()
    contamination = row.get('contamination', 0)
    
    if 'complete' in quality:
        return 'High-quality'
    elif 'high' in quality:
        return 'High-quality'
    elif 'medium' in quality:
        return 'Medium-quality'
    elif 'low' in quality:
        return 'Low-quality'
    elif contamination > 0:
        return 'Contaminated'
    else:
        return 'Not-determined'

def categorize_ccd_action(row) -> str:
    """Categorize CCD action taken."""
    if pd.isna(row.get('action')):
        return 'Not-analyzed'
    
    action = row['action'].lower()
    if 'split' in action:
        return 'Split'
    elif 'preserve' in action:
        return 'Preserved'
    elif 'flag' in action:
        return 'Flagged'
    else:
        return 'Other'

def analyze_overlap(checkv_data: Dict, ccd_data: Dict) -> Dict:
    """Analyze overlap between CheckV and CCD results."""
    analysis = {}
    
    # Get contig sets
    checkv_contigs = set()
    ccd_contigs = set()
    
    if 'quality' in checkv_data:
        checkv_contigs.update(checkv_data['quality']['contig_id'].astype(str))
    
    if 'decisions' in ccd_data:
        ccd_contigs.update(ccd_data['decisions']['contig_id'].astype(str))
    
    analysis['total_checkv'] = len(checkv_contigs)
    analysis['total_ccd'] = len(ccd_contigs)
    analysis['overlap'] = len(checkv_contigs.intersection(ccd_contigs))
    analysis['checkv_only'] = len(checkv_contigs - ccd_contigs)
    analysis['ccd_only'] = len(ccd_contigs - checkv_contigs)
    
    # Problematic contigs
    checkv_problems = set()
    ccd_problems = set()
    
    # CheckV problematic contigs
    if 'contamination' in checkv_data:
        contaminated = checkv_data['contamination']
        if 'region_types' in contaminated.columns:
            # Contigs with mixed regions
            mixed_contigs = contaminated[
                contaminated['region_types'].str.contains('host', na=False)
            ]['contig_id'].astype(str)
            checkv_problems.update(mixed_contigs)
    
    if 'quality' in checkv_data:
        quality = checkv_data['quality']
        # Low quality or contaminated contigs
        poor_quality = quality[
            (quality['checkv_quality'].str.contains('Low-quality|Not-determined', na=False)) |
            (quality['contamination'] > 0)
        ]['contig_id'].astype(str)
        checkv_problems.update(poor_quality)
    
    # CCD problematic contigs
    if 'decisions' in ccd_data:
        decisions = ccd_data['decisions']
        # Split or flagged contigs
        ccd_problems.update(
            decisions[decisions['action'].isin(['split', 'flag'])]['contig_id'].astype(str)
        )
    
    analysis['checkv_problems'] = len(checkv_problems)
    analysis['ccd_problems'] = len(ccd_problems)
    analysis['both_problems'] = len(checkv_problems.intersection(ccd_problems))
    
    return analysis

def create_comparison_table(checkv_data: Dict, ccd_data: Dict) -> pd.DataFrame:
    """Create detailed comparison table."""
    # Start with all contigs from both tools
    all_contigs = set()
    
    if 'quality' in checkv_data:
        all_contigs.update(checkv_data['quality']['contig_id'].astype(str))
    
    if 'decisions' in ccd_data:
        all_contigs.update(ccd_data['decisions']['contig_id'].astype(str))
    
    comparison = []
    
    for contig_id in all_contigs:
        row = {'contig_id': contig_id}
        
        # CheckV data
        if 'quality' in checkv_data:
            checkv_row = checkv_data['quality'][
                checkv_data['quality']['contig_id'].astype(str) == contig_id
            ]
            if not checkv_row.empty:
                checkv_row = checkv_row.iloc[0]
                row['checkv_quality'] = checkv_row.get('checkv_quality', 'N/A')
                row['checkv_completeness'] = checkv_row.get('completeness', 0)
                row['checkv_contamination'] = checkv_row.get('contamination', 0)
                row['contig_length'] = checkv_row.get('contig_length', 0)
                row['checkv_category'] = categorize_checkv_quality(checkv_row)
            else:
                row['checkv_quality'] = 'Not-analyzed'
                row['checkv_completeness'] = 0
                row['checkv_contamination'] = 0
                row['contig_length'] = 0
                row['checkv_category'] = 'Not-analyzed'
        
        # CCD data
        if 'decisions' in ccd_data:
            ccd_row = ccd_data['decisions'][
                ccd_data['decisions']['contig_id'].astype(str) == contig_id
            ]
            if not ccd_row.empty:
                ccd_row = ccd_row.iloc[0]
                row['ccd_action'] = ccd_row.get('action', 'N/A')
                row['ccd_confidence'] = ccd_row.get('confidence', 0)
                row['ccd_reason'] = ccd_row.get('reason', 'N/A')
                row['ccd_category'] = categorize_ccd_action(ccd_row)
                if 'original_length' in ccd_row:
                    row['original_length'] = ccd_row['original_length']
            else:
                row['ccd_action'] = 'Not-analyzed'
                row['ccd_confidence'] = 0
                row['ccd_reason'] = 'N/A'
                row['ccd_category'] = 'Not-analyzed'
        
        comparison.append(row)
    
    return pd.DataFrame(comparison)

def generate_summary_report(comparison_df: pd.DataFrame, analysis: Dict, output_dir: str):
    """Generate comprehensive summary report."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Summary statistics
    print("\n" + "="*60)
    print("📊 CHECKV vs CHIMERIC DETECTIVE COMPARISON")
    print("="*60)
    
    print(f"\n🔍 Dataset Overview:")
    print(f"  Total unique contigs: {len(comparison_df)}")
    print(f"  CheckV analyzed: {analysis['total_checkv']}")
    print(f"  CCD analyzed: {analysis['total_ccd']}")
    print(f"  Both tools analyzed: {analysis['overlap']}")
    
    print(f"\n⚠️  Problematic Contigs:")
    print(f"  CheckV flagged issues: {analysis['checkv_problems']}")
    print(f"  CCD detected chimeras: {analysis['ccd_problems']}")
    print(f"  Both tools flagged: {analysis['both_problems']}")
    
    # Agreement analysis
    both_analyzed = comparison_df[
        (comparison_df['checkv_category'] != 'Not-analyzed') & 
        (comparison_df['ccd_category'] != 'Not-analyzed')
    ]
    
    if len(both_analyzed) > 0:
        print(f"\n🤝 Agreement Analysis (n={len(both_analyzed)} contigs):")
        
        # High confidence problems (both tools agree)
        high_conf_problems = both_analyzed[
            (both_analyzed['checkv_category'].isin(['Low-quality', 'Contaminated'])) &
            (both_analyzed['ccd_category'] == 'Split')
        ]
        print(f"  High confidence issues: {len(high_conf_problems)} contigs")
        
        # Conservative splitting (CCD preserves what CheckV likes)
        conservative = both_analyzed[
            (both_analyzed['checkv_category'].isin(['High-quality', 'Medium-quality'])) &
            (both_analyzed['ccd_category'] == 'Preserved')
        ]
        print(f"  Conservative preservation: {len(conservative)} contigs")
        
        # Disagreements
        disagreements = both_analyzed[
            ((both_analyzed['checkv_category'].isin(['High-quality', 'Medium-quality'])) &
             (both_analyzed['ccd_category'] == 'Split')) |
            ((both_analyzed['checkv_category'].isin(['Low-quality', 'Contaminated'])) &
             (both_analyzed['ccd_category'] == 'Preserved'))
        ]
        print(f"  Disagreements: {len(disagreements)} contigs")
    
    # Category breakdown
    print(f"\n📋 Category Breakdown:")
    print("\nCheckV Categories:")
    checkv_cats = comparison_df['checkv_category'].value_counts()
    for cat, count in checkv_cats.items():
        print(f"  {cat}: {count}")
    
    print("\nCCD Categories:")
    ccd_cats = comparison_df['ccd_category'].value_counts()
    for cat, count in ccd_cats.items():
        print(f"  {cat}: {count}")
    
    # Write detailed TSV
    comparison_file = output_path / "checkv_ccd_comparison.tsv"
    comparison_df.to_csv(comparison_file, sep='\t', index=False)
    print(f"\n📄 Detailed comparison saved: {comparison_file}")
    
    # High-value targets for manual inspection
    manual_check = comparison_df[
        ((comparison_df['checkv_category'] == 'High-quality') & 
         (comparison_df['ccd_category'] == 'Split')) |
        ((comparison_df['checkv_category'] == 'Contaminated') & 
         (comparison_df['ccd_category'] == 'Preserved')) |
        ((comparison_df['checkv_contamination'] > 0) & 
         (comparison_df['ccd_confidence'] > 0.8))
    ]
    
    if len(manual_check) > 0:
        manual_file = output_path / "manual_inspection_targets.tsv"
        manual_check.to_csv(manual_file, sep='\t', index=False)
        print(f"🔍 Manual inspection targets: {manual_file} ({len(manual_check)} contigs)")
    
    # Summary statistics file
    summary_stats = {
        'total_contigs': len(comparison_df),
        'checkv_analyzed': analysis['total_checkv'],
        'ccd_analyzed': analysis['total_ccd'],
        'both_analyzed': analysis['overlap'],
        'checkv_problems': analysis['checkv_problems'],
        'ccd_problems': analysis['ccd_problems'],
        'both_flagged': analysis['both_problems'],
        'high_confidence_issues': len(high_conf_problems) if len(both_analyzed) > 0 else 0,
        'disagreements': len(disagreements) if len(both_analyzed) > 0 else 0
    }
    
    summary_file = output_path / "comparison_summary.tsv"
    pd.DataFrame([summary_stats]).to_csv(summary_file, sep='\t', index=False)
    print(f"📊 Summary statistics: {summary_file}")

def create_visualizations(comparison_df: pd.DataFrame, output_dir: str):
    """Create comparison visualizations."""
    output_path = Path(output_dir)
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('CheckV vs Chimeric Detective Comparison', fontsize=16, fontweight='bold')
    
    # 1. Category overlap heatmap
    ax1 = axes[0, 0]
    pivot_data = comparison_df.groupby(['checkv_category', 'ccd_category']).size().unstack(fill_value=0)
    sns.heatmap(pivot_data, annot=True, fmt='d', cmap='Blues', ax=ax1)
    ax1.set_title('Tool Agreement Matrix')
    ax1.set_xlabel('CCD Category')
    ax1.set_ylabel('CheckV Category')
    
    # 2. Quality vs confidence scatter
    ax2 = axes[0, 1]
    both_analyzed = comparison_df[
        (comparison_df['checkv_completeness'] > 0) & 
        (comparison_df['ccd_confidence'] > 0)
    ]
    if len(both_analyzed) > 0:
        scatter = ax2.scatter(both_analyzed['checkv_completeness'], 
                             both_analyzed['ccd_confidence'],
                             c=both_analyzed['contig_length'], 
                             cmap='viridis', alpha=0.6)
        ax2.set_xlabel('CheckV Completeness')
        ax2.set_ylabel('CCD Confidence')
        ax2.set_title('Quality vs Confidence')
        plt.colorbar(scatter, ax=ax2, label='Contig Length')
    
    # 3. Length distribution by category
    ax3 = axes[1, 0]
    categories = ['High-quality', 'Medium-quality', 'Low-quality', 'Split', 'Preserved']
    plot_data = []
    labels = []
    
    for cat in categories:
        if cat in ['High-quality', 'Medium-quality', 'Low-quality']:
            data = comparison_df[comparison_df['checkv_category'] == cat]['contig_length']
            labels.append(f'CheckV {cat}')
        else:
            data = comparison_df[comparison_df['ccd_category'] == cat]['contig_length']
            labels.append(f'CCD {cat}')
        
        if len(data) > 0:
            plot_data.append(data.dropna())
    
    if plot_data:
        ax3.boxplot(plot_data, labels=[l.replace(' ', '\n') for l in labels])
        ax3.set_ylabel('Contig Length (bp)')
        ax3.set_title('Length Distribution by Category')
        ax3.tick_params(axis='x', rotation=45)
    
    # 4. Problem detection comparison
    ax4 = axes[1, 1]
    problem_data = {
        'CheckV Only': analysis['checkv_problems'] - analysis['both_problems'],
        'Both Tools': analysis['both_problems'],
        'CCD Only': analysis['ccd_problems'] - analysis['both_problems']
    }
    
    bars = ax4.bar(problem_data.keys(), problem_data.values(), 
                   color=['lightblue', 'orange', 'lightgreen'])
    ax4.set_ylabel('Number of Contigs')
    ax4.set_title('Problem Detection Overlap')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = output_path / "checkv_ccd_comparison.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"📈 Visualizations saved: {plot_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare CheckV and Chimeric Detective results')
    parser.add_argument('--checkv-dir', required=True, help='CheckV output directory')
    parser.add_argument('--ccd-dir', required=True, help='Chimeric Detective output directory')
    parser.add_argument('--output-dir', default='comparison_results', help='Output directory for comparison results')
    
    args = parser.parse_args()
    
    print("🔍 Loading CheckV results...")
    checkv_data = load_checkv_results(args.checkv_dir)
    
    print("\n🔍 Loading Chimeric Detective results...")
    ccd_data = load_ccd_results(args.ccd_dir)
    
    if not checkv_data and not ccd_data:
        print("❌ No valid results found in either directory!")
        sys.exit(1)
    
    print("\n📊 Analyzing overlap...")
    global analysis  # Make available for plotting function
    analysis = analyze_overlap(checkv_data, ccd_data)
    
    print("\n📋 Creating comparison table...")
    comparison_df = create_comparison_table(checkv_data, ccd_data)
    
    print("\n📄 Generating summary report...")
    generate_summary_report(comparison_df, analysis, args.output_dir)
    
    print("\n📈 Creating visualizations...")
    try:
        create_visualizations(comparison_df, args.output_dir)
    except Exception as e:
        print(f"⚠️  Visualization creation failed: {e}")
        print("   (Continuing with text-based analysis)")
    
    print(f"\n✅ Comparison complete! Results in: {args.output_dir}")

if __name__ == "__main__":
    main()