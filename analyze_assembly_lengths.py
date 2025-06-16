#!/usr/bin/env python3
"""
Analyze assembly contig lengths to guide minimum length threshold selection.
"""

import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def analyze_assembly_lengths(assembly_file):
    """Analyze contig length distribution."""
    lengths = []
    
    print(f"Analyzing assembly: {assembly_file}")
    
    for record in SeqIO.parse(assembly_file, "fasta"):
        lengths.append(len(record.seq))
    
    lengths = np.array(lengths)
    
    print(f"\n📊 Assembly Statistics:")
    print(f"Total contigs: {len(lengths):,}")
    print(f"Total length: {lengths.sum():,} bp")
    print(f"Mean length: {lengths.mean():.1f} bp")
    print(f"Median length: {np.median(lengths):.1f} bp")
    print(f"N50: {calculate_n50(lengths):,} bp")
    
    # Analyze different thresholds
    thresholds = [500, 1000, 1500, 2000, 3000, 5000]
    
    print(f"\n🎯 Impact of Different Minimum Length Thresholds:")
    print(f"{'Threshold':<10} {'Contigs':<10} {'% Kept':<10} {'Total bp':<15} {'% Length':<10}")
    print("-" * 60)
    
    for threshold in thresholds:
        kept_contigs = lengths[lengths >= threshold]
        pct_contigs = (len(kept_contigs) / len(lengths)) * 100
        pct_length = (kept_contigs.sum() / lengths.sum()) * 100
        
        print(f"{threshold:<10} {len(kept_contigs):<10,} {pct_contigs:<10.1f} {kept_contigs.sum():<15,} {pct_length:<10.1f}")
    
    # Specific analysis for chimera detection relevance
    print(f"\n🔬 Chimera Detection Relevance:")
    long_contigs = lengths[lengths >= 2000]
    very_long_contigs = lengths[lengths >= 5000]
    
    print(f"Contigs ≥2000bp: {len(long_contigs):,} ({len(long_contigs)/len(lengths)*100:.1f}%)")
    print(f"Contigs ≥5000bp: {len(very_long_contigs):,} ({len(very_long_contigs)/len(lengths)*100:.1f}%)")
    print(f"These longer contigs contain {long_contigs.sum()/lengths.sum()*100:.1f}% of total sequence")
    
    return lengths

def calculate_n50(lengths):
    """Calculate N50 statistic."""
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length / 2
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= target_length:
            return length
    return sorted_lengths[-1]

def plot_length_distribution(lengths, output_file="contig_lengths.png"):
    """Plot contig length distribution."""
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Histogram
    plt.subplot(2, 2, 1)
    plt.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Contig Length (bp)')
    plt.ylabel('Count')
    plt.title('Contig Length Distribution')
    plt.yscale('log')
    
    # Subplot 2: Log scale histogram
    plt.subplot(2, 2, 2)
    plt.hist(lengths, bins=np.logspace(2, 6, 50), alpha=0.7, edgecolor='black')
    plt.xlabel('Contig Length (bp)')
    plt.ylabel('Count')
    plt.title('Contig Length Distribution (Log Scale)')
    plt.xscale('log')
    plt.yscale('log')
    
    # Subplot 3: Cumulative plot
    plt.subplot(2, 2, 3)
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative_pct = np.cumsum(sorted_lengths) / np.sum(sorted_lengths) * 100
    plt.plot(range(len(sorted_lengths)), cumulative_pct)
    plt.xlabel('Contig Rank')
    plt.ylabel('Cumulative % of Total Length')
    plt.title('Cumulative Length Distribution')
    
    # Subplot 4: Length thresholds
    plt.subplot(2, 2, 4)
    thresholds = [500, 1000, 1500, 2000, 3000, 5000, 10000]
    contigs_kept = [len(lengths[lengths >= t]) for t in thresholds]
    
    plt.bar(range(len(thresholds)), contigs_kept, alpha=0.7)
    plt.xticks(range(len(thresholds)), [f"{t/1000:.1f}k" for t in thresholds], rotation=45)
    plt.xlabel('Minimum Length Threshold')
    plt.ylabel('Contigs Remaining')
    plt.title('Impact of Length Filtering')
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n📈 Length distribution plot saved: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_assembly_lengths.py <assembly.fasta>")
        sys.exit(1)
    
    assembly_file = sys.argv[1]
    lengths = analyze_assembly_lengths(assembly_file)
    
    # Only plot if matplotlib is available
    try:
        plot_length_distribution(lengths)
    except ImportError:
        print("Note: matplotlib not available, skipping plots")