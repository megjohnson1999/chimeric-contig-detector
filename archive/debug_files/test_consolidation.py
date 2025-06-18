#!/usr/bin/env python3
"""Test script to verify breakpoint consolidation is working."""

import sys
import json
from pathlib import Path

# Run the detector on test data
test_assembly = "../test_data_large/large_test_assembly.fasta"
test_reads_dir = "../test_data_large/reads"
output_dir = "test_consolidation_output"

# Check if test data exists
if not Path(test_assembly).exists():
    print(f"Error: Test assembly not found at {test_assembly}")
    sys.exit(1)

# Run chimeric detective
import subprocess

cmd = [
    "chimeric_detective",
    "--assembly", test_assembly,
    "--reads-dir", test_reads_dir,
    "--reads-pattern", "*_R{1,2}.fastq.gz",
    "--log-level", "INFO",
    "-o", output_dir
]

print(f"Running: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print(f"Error running chimeric_detective:")
    print(result.stderr)
    sys.exit(1)

# Check the results
results_file = Path(output_dir) / "chimeric_detective_results.json"
if results_file.exists():
    with open(results_file) as f:
        results = json.load(f)
    
    print(f"\n✅ Analysis completed successfully!")
    print(f"Total contigs: {results['summary']['total_contigs']}")
    print(f"Total analyses: {results['summary']['total_analyses']}")
    print(f"Average analyses per contig: {results['summary']['total_analyses'] / results['summary']['total_contigs']:.1f}")
    
    # Count analyses per contig
    contig_counts = {}
    for analysis in results['analyses']:
        contig_id = analysis['contig_id']
        contig_counts[contig_id] = contig_counts.get(contig_id, 0) + 1
    
    print(f"\nAnalyses per contig:")
    for contig_id, count in sorted(contig_counts.items()):
        print(f"  {contig_id}: {count} chimeric candidates")
        
    # Check for breakpoint regions and supporting candidates
    print(f"\nChecking for new consolidation features:")
    has_regions = False
    has_supporting = False
    for analysis in results['analyses']:
        if 'breakpoint_region' in analysis:
            has_regions = True
        if 'supporting_candidates' in analysis:
            has_supporting = True
    
    print(f"  Breakpoint regions: {'✅ Present' if has_regions else '❌ Missing'}")
    print(f"  Supporting candidates: {'✅ Present' if has_supporting else '❌ Missing'}")
    
else:
    print(f"Error: Results file not found at {results_file}")