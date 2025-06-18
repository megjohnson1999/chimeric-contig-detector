# Simplified GC-Only Chimera Detection

This is a simplified version of Chimeric Detective that uses only GC content analysis to detect chimeric contigs. It's designed to be fast, simple, and easy to understand.

## Quick Start

```bash
# Install (if not already installed)
pip install -e .

# Run simplified detection
chimeric_detective_simple assembly.fasta -o output_dir

# Or use the Python script directly
python examples/simple_gc_detection.py assembly.fasta
```

## Features

- **Single detection method**: Uses only GC content changes to identify breakpoints
- **Fast processing**: No alignment or k-mer analysis required
- **Simple output**: Clear visualization of GC content and detected breakpoints
- **Minimal dependencies**: Works with just the assembly file

## How It Works

1. **Sliding window analysis**: Calculates GC content in windows along each contig
2. **Change detection**: Identifies significant changes in GC content between adjacent windows
3. **Classification**: Determines if changes are likely technical artifacts
4. **Resolution**: Optionally splits contigs at high-confidence breakpoints

## Parameters

- `--gc-threshold`: Minimum GC difference to consider (default: 0.1)
- `--window-size`: Size of sliding window (default: 1000)
- `--step-size`: Step size for sliding window (default: 500)
- `--min-confidence`: Minimum confidence for splitting (default: 0.7)
- `--no-split`: Only detect, don't split contigs
- `--no-visualize`: Skip HTML report generation

## Example Output

The tool generates:
- List of detected chimeras with confidence scores
- Optional: Split contigs in FASTA format
- Optional: HTML report with GC content plots

## When to Use This Version

Use the simplified version when:
- You want fast, straightforward chimera detection
- You don't have alignment data available
- You're primarily concerned with GC content discontinuities
- You want to understand the basic detection algorithm

## Differences from Full Version

| Feature | Simple Version | Full Version |
|---------|---------------|--------------|
| Detection Methods | GC only | GC, coverage, k-mer, orientation |
| Input Requirements | Assembly only | Assembly + reads/alignment |
| Processing Time | Fast | Slower but more comprehensive |
| Accuracy | Good for obvious chimeras | Better for subtle chimeras |
| Complexity | Simple, easy to modify | Complex, more parameters |

## Python API

```python
from chimeric_detective import SimpleChimeraDetector, SimpleChimeraAnalyzer

# Detect chimeras
detector = SimpleChimeraDetector(gc_content_threshold=0.1)
candidates = detector.detect_chimeras("assembly.fasta")

# Analyze results
analyzer = SimpleChimeraAnalyzer()
classifications = analyzer.analyze_candidates(candidates)

# Filter high confidence
high_confidence = analyzer.filter_high_confidence(classifications, min_confidence=0.7)
```