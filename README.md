# Chimeric Detective

[![CI](https://github.com/megjohnson1999/chimeric-contig-detector/workflows/CI/badge.svg)](https://github.com/megjohnson1999/chimeric-contig-detector/actions)
[![codecov](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector/branch/main/graph/badge.svg)](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector)
[![PyPI version](https://badge.fury.io/py/chimeric-detective.svg)](https://badge.fury.io/py/chimeric-detective)
[![Python](https://img.shields.io/pypi/pyversions/chimeric-detective.svg)](https://pypi.org/project/chimeric-detective/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive command-line tool for detecting, analyzing, explaining, and resolving chimeric contigs in viral metagenomic assemblies.

## Features

- **Chimera Detection**: Multiple methods including coverage discontinuities, taxonomic transitions, and sequence composition shifts
- **Precise Breakpoint Identification**: Accurately locate chimeric junctions
- **Intelligent Classification**: Distinguish between technical artifacts and biological recombination
- **Automated Resolution**: Split technical chimeras while preserving biological events
- **Interactive Visualization**: Generate HTML reports with detailed visualizations
- **Flexible Input**: Support for FASTA assemblies with FASTQ reads or BAM alignments

## Installation

### Option 1: Conda (Recommended)

Conda will automatically install all dependencies including external bioinformatics tools:

```bash
# Clone the repository
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector

# Create and activate conda environment
conda env create -f environment.yml
conda activate chimeric-detective

# Install the package
pip install -e .
```

### Option 2: Pip + Manual Dependencies

```bash
# Install external tools manually (Ubuntu/Debian)
sudo apt-get install bwa minimap2 samtools ncbi-blast+

# Or on macOS with Homebrew
brew install bwa minimap2 samtools blast

# Install Python package
pip install chimeric-detective
```

### Option 3: From Source

```bash
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector
pip install -e .

# You'll need to install external tools separately (see Option 2)
```

### Dependencies

**Python packages** (installed automatically):
- biopython, pysam, numpy, pandas, scipy
- scikit-learn, matplotlib, seaborn, plotly
- click, tqdm, jinja2

**External tools** (required for functionality):
- BWA or minimap2 (for read alignment)
- samtools (for BAM file processing)
- BLAST (optional, for taxonomic classification)

## Quick Start

### Basic usage with BAM file
```bash
chimeric_detective --assembly viral_assembly.fasta --bam reads_aligned.bam --out results_dir
```

### Usage with raw reads (paired-end)
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads1 forward_reads.fastq.gz \
                  --reads2 reverse_reads.fastq.gz \
                  --out results_dir
```

### Usage with single-end reads
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads single_reads.fastq \
                  --out results_dir
```

## Output

The tool creates a structured output directory:

```
results_dir/
├── cleaned_assembly.fasta        # Assembly with technical chimeras split
├── chimeric_contigs/             # Individual files for each chimeric contig
├── chimeric_detective_report.html # Interactive HTML report
├── chimeric_detective_results.json # Machine-readable results
├── splitting_decisions.tsv       # Table of all modifications made
└── figures/                      # Static visualizations
```

## Benchmarking and Evaluation

Chimeric Detective includes a comprehensive benchmarking suite to compare performance against other contig-level quality assessment tools.

### Quick Benchmark

For rapid evaluation against other tools:

```bash
cd benchmarking/
./quick_benchmark.sh -a your_assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz
```

This compares Chimeric Detective against:
- **CheckV** - Viral contamination detection
- **VirSorter2** - Low-confidence viral sequence identification  
- **QUAST** - General assembly quality metrics

### Comprehensive Benchmark

For detailed analysis with synthetic test data:

```bash
# Create synthetic data with known chimeras
python benchmarking/benchmark_chimeric_detective.py --create-test-data -o benchmark_results

# Or with your real data
python benchmarking/benchmark_chimeric_detective.py \
    -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz \
    -o benchmark_results
```

### Installing Comparison Tools

```bash
# Install comparison tools via conda
conda install -c bioconda checkv virsorter=2 quast bwa samtools
```

### Benchmark Output

The benchmarking generates:
- **Summary table** (TSV format) with detection counts and runtimes
- **Detailed report** (Markdown) with tool-specific analysis
- **Individual tool outputs** for manual inspection

### Interpreting Results

| Tool | Purpose | Interpretation |
|------|---------|----------------|
| Chimeric Detective | Multi-method chimera detection | Comprehensive analysis with breakpoint resolution |
| CheckV | Viral contamination | Host gene contamination in viral contigs |
| VirSorter2 | Viral identification | Low-confidence sequences (potential chimeras) |
| QUAST | Assembly quality | General misassemblies and structural issues |

**Important**: Higher chimera counts don't necessarily indicate better performance. Consider:
- False positive rates vs. sensitivity
- Biological relevance of detected chimeras
- Tool-specific definitions of "chimeric"
- Manual validation of results

### Example Benchmark Results

```
Tool                Status   Runtime(s)  Chimeras_Detected  Notes
Chimeric_Detective  SUCCESS  45          12                 Comprehensive analysis
CheckV              SUCCESS  120         3                  Conservative, viral-specific  
VirSorter2          SUCCESS  90          8                  Low-confidence sequences
QUAST               SUCCESS  15          5                  General misassemblies
```

For detailed benchmarking documentation, see [`benchmarking/README.md`](benchmarking/README.md).

## Advanced Usage

### Multi-sample Processing

Process multiple samples in parallel:

```bash
# Analyze each sample separately
chimeric_detective -a assembly.fasta --reads-dir samples/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" -o results/

# Merge all samples for higher coverage
chimeric_detective -a assembly.fasta --reads-dir samples/ \
                  --multi-sample-mode merged -o results/

# Parallel processing with custom workers
chimeric_detective -a assembly.fasta --reads-dir samples/ \
                  --max-workers 8 --parallel -o results/
```

### Sensitivity Settings

```bash
# High sensitivity (more chimeras detected, potential false positives)
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity high

# Low sensitivity (conservative, fewer false positives)
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity low
```

### Custom Thresholds

```bash
chimeric_detective -a assembly.fasta -b reads.bam -o results/ \
                  --min-coverage 10.0 \
                  --coverage-fold-change 3.0 \
                  --gc-content-threshold 0.15 \
                  --confidence-threshold 0.7
```

## License

MIT License - see LICENSE file for details.