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

## License

MIT License - see LICENSE file for details.