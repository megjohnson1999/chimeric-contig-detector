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

### Single Sample Analysis

#### With BAM file
```bash
chimeric_detective --assembly viral_assembly.fasta --bam reads_aligned.bam --out results_dir
```

#### With paired-end reads
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads1 forward_reads.fastq.gz \
                  --reads2 reverse_reads.fastq.gz \
                  --out results_dir
```

#### With single-end reads
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads single_reads.fastq \
                  --out results_dir
```

### Multi-Sample Analysis

#### Directory with paired-end reads
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads-dir /path/to/reads_directory/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --out results_dir
```

#### Directory with single-end reads  
```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads-dir /path/to/reads_directory/ \
                  --reads-pattern "*.fastq.gz" \
                  --out results_dir
```

#### Custom file naming patterns

The `--reads-pattern` option uses `{1,2}` to indicate paired-end read files. Common patterns:

```bash
# Standard Illumina naming: sample_R1.fastq.gz, sample_R2.fastq.gz
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --out results/

# Underscore separated: sample_1.fq.gz, sample_2.fq.gz  
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*_{1,2}.fq.gz" \
                  --out results/

# Dot separated: sample.R1.fastq, sample.R2.fastq
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*.R{1,2}.fastq" \
                  --out results/

# Single-end reads (no {1,2}): sample.fastq.gz
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*.fastq.gz" \
                  --out results/
```

### Common Directory Structures

Your reads directory might look like:

```
# Paired-end example
reads_dir/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz  
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── sample3_R1.fastq.gz
└── sample3_R2.fastq.gz

# Single-end example
reads_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── sample3.fastq.gz

# Mixed compression
reads_dir/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq     # Uncompressed files also work
├── sample2_R2.fastq
```

### Multi-Sample Processing Modes

#### Analyze each sample separately (default)
```bash
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --multi-sample-mode separate \
                  --out results/
```

#### Merge all samples for higher coverage
```bash
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --multi-sample-mode merged \
                  --out results/
```

#### Parallel processing with custom workers
```bash
chimeric_detective --assembly assembly.fasta \
                  --reads-dir reads/ \
                  --max-workers 8 \
                  --parallel \
                  --out results/
```

## Output

### Single Sample Output

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

### Multi-Sample Output

When processing multiple samples, the output structure depends on the processing mode:

#### Separate Mode (default)
```
results_dir/
├── sample1/                      # Individual sample results
│   ├── cleaned_assembly.fasta
│   ├── chimeric_detective_report.html
│   └── ...
├── sample2/
│   ├── cleaned_assembly.fasta  
│   ├── chimeric_detective_report.html
│   └── ...
├── multi_sample_summary.tsv     # Summary across all samples
├── multi_sample_report.html     # Combined report
└── processing_log.txt           # Processing details
```

#### Merged Mode
```
results_dir/
├── merged_analysis/             # Combined analysis results
│   ├── cleaned_assembly.fasta   # Assembly with all samples merged
│   ├── chimeric_detective_report.html
│   └── ...
└── sample_contributions.tsv     # Which samples contributed to each detection
```

### Key Output Files

| File | Description |
|------|-------------|
| `cleaned_assembly.fasta` | Assembly with technical chimeras split |
| `chimeric_detective_report.html` | Interactive visualization and analysis |
| `chimeric_detective_results.json` | Machine-readable results for downstream analysis |
| `splitting_decisions.tsv` | Detailed table of all modifications made |
| `multi_sample_summary.tsv` | Cross-sample comparison (multi-sample only) |

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

### Processing Mode Details

#### When to Use Each Mode

- **Separate mode**: Analyze each sample independently (default)
  - Best for: Comparing contamination across samples, sample-specific analysis
  - Output: Individual results for each sample + combined summary

- **Merged mode**: Combine all reads for higher coverage analysis  
  - Best for: Low-coverage samples, co-assembly analysis
  - Output: Single analysis with all samples combined

- **Batch mode**: Memory-efficient processing of many samples
  - Best for: Large numbers of samples with limited RAM
  - Output: Processed in chunks to manage memory usage

```bash
# Memory-efficient batch processing
chimeric_detective -a assembly.fasta --reads-dir samples/ \
                  --multi-sample-mode batch \
                  --batch-size 5 \
                  --max-workers 4 -o results/
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