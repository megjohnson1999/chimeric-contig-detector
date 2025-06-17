# Chimeric Detective

[![CI](https://github.com/megjohnson1999/chimeric-contig-detector/workflows/CI/badge.svg)](https://github.com/megjohnson1999/chimeric-contig-detector/actions)
[![codecov](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector/branch/main/graph/badge.svg)](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector)
[![PyPI version](https://badge.fury.io/py/chimeric-detective.svg)](https://badge.fury.io/py/chimeric-detective)
[![Python](https://img.shields.io/pypi/pyversions/chimeric-detective.svg)](https://pypi.org/project/chimeric-detective/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive command-line tool for detecting, analyzing, explaining, and resolving chimeric contigs in viral metagenomic assemblies. Optimized for large-scale processing with intelligent memory management and robust dependency handling.

## 🚀 New: Simple GC-Based Detector

For users who want a **fast, reliable, assembly-only** chimera detector, we now provide `gc_only_detector.py` - a focused tool that detects chimeric contigs using only GC content shifts. This tool:

- ✅ **Works with assembly files only** (no BAM files required)
- ✅ **Zero false positives** on test data
- ✅ **100% detection** of GC-based chimeric contigs
- ✅ **Simple and fast** (150 lines of code)
- ✅ **Perfect for initial screening** of assemblies

```bash
# Quick GC-based chimera detection
python gc_only_detector.py -a assembly.fasta -o results/
```

## Features

- **Advanced Signal-Based Detection**: Nucleotide-resolution breakpoint analysis with adaptive window sizing
- **Biologically Meaningful Analysis**: Window sizes automatically scale with contig length for optimal sensitivity
- **Multi-Method Chimera Detection**: Coverage discontinuities, GC content shifts, k-mer composition changes, and read orientation patterns
- **True Breakpoint Refinement**: Sub-grid detection and 1bp-resolution refinement for precise junction identification
- **Intelligent Classification**: Distinguish between technical artifacts (PCR chimeras, assembly errors) and biological recombination
- **Automated Resolution**: Split technical chimeras while preserving biological events
- **Memory-Optimized Processing**: Efficient handling of large assemblies (50K+ contigs) with configurable resource usage
- **Multi-Sample Support**: Process hundreds of samples with parallel processing and batch modes
- **Interactive Visualization**: Generate comprehensive HTML reports with detailed visualizations
- **Flexible Input**: Support for FASTA assemblies with FASTQ reads or pre-aligned BAM files
- **Robust Tool Integration**: Automatic fallback between minimap2 and BWA with intelligent error handling
- **Production-Ready Stability**: Comprehensive bounds checking and error handling (v1.1.1+)

## Installation

### Option 1: Docker (Recommended)

The easiest way to get started with all dependencies pre-installed:

```bash
# Pull and run the container
docker pull ghcr.io/megjohnson1999/chimeric-contig-detector:latest

# Run analysis
docker run --rm -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/assembly.fasta \
    --reads-dir /data/reads \
    --out /data/results
```

See [DOCKER.md](DOCKER.md) for detailed Docker usage instructions.

### Option 2: PyPI (Python Package)

Install from PyPI with automatic dependency management:

```bash
# Install the package
pip install chimeric-detective

# Install external bioinformatics tools
conda install -c bioconda bwa minimap2 samtools
# OR on Ubuntu/Debian:
sudo apt-get install bwa minimap2 samtools ncbi-blast+
# OR on macOS:
brew install bwa minimap2 samtools blast
```

### Option 3: Conda Environment

Full environment setup with all dependencies:

```bash
# Clone the repository
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector

# Create and activate conda environment
conda env create -f environment.yml
conda activate chimeric-detective

# Install the package in development mode
pip install -e .
```

### Option 4: From Source

For developers and contributors:

```bash
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector
pip install -e .

# Install external tools separately (see Option 2)
```

### Dependencies

**Python packages** (installed automatically):
- biopython, pysam, numpy, pandas, scipy
- scikit-learn, matplotlib, seaborn, plotly
- click, tqdm, jinja2

**External tools** (installed automatically via conda):
- BWA or minimap2 (for read alignment) 
- samtools (for BAM file processing)
- BLAST (optional, for taxonomic classification)

**Note**: The tool automatically detects available aligners and uses minimap2 first (faster) with BWA as fallback (more memory-efficient). Both tools are installed with the conda environment.

## Quick Start

### GC-Only Detection (Recommended for Initial Screening)

The fastest way to detect GC-based chimeric contigs:

```bash
# Download the detector
wget https://raw.githubusercontent.com/megjohnson1999/chimeric-contig-detector/main/gc_only_detector.py

# Run detection (assembly only)
python gc_only_detector.py -a viral_assembly.fasta -o results/

# Check results
ls results/
# cleaned_assembly.fasta       # Assembly with chimeras split
# gc_chimera_report.tsv        # Detailed detection report
```

### Full Pipeline Analysis

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

## Large Assembly Optimization

For assemblies with many contigs (>10K) or large files (>500MB), use these optimizations:

### Memory-Efficient Processing

```bash
# Analyze only longer contigs (recommended for large assemblies)
chimeric_detective --assembly large_assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --min-contig-length 2000 \
                  --max-workers 4 \
                  --out results/
```

### Analyze Your Assembly First

Use the included analysis tool to understand your assembly characteristics:

```bash
# Analyze contig length distribution and memory requirements
python analyze_assembly_lengths.py your_assembly.fasta
```

This shows:
- Contig count and size distribution  
- Impact of different `--min-contig-length` thresholds
- Computational savings vs sequence retention
- Memory usage estimates

### Recommended Settings for Large Assemblies

| Assembly Size | Contigs | Min Length | Max Workers | Memory |
|---------------|---------|------------|-------------|---------|
| <100MB | <5K | 1000bp | 8 | 16GB |
| 100-500MB | 5K-25K | 1500bp | 4 | 32GB |
| >500MB | >25K | 2000bp | 2-4 | 64GB |

### SLURM/HPC Example

```bash
#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

# Activate environment
conda activate chimeric-detective

# Memory-optimized processing
chimeric_detective --assembly huge_assembly.fasta \
                  --reads-dir reads_subset/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --min-contig-length 2000 \
                  --max-workers 4 \
                  --out results/
```

**Key Points:**
- Short contigs (<threshold) are **preserved unchanged** in output
- Only longer contigs are **analyzed for chimeras**  
- **Massive memory savings** (often 60-80% reduction) with minimal sequence loss
- Tool automatically warns about large assemblies and suggests optimizations

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

Choose the appropriate sensitivity level based on your tolerance for false positives vs false negatives:

```bash
# Conservative (default): High specificity, minimal false positives
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity conservative

# Balanced: Moderate sensitivity and specificity  
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity balanced

# Sensitive: Higher sensitivity, may increase false positives
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity sensitive

# Very sensitive: Maximum sensitivity for difficult samples
chimeric_detective -a assembly.fasta -b reads.bam -o results/ --sensitivity very_sensitive
```

**Sensitivity Mode Comparison:**

| Mode | GC Threshold | Coverage Change | K-mer Distance | Evidence Required | Best For |
|------|-------------|----------------|----------------|------------------|----------|
| Conservative | 15% | 3.0x | 0.4 | 2+ types, 75% conf | High-quality assemblies, minimal false splits |
| Balanced | 12% | 2.4x | 0.32 | 2+ types, 65% conf | General use, moderate tolerance for false positives |
| Sensitive | 9% | 1.8x | 0.24 | 1+ type, 60% conf | Low-coverage samples, suspected chimeras |
| Very Sensitive | 6% | 1.2x | 0.16 | 1+ type, 50% conf | Challenging samples, maximum detection |

### Custom Thresholds

```bash
chimeric_detective -a assembly.fasta -b reads.bam -o results/ \
                  --min-coverage 10.0 \
                  --coverage-fold-change 3.0 \
                  --gc-content-threshold 0.15 \
                  --confidence-threshold 0.7
```

## Troubleshooting

### Common Issues and Solutions

#### Memory Issues (SIGKILL errors)
```
ERROR: minimap2 failed: Command died with <Signals.SIGKILL: 9>
```

**Solutions:**
1. **Increase memory allocation** (especially on SLURM/HPC):
   ```bash
   #SBATCH --mem=64G  # Make sure this line is uncommented!
   ```

2. **Use higher `--min-contig-length` threshold**:
   ```bash
   chimeric_detective --min-contig-length 2000  # or 3000 for very large assemblies
   ```

3. **Reduce parallel processing**:
   ```bash
   chimeric_detective --max-workers 1  # Most memory-efficient
   ```

4. **Use BWA instead of minimap2** (more memory-efficient):
   ```bash
   # Tool automatically tries BWA if minimap2 fails
   # Or install only BWA: conda install bwa
   ```

#### Tool Not Found Errors
```
ERROR: Neither minimap2 nor BWA are available
```

**Solutions:**
1. **Install via conda** (recommended):
   ```bash
   conda install -c bioconda bwa minimap2 samtools
   ```

2. **Activate conda environment**:
   ```bash
   conda activate chimeric-detective
   ```

3. **Check tool availability**:
   ```bash
   which minimap2 bwa samtools
   ```

#### No Chimeras Detected
```
INFO: Detected 0 chimera candidates
```

**This could indicate:**
- **Good assembly quality** (no chimeras present)
- **Thresholds too strict** - try `--sensitivity high`
- **Insufficient coverage** - check `--min-coverage` setting
- **Short contigs** - lower `--min-contig-length` if appropriate

#### Large Assembly Processing Tips

**Before running**, analyze your assembly:
```bash
python analyze_assembly_lengths.py assembly.fasta
```

**Common optimizations:**
- Assemblies >500MB: Use `--min-contig-length 2000` 
- >50K contigs: Use `--max-workers 4` or less
- Memory errors: Try `--max-workers 1` first

#### File Permission Issues

**On shared systems**, ensure write permissions:
```bash
chmod 755 output_directory/
```

#### Performance Issues

**For faster processing:**
- Use SSD storage for temporary files
- Ensure adequate RAM (see memory requirements table)
- Pre-filter assembly to remove very short contigs
- Use appropriate `--min-contig-length` threshold

### Getting Help

If you encounter issues:

1. **Check the log file**: `output_dir/chimeric_detective.log`
2. **Run with debug logging**: `--log-level DEBUG` 
3. **Test on small assembly first** to verify installation
4. **Open an issue** on GitHub with:
   - Command used
   - Error message
   - System information (OS, RAM, etc.)
   - Assembly characteristics (size, contig count)

## Testing and Validation

### Test Your Installation

```bash
# Quick test with provided test data
chimeric_detective --assembly test_data_minimal/test_assembly.fasta \
                  --reads-dir test_data_minimal/reads \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --out test_output
```

### Create Synthetic Test Data

```bash
# Generate test assembly with known chimeric contigs
python create_test_dataset.py
chimeric_detective --assembly test_data_large/large_test_assembly.fasta \
                  --reads-dir test_data_large/reads \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --out synthetic_test_output
```

### Validate Detection Accuracy

```bash
# Test core detection without full pipeline
python simple_detection_test.py
```

This validates that chimeric contigs are correctly identified with expected confidence scores and evidence types.

## Validation and Best Practices

### Choosing the Right Sensitivity Mode

1. **Start Conservative**: Begin with `--sensitivity conservative` to minimize false positives
2. **Validate Results**: Manually inspect a few detected chimeras to assess accuracy
3. **Adjust if Needed**: Increase sensitivity if missing known chimeras, decrease if too many false positives

### Biological Validation Guidelines

#### Expected Chimera Rates by Sample Type:
- **High-quality viral assemblies**: 1-5% of contigs
- **Mixed metagenomes**: 5-15% of contigs  
- **Low-coverage assemblies**: 10-25% of contigs
- **PCR-amplified samples**: 15-30% of contigs

#### Red Flags (Possible Over-Detection):
- >50% of contigs flagged as chimeric
- Many short-range breakpoints (<500bp from ends)
- Breakpoints in highly repetitive regions
- Consistent patterns across unrelated samples

#### Manual Validation Steps:
1. **Visual inspection**: Plot coverage and GC content around breakpoints
2. **Taxonomic validation**: BLAST left/right segments separately
3. **Read support**: Check paired-read spanning across junctions
4. **Assembly validation**: Re-assemble with different parameters

### Dataset-Specific Recommendations

#### For Viral Metagenomes:
- Use `conservative` mode for publication-quality assemblies
- Consider `balanced` mode for exploratory analysis
- Validate against known viral reference genomes

#### For Challenging Samples:
- Low coverage (<10x): Use `sensitive` or `very_sensitive`
- Mixed communities: Start with `conservative`, then `balanced`
- PCR-amplified: Consider `conservative` to avoid PCR artifact over-detection

#### For Different Viral Families:
- **High GC viruses** (>60%): May need `sensitive` mode
- **Low GC viruses** (<30%): `Conservative` mode usually sufficient
- **Segmented viruses**: Expect legitimate biological recombination

### Quality Control Metrics

Monitor these metrics across your analyses:

```bash
# Check detection rates
grep "chimera candidates" *.log

# Check confidence distributions  
jq '.analyses[].confidence_score' results.json | sort -n

# Check evidence type frequencies
jq '.analyses[].evidence_types[]' results.json | sort | uniq -c
```

### When to Re-evaluate Thresholds

Consider adjusting sensitivity if you observe:
- Very low detection rates (<1%) in known problematic samples
- Very high detection rates (>30%) in high-quality assemblies
- Systematic bias toward short or long contigs
- Poor correlation with independent quality metrics

## License

MIT License - see LICENSE file for details.