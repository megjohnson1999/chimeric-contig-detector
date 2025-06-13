# Chimeric Detective Tutorial

This tutorial will guide you through using Chimeric Detective to detect and resolve chimeric contigs in viral metagenomic assemblies.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Input Preparation](#input-preparation)
4. [Basic Usage Examples](#basic-usage-examples)
5. [Advanced Usage](#advanced-usage)
6. [Understanding Results](#understanding-results)
7. [Troubleshooting](#troubleshooting)

## Installation

### From PyPI (Recommended)
```bash
pip install chimeric-detective
```

### From Source
```bash
git clone https://github.com/yourusername/chimeric-detective.git
cd chimeric-detective
pip install -e .
```

### External Dependencies

Chimeric Detective requires some external bioinformatics tools:

```bash
# Install via conda (recommended)
conda install -c bioconda bwa minimap2 samtools blast

# Or via package manager (Ubuntu/Debian)
sudo apt-get install bwa minimap2 samtools ncbi-blast+
```

## Quick Start

The fastest way to get started is with a BAM file containing pre-aligned reads:

```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --bam aligned_reads.bam \
                  --out results/
```

## Input Preparation

### Assembly File
Your assembly should be in FASTA format. Chimeric Detective works best with:
- Viral metagenomic assemblies
- Contigs ≥ 1000 bp (configurable)
- Good quality assemblies with reasonable coverage

### Read Data
You can provide reads in several formats:

1. **Pre-aligned BAM file** (fastest)
2. **Raw FASTQ files** (single or paired-end)
3. **Multiple samples in a directory**

#### Preparing Read Alignments
If you don't have a BAM file, Chimeric Detective can align reads for you:

```bash
# Paired-end reads
chimeric_detective --assembly assembly.fasta \
                  --reads1 sample_R1.fastq.gz \
                  --reads2 sample_R2.fastq.gz \
                  --out results/

# Single-end reads  
chimeric_detective --assembly assembly.fasta \
                  --reads sample.fastq.gz \
                  --out results/
```

## Basic Usage Examples

### Example 1: Basic Analysis with BAM File

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --bam aligned_reads.bam \
                  --out chimera_analysis/
```

This will:
- Detect chimeric contigs using default parameters
- Classify each chimera type
- Split technical artifacts
- Generate an interactive HTML report

### Example 2: Analysis with Raw Reads

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --reads1 sample_R1.fastq.gz \
                  --reads2 sample_R2.fastq.gz \
                  --out chimera_analysis/ \
                  --threads 8
```

### Example 3: Multiple Samples

```bash
# Analyze each sample separately (default)
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --out chimera_analysis/

# Merge all samples for combined analysis
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --multi-sample-mode merged \
                  --out chimera_analysis/

# Process with parallel workers
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --max-workers 8 \
                  --parallel \
                  --out chimera_analysis/
```

## Multi-Sample Processing

Chimeric Detective provides powerful capabilities for processing multiple samples efficiently. This is particularly useful for comparative studies or when you have reads from multiple sequencing runs.

### Processing Modes

#### 1. Separate Analysis (Default)
Each sample is analyzed independently, producing individual results:

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --multi-sample-mode separate \
                  --out results/
```

**Output structure:**
```
results/
├── sample_1/
│   ├── cleaned_assembly.fasta
│   ├── chimeric_detective_report.html
│   └── ...
├── sample_2/
│   ├── cleaned_assembly.fasta
│   ├── chimeric_detective_report.html
│   └── ...
├── multi_sample_report.html     # Combined overview
├── multi_sample_summary.tsv     # Summary table
└── combined_results.json        # All results
```

#### 2. Merged Analysis
All reads are combined before analysis, useful for increasing coverage:

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --multi-sample-mode merged \
                  --out results/
```

#### 3. Batch Processing
Memory-efficient processing for many samples:

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --multi-sample-mode batch \
                  --batch-size 5 \
                  --out results/
```

### File Naming Patterns

Chimeric Detective supports various file naming conventions:

| Pattern | Example Files | Description |
|---------|---------------|-------------|
| `*_R{1,2}.fastq.gz` | `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz` | Standard Illumina naming |
| `*_{1,2}.fq.gz` | `sample1_1.fq.gz`, `sample1_2.fq.gz` | Alternative numbering |
| `*.R{1,2}.fastq` | `sample1.R1.fastq`, `sample1.R2.fastq` | With dots |
| `*_F.fastq.gz` | `sample1_F.fastq.gz` | Single-end with suffix |

### Parallel Processing

```bash
# Use multiple CPU cores
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --max-workers 8 \
                  --parallel \
                  --out results/

# Disable parallel processing (sequential)
chimeric_detective --assembly viral_contigs.fasta \
                  --reads-dir /path/to/reads/ \
                  --no-parallel \
                  --out results/
```

### Multi-Sample Examples

#### Example 1: Comparative Study
You have samples from different conditions and want to compare chimera patterns:

```bash
# Directory structure:
# reads/
# ├── control_1_R1.fastq.gz
# ├── control_1_R2.fastq.gz
# ├── control_2_R1.fastq.gz
# ├── control_2_R2.fastq.gz
# ├── treatment_1_R1.fastq.gz
# ├── treatment_1_R2.fastq.gz
# └── treatment_2_R2.fastq.gz

chimeric_detective --assembly viral_assembly.fasta \
                  --reads-dir reads/ \
                  --reads-pattern "*_R{1,2}.fastq.gz" \
                  --multi-sample-mode separate \
                  --max-workers 4 \
                  --out comparative_analysis/
```

#### Example 2: Low Coverage Samples
Merge multiple low-coverage samples to improve detection:

```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads-dir low_coverage_reads/ \
                  --multi-sample-mode merged \
                  --min-coverage 2.0 \
                  --sensitivity high \
                  --out merged_analysis/
```

#### Example 3: Large-Scale Study
Process many samples efficiently:

```bash
chimeric_detective --assembly viral_assembly.fasta \
                  --reads-dir large_study_reads/ \
                  --multi-sample-mode batch \
                  --batch-size 10 \
                  --max-workers 16 \
                  --out large_scale_results/
```

### Multi-Sample Report Features

The multi-sample report provides:

1. **Cross-sample comparison** of chimera detection rates
2. **Interactive visualizations** showing patterns across samples
3. **Summary statistics** for each sample and overall
4. **Links to individual reports** for detailed analysis
5. **Downloadable summary tables** for further analysis

### Performance Considerations

- **Separate mode**: Best for comparing samples, requires more disk space
- **Merged mode**: Best for low-coverage samples, uses less compute time
- **Batch mode**: Best for memory-limited systems with many samples
- **Parallel processing**: Scales with available CPU cores
- **Max workers**: Set to number of CPU cores, but consider memory usage

### Troubleshooting Multi-Sample Issues

#### No samples found
```bash
# Check your pattern matches files
ls /path/to/reads/*_R{1,2}.fastq.gz

# Try different patterns
--reads-pattern "*_{1,2}.fq.gz"
--reads-pattern "*.R{1,2}.fastq"
```

#### Memory issues
```bash
# Use batch processing
--multi-sample-mode batch --batch-size 3

# Reduce parallel workers
--max-workers 2

# Filter short contigs
--min-contig-length 2000
```

#### Inconsistent results across samples
```bash
# Use consistent parameters
--min-coverage 5.0 --confidence-threshold 0.7

# Check sample quality first
--generate-report  # Review individual reports
```

## Advanced Usage

### Custom Parameters

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --bam aligned_reads.bam \
                  --out results/ \
                  --min-contig-length 500 \
                  --min-coverage 3.0 \
                  --coverage-fold-change 3.0 \
                  --confidence-threshold 0.7 \
                  --sensitivity high \
                  --threads 16
```

### With Reference Database

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --bam aligned_reads.bam \
                  --reference viral_reference_db.fasta \
                  --out results/
```

### Custom Splitting Behavior

```bash
chimeric_detective --assembly viral_contigs.fasta \
                  --bam aligned_reads.bam \
                  --out results/ \
                  --no-split-pcr \
                  --preserve-biological \
                  --min-split-length 1000
```

## Understanding Results

### Output Directory Structure

```
results/
├── cleaned_assembly.fasta          # Assembly with technical chimeras split
├── chimeric_contigs/               # Individual chimeric contig files
│   ├── contig_123_chimeric.fasta
│   ├── contig_123_split_A.fasta
│   └── contig_123_split_B.fasta
├── chimeric_detective_report.html  # Interactive HTML report
├── chimeric_detective_results.json # Machine-readable results
├── splitting_decisions.tsv         # Table of all modifications
├── figures/                        # Static visualizations
│   ├── chimera_types.png
│   └── confidence_distribution.png
└── chimeric_detective.log          # Detailed log file
```

### Key Output Files

1. **`cleaned_assembly.fasta`**: Your original assembly with technical chimeras split at breakpoints
2. **`chimeric_detective_report.html`**: Interactive report with visualizations and detailed explanations
3. **`splitting_decisions.tsv`**: Tab-separated table of all decisions made
4. **`chimeric_detective_results.json`**: Complete machine-readable results

### Understanding the HTML Report

The interactive HTML report contains:

1. **Summary Statistics**: Overview of detected chimeras and decisions
2. **Interactive Visualizations**: 
   - Chimera type distribution
   - Confidence score distributions  
   - Evidence type overview
3. **Detailed Results Table**: Complete analysis for each chimeric contig
4. **Individual Contig Plots**: Coverage, GC content, and evidence for each chimera

### Interpreting Chimera Types

- **Technical Artifact**: Assembly errors that should be split
- **PCR Chimera**: Amplification artifacts that should typically be split
- **Biological Recombination**: Genuine recombination events that should be preserved
- **Provirus Integration**: Virus-host integration sites
- **Horizontal Gene Transfer**: Gene transfer events

### Confidence Scores

- **High (>0.8)**: Very confident in classification
- **Medium (0.5-0.8)**: Moderately confident 
- **Low (<0.5)**: Low confidence, manual review recommended

## Parameter Tuning

### Sensitivity Levels

```bash
# High sensitivity (detect more chimeras, more false positives)
--sensitivity high

# Medium sensitivity (default, balanced)
--sensitivity medium  

# Low sensitivity (detect fewer chimeras, fewer false positives)
--sensitivity low
```

### Coverage Parameters

```bash
# Minimum coverage for analysis
--min-coverage 10.0

# Minimum fold change to consider
--coverage-fold-change 2.5
```

### Composition Parameters

```bash
# GC content difference threshold
--gc-content-threshold 0.15

# K-mer distance threshold  
--kmer-distance-threshold 0.4
```

## Common Workflows

### Workflow 1: Standard Viral Metagenome Analysis

```bash
# Step 1: Basic analysis
chimeric_detective --assembly viral_contigs.fasta \
                  --reads1 sample_R1.fastq.gz \
                  --reads2 sample_R2.fastq.gz \
                  --out initial_analysis/

# Step 2: Review results in HTML report

# Step 3: Re-run with adjusted parameters if needed
chimeric_detective --assembly viral_contigs.fasta \
                  --reads1 sample_R1.fastq.gz \
                  --reads2 sample_R2.fastq.gz \
                  --sensitivity high \
                  --confidence-threshold 0.7 \
                  --out refined_analysis/
```

### Workflow 2: High-Quality Assembly Analysis

```bash
# For high-quality assemblies, use stricter parameters
chimeric_detective --assembly high_quality_assembly.fasta \
                  --bam aligned_reads.bam \
                  --min-coverage 10.0 \
                  --coverage-fold-change 3.0 \
                  --confidence-threshold 0.8 \
                  --out results/
```

## Troubleshooting

### Common Issues

#### 1. No Chimeras Detected
- Check that your assembly has sufficient coverage
- Try lowering `--sensitivity` to `high`
- Reduce `--min-coverage` threshold
- Verify reads are properly aligned

#### 2. Too Many False Positives
- Increase `--confidence-threshold`
- Use `--sensitivity low`
- Increase `--coverage-fold-change` threshold

#### 3. External Tool Errors
```bash
# Check if tools are installed
which bwa
which minimap2  
which samtools

# Install missing tools
conda install -c bioconda bwa minimap2 samtools
```

#### 4. Memory Issues
- Use `--min-contig-length` to filter short contigs
- Process subsets of your assembly
- Increase available memory

### Getting Help

```bash
# View all available options
chimeric_detective --help

# Check version
chimeric_detective --version
```

### Log Files

Check the log file (`chimeric_detective.log`) for detailed information about:
- Processing steps
- Warning messages
- Error details

## Example Interpretation

Let's interpret a typical result:

```
Contig NODE_1234 shows clear evidence of being a technical chimera 
resulting from misassembly. At position 3,456, there is a sharp 8.3x 
drop in read coverage and a significant shift in GC content (58% to 42%). 
The 5' region (0-3,456bp) shows strong alignment to Caudovirales phages, 
while the 3' region (3,457-5,892bp) aligns to bacterial host sequences 
(Escherichia coli). The absence of paired reads spanning this junction 
further supports this being an assembly artifact rather than biological 
recombination. Confidence: HIGH.
```

**Interpretation**: 
- This is a clear technical artifact that should be split
- Evidence includes coverage drop, GC shift, taxonomic change, and lack of spanning reads
- High confidence means you can trust this classification

## Best Practices

1. **Always review the HTML report** before using cleaned assemblies
2. **Check log files** for warnings or errors
3. **Validate high-impact decisions** manually when possible
4. **Use reference databases** when available for better classification
5. **Keep intermediate files** (`--keep-intermediates`) for troubleshooting
6. **Document your parameters** for reproducibility

## Performance Tips

- Use pre-aligned BAM files when possible
- Use multiple threads (`--threads`)
- Filter short contigs (`--min-contig-length`)
- Use appropriate sensitivity settings for your data

For more detailed information, see the [full documentation](https://chimeric-detective.readthedocs.io).