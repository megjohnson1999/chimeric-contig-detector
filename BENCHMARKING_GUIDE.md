# Chimeric Detective Benchmarking Guide

This guide walks you through evaluating Chimeric Detective's performance against other contig-level quality assessment tools.

## 🚀 Quick Start (5 minutes)

### 1. Install Comparison Tools

```bash
# Install tools via conda (recommended)
conda install -c bioconda checkv virsorter=2 quast bwa samtools

# Or install individually if needed
conda install -c bioconda checkv
conda install -c bioconda virsorter=2  
conda install -c bioconda quast
conda install -c bioconda bwa samtools
```

### 2. Run Quick Benchmark

```bash
cd benchmarking/
./quick_benchmark.sh -a your_assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz
```

**Output**: Results in timestamped directory with summary table and detailed report.

### 3. View Results

```bash
# Check summary
cat quick_benchmark_*/benchmark_summary.tsv

# View detailed report  
open quick_benchmark_*/benchmark_report.md  # macOS
xdg-open quick_benchmark_*/benchmark_report.md  # Linux
```

## 📊 Comprehensive Benchmarking

### Create Synthetic Test Data

Perfect for validating detection accuracy:

```bash
python benchmarking/benchmark_chimeric_detective.py \
    --create-test-data \
    --test-contigs 50 \
    --chimera-fraction 0.3 \
    -o synthetic_benchmark
```

This creates:
- 50 contigs total
- 30% are artificially chimeric (known ground truth)
- Documented breakpoint locations

### Benchmark Real Data

```bash
python benchmarking/benchmark_chimeric_detective.py \
    -a viral_assembly.fasta \
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -o real_data_benchmark \
    -t 8
```

### Run Specific Tools Only

```bash
# Compare only against CheckV and QUAST
python benchmarking/benchmark_chimeric_detective.py \
    -a assembly.fasta -1 reads.fastq.gz \
    --tools chimeric_detective checkv quast \
    -o targeted_benchmark
```

## 🔍 Understanding Results

### Summary Table Columns

| Column | Description |
|--------|-------------|
| **Tool** | Software name |
| **Status** | SUCCESS/FAILED |
| **Runtime(s)** | Execution time in seconds |
| **Chimeras_Detected** | Number of problematic contigs found |
| **Notes** | Tool-specific interpretation |

### Tool-Specific Interpretations

#### Chimeric Detective
- **What it detects**: Multi-method chimera identification
- **Output metric**: True chimeric contigs with breakpoint analysis
- **Strength**: Comprehensive analysis with resolution strategies

#### CheckV  
- **What it detects**: Host gene contamination in viral contigs
- **Output metric**: Contigs with host_genes > 0
- **Strength**: Conservative, viral-specific quality control

#### VirSorter2
- **What it detects**: Low-confidence viral sequences
- **Output metric**: Sequences with viral scores < 0.7
- **Strength**: Identifies ambiguous viral content

#### QUAST
- **What it detects**: General assembly misassemblies
- **Output metric**: Structural assembly problems
- **Strength**: Broad assembly quality assessment

### Example Result Interpretation

```
Tool                Status   Runtime(s)  Chimeras_Detected  Notes
Chimeric_Detective  SUCCESS  45          12                 Comprehensive analysis
CheckV              SUCCESS  120         3                  Conservative detection
VirSorter2          SUCCESS  90          8                  Low-confidence sequences  
QUAST               SUCCESS  15          5                  General misassemblies
```

**Analysis**:
- **Chimeric Detective** found 12 chimeras using multiple detection methods
- **CheckV** found 3 with clear host contamination (high confidence)
- **VirSorter2** found 8 ambiguous viral sequences (may include chimeras)
- **QUAST** found 5 general assembly problems

**Validation approaches**:
1. **Intersection analysis**: Which chimeras are detected by multiple tools?
2. **Manual inspection**: Review specific sequences flagged by each tool
3. **Biological assessment**: Do detected chimeras make biological sense?

## 🎯 Best Practices

### For Method Development

1. **Use synthetic data first** to validate detection accuracy
2. **Test multiple sensitivity levels** to understand trade-offs
3. **Manual validation** of a subset for ground truth
4. **Cross-tool comparison** to identify consensus chimeras

### For Production Use

1. **Start with default settings** and review results
2. **Compare with CheckV** for viral contamination baseline
3. **Manual inspection** of high-confidence detections
4. **Consider biological context** when interpreting results

### Performance Optimization

```bash
# For large datasets
./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz -t 16

# Memory-efficient processing
python benchmark_chimeric_detective.py \
    -a assembly.fasta -1 reads.fastq.gz \
    --tools chimeric_detective checkv \
    -t 4
```

## 🔧 Troubleshooting

### Common Issues

#### "Tool not found" errors
```bash
# Check what's installed
conda list | grep -E "(checkv|virsorter|quast|bwa|samtools)"

# Install missing tools
conda install -c bioconda [missing_tool]
```

#### Memory issues
```bash
# Reduce thread count
./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz -t 2

# Filter short contigs
chimeric_detective -a assembly.fasta -1 reads.fastq.gz --min-contig-length 2000
```

#### BAM file creation fails
```bash
# Check BWA installation
bwa
samtools --version

# Pre-create BAM file
bwa index assembly.fasta
bwa mem -t 8 assembly.fasta reads_R1.fastq.gz reads_R2.fastq.gz | \
    samtools sort -@ 8 -o aligned.bam -
samtools index aligned.bam

# Use in benchmark
./quick_benchmark.sh -a assembly.fasta -b aligned.bam
```

### Getting Help

1. **Check tool-specific logs** in output directories
2. **Review stderr.log files** for detailed error messages
3. **Test with smaller datasets** to isolate issues
4. **Consult individual tool documentation** for specific parameters

## 📈 Validation Strategies

### For Published Datasets

```bash
# Download published viral metagenome
# Run benchmark
python benchmark_chimeric_detective.py -a published_assembly.fasta -1 reads.fastq.gz

# Compare with literature findings
# Document detection rates and novel findings
```

### For Mock Communities

```bash
# Create known chimeric sequences
# Mix with clean viral genomes  
# Run benchmark and validate against known truth
```

### Cross-Validation

```bash
# Run on multiple related samples
for sample in sample1 sample2 sample3; do
    ./quick_benchmark.sh -a ${sample}_assembly.fasta -1 ${sample}_reads.fastq.gz
done

# Compare consistency across samples
```

## 📚 Further Reading

- [Benchmarking README](benchmarking/README.md) - Detailed technical documentation
- [CheckV Paper](https://www.nature.com/articles/s41587-020-00774-7) - Viral quality assessment
- [VirSorter2 Paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y) - Viral identification
- [QUAST Paper](https://academic.oup.com/bioinformatics/article/29/8/1072/228832) - Assembly quality metrics

---

*For technical support or feature requests, please open an issue on [GitHub](https://github.com/megjohnson1999/chimeric-contig-detector/issues).*