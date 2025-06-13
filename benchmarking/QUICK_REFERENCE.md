# Benchmarking Quick Reference

## 🚀 One-Minute Setup

```bash
# Install comparison tools
conda install -c bioconda checkv virsorter=2 quast bwa samtools

# Run benchmark  
cd benchmarking/
./quick_benchmark.sh -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz
```

## 📊 Commands Cheat Sheet

| Task | Command |
|------|---------|
| **Quick benchmark** | `./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz` |
| **With BAM file** | `./quick_benchmark.sh -a assembly.fasta -b aligned.bam` |
| **Custom threads** | `./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz -t 8` |
| **Synthetic data** | `python benchmark_chimeric_detective.py --create-test-data` |
| **Specific tools** | `python benchmark_chimeric_detective.py --tools chimeric_detective checkv` |

## 🔍 Understanding Output

| Tool | Detects | Good For |
|------|---------|----------|
| **Chimeric Detective** | Multi-method chimeras | Comprehensive detection |
| **CheckV** | Host contamination | Conservative validation |
| **VirSorter2** | Low-confidence viral | Ambiguous sequences |
| **QUAST** | Assembly errors | General quality |

## ⚡ Common Usage Patterns

### Development & Testing
```bash
# Test with synthetic data
python benchmark_chimeric_detective.py --create-test-data --test-contigs 20

# Quick validation
./quick_benchmark.sh -a test_assembly.fasta -1 test_reads.fastq.gz -t 2
```

### Production Analysis
```bash
# Full comparison
./quick_benchmark.sh -a real_assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -t 16

# Focus on viral quality
python benchmark_chimeric_detective.py --tools chimeric_detective checkv virsorter2
```

## 🎯 Result Interpretation

### Good Results
```
Tool                Status   Chimeras_Detected
Chimeric_Detective  SUCCESS  8
CheckV              SUCCESS  2  
VirSorter2          SUCCESS  5
QUAST               SUCCESS  3
```
**→ Reasonable detection with conservative CheckV baseline**

### Concerning Results  
```
Tool                Status   Chimeras_Detected
Chimeric_Detective  SUCCESS  45
CheckV              SUCCESS  1
VirSorter2          SUCCESS  3  
QUAST               SUCCESS  2
```
**→ High detection rate - check sensitivity settings or validate manually**

## 🔧 Troubleshooting

| Problem | Solution |
|---------|----------|
| "Tool not found" | `conda install -c bioconda [tool_name]` |
| Memory error | Use `-t 2` (fewer threads) |
| BAM creation fails | Pre-create: `bwa mem ... \| samtools sort` |
| No chimeras detected | Try `--sensitivity high` |

## 📁 Output Files

```
benchmark_results/
├── benchmark_summary.tsv     # Key metrics table
├── benchmark_report.md       # Detailed analysis  
├── chimeric_detective_results/ # Our tool output
├── checkv_results/           # CheckV output
├── virsorter2_results/       # VirSorter2 output
└── quast_results/            # QUAST output
```

---
**💡 Pro Tip**: Start with quick benchmark, then drill down with comprehensive analysis for detailed validation.