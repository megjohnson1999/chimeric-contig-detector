# Chimeric Detective Benchmarking Suite

This directory contains tools for benchmarking Chimeric Detective against other chimera detection and assembly quality assessment tools.

## Quick Start

### Prerequisites

Install comparison tools using conda/mamba:

```bash
# Install conda/mamba if not available
conda install -c conda-forge mamba

# Install bioinformatics tools
mamba install -c bioconda checkv virsorter=2 quast bwa samtools pandas
```

### Quick Benchmark

Run a simple comparison with minimal setup:

```bash
# Basic usage with paired-end reads
./quick_benchmark.sh -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz

# With single-end reads
./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz

# With existing BAM file
./quick_benchmark.sh -a assembly.fasta -b aligned_reads.bam

# Custom output directory and threads
./quick_benchmark.sh -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o my_benchmark -t 8
```

### Comprehensive Benchmark

Use the Python script for more detailed analysis:

```bash
# With synthetic test data
python benchmark_chimeric_detective.py --create-test-data -o comprehensive_benchmark

# With real data
python benchmark_chimeric_detective.py -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results

# Run specific tools only
python benchmark_chimeric_detective.py -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz --tools chimeric_detective checkv virsorter2
```

## Tools Compared

### Primary Chimera Detection Tools

1. **Chimeric Detective** (our tool)
   - Multi-method detection (coverage, GC content, k-mer composition)
   - Machine learning classification for viral contigs
   - Automated resolution and cleaning
   - Interactive HTML reports

2. **CheckV**
   - Specialized for viral sequences
   - Detects host contamination in viral contigs
   - Quality assessment for viral genomes

3. **VirSorter2**
   - Viral contig identification and quality assessment
   - Can identify chimeric/contaminated sequences

### Assembly Quality Assessment Tools

4. **QUAST**
   - General assembly quality metrics
   - Misassembly detection for contigs
   - Comprehensive contig statistics

5. **metaQUAST**
   - Metagenomic-specific version of QUAST
   - Reference-free contig analysis
   - Better for complex metagenomic assemblies

6. **Bandage**
   - Assembly graph visualization
   - Can identify problematic contig connections
   - Useful for manual chimera inspection

## Output Files

### Quick Benchmark
- `benchmark_summary.tsv` - Tabular summary of all results
- `benchmark_report.md` - Detailed markdown report
- Tool-specific directories with individual outputs

### Comprehensive Benchmark
- `benchmark_summary.tsv` - Summary table
- `benchmark_report.md` - Detailed report
- `benchmark.log` - Execution log
- Individual tool output directories

## Interpreting Results

### Key Metrics

1. **Chimeras Detected**: Number of potential chimeric sequences identified
2. **Runtime**: Execution time in seconds
3. **Success Rate**: Whether the tool completed without errors

### Important Considerations

- **Higher counts ≠ better detection**: Different tools have different sensitivity/specificity trade-offs
- **False positives**: Some tools may flag legitimate biological sequences as chimeric
- **Tool-specific definitions**: Each tool may define "chimera" slightly differently
- **Dataset dependency**: Performance varies significantly with data type and quality

### Validation Approaches

1. **Manual inspection**: Review individual chimeric sequences
2. **Cross-validation**: Look for consensus across multiple tools
3. **Known controls**: Use datasets with known chimeric content
4. **Biological validation**: Check if detected chimeras make biological sense

## Example Interpretation

```
Tool                Status   Runtime(s)  Chimeras_Detected  Notes
Chimeric_Detective  SUCCESS  45          12                 Comprehensive analysis
VSEARCH_UCHIME      SUCCESS  8           18                 Fast, may have false positives
CheckV              SUCCESS  120         3                  Conservative, viral-specific
QUAST               SUCCESS  15          5                  General misassemblies
```

**Interpretation**: 
- VSEARCH detected the most chimeras but might include false positives
- CheckV was most conservative, focusing on clear contamination
- Chimeric Detective provided middle-ground detection with detailed analysis
- QUAST detected general assembly issues beyond just chimeras

## Troubleshooting

### Common Issues

1. **Tool not found errors**
   ```bash
   # Install missing tools
   mamba install -c bioconda [tool_name]
   ```

2. **Memory issues**
   ```bash
   # Reduce thread count
   ./quick_benchmark.sh -a assembly.fasta -1 reads.fastq.gz -t 2
   ```

3. **BAM file creation fails**
   ```bash
   # Check if BWA and samtools are installed
   conda install -c bioconda bwa samtools
   ```

4. **Python dependencies missing**
   ```bash
   pip install pandas matplotlib seaborn
   ```

### Performance Tips

- Use fewer threads if running out of memory
- Pre-create BAM files for multiple runs
- Filter short contigs to reduce processing time
- Use subsampled datasets for initial testing

## Creating Test Datasets

### Synthetic Data

```bash
# Create synthetic test data with known chimeras
python benchmark_chimeric_detective.py --create-test-data --test-contigs 50 --chimera-fraction 0.2 -o synthetic_test
```

### Real Data Validation

1. **Use published datasets** with known chimeric content
2. **Manual curation** of a subset for ground truth
3. **Spike-in controls** with artificially created chimeras
4. **Cross-reference** with literature reports

## Citation

If you use this benchmarking suite in your research, please cite:

```
Chimeric Detective: A comprehensive tool for detecting and resolving 
chimeric contigs in viral metagenomic assemblies
[Your citation information here]
```

## Contributing

To add new tools to the benchmark:

1. Add tool configuration to `benchmark_chimeric_detective.py`
2. Implement parsing function for tool output
3. Add installation instructions to this README
4. Test with various datasets

## Support

For issues or questions:
- GitHub Issues: https://github.com/megjohnson1999/chimeric-contig-detector/issues
- Documentation: See main project README