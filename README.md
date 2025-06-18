# Chimeric Detective

[![CI](https://github.com/megjohnson1999/chimeric-contig-detector/workflows/CI/badge.svg)](https://github.com/megjohnson1999/chimeric-contig-detector/actions)
[![codecov](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector/branch/main/graph/badge.svg)](https://codecov.io/gh/megjohnson1999/chimeric-contig-detector)
[![PyPI version](https://badge.fury.io/py/chimeric-detective.svg)](https://badge.fury.io/py/chimeric-detective)
[![Python](https://img.shields.io/pypi/pyversions/chimeric-detective.svg)](https://pypi.org/project/chimeric-detective/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Focused chimera detection tools for viral metagenomic assemblies**

Chimeric Detective provides two specialized approaches for detecting chimeric contigs, each optimized for different use cases and data availability.

## 🎯 Two Focused Approaches

### 1. **GC-Content Based Detection** (`chimeric_detective_simple`)
**Best for: Quick screening, assembly-only workflows**

- ✅ **Assembly file only** - No alignment required
- ✅ **Fast processing** - Analyze large assemblies quickly  
- ✅ **Simple and reliable** - Easy to understand and debug
- ✅ **GC shift detection** - Identifies clear compositional changes

```bash
# Quick GC-based detection
chimeric_detective_simple assembly.fasta -o results/
```

### 2. **Read-Pair Based Detection** (`chimeric_detective_readpair`)
**Best for: Detailed analysis, co-assembled viral datasets**

- 🎯 **Statistically rigorous** - Multiple configurable detection methods
- 🎯 **Co-assembly optimized** - Works with variable coverage depths
- 🎯 **Highly configurable** - Tune parameters for your data
- 🎯 **Rich output formats** - JSON, TSV, BED, VCF support

```bash
# Comprehensive read-pair analysis  
chimeric_detective_readpair alignment.bam assembly.fasta -o results/
```

## Quick Start

### Installation

```bash
# Install from PyPI
pip install chimeric-detective

# Or install development version
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector
pip install -e .
```

### Basic Usage

```bash
# GC-based detection (assembly only)
chimeric_detective_simple assembly.fasta

# Read-pair detection (requires BAM file)
chimeric_detective_readpair alignment.bam assembly.fasta

# Custom output directory and formats
chimeric_detective_readpair alignment.bam assembly.fasta \\
    -o my_results/ \\
    --output-formats json tsv bed
```

## When to Use Which Approach

| Use Case | Recommended Approach | Reason |
|----------|---------------------|---------|
| **Initial assembly screening** | GC-based | Fast, no alignment needed |
| **High-throughput pipelines** | GC-based | Minimal computational requirements |
| **Co-assembled viral datasets** | Read-pair | Handles variable coverage well |
| **Detailed breakpoint analysis** | Read-pair | Statistical rigor and rich output |
| **Publication-quality results** | Read-pair | Comprehensive evidence and visualization |
| **Limited computational resources** | GC-based | Lower memory and CPU requirements |

## Key Features

### GC-Based Detection
- **Sliding window analysis** of GC content changes
- **Adaptive window sizing** based on sequence length
- **Confidence scoring** for detected breakpoints
- **HTML visualization** with interactive plots
- **Minimal dependencies** and fast processing

### Read-Pair Detection  
- **Baseline establishment** from proper read pairs
- **Statistical anomaly detection** (MAD, IQR, z-score, KDE)
- **Multiple evidence types**: proper pair drops, insert size shifts, discordant spikes
- **Configurable parameters** via CLI or config files
- **Multiple output formats** for different downstream uses
- **Interactive visualizations** for manual inspection

## Configuration Examples

### GC-Based Detection
```bash
chimeric_detective_simple assembly.fasta \\
    --gc-threshold 0.15 \\
    --window-size 2000 \\
    --step-size 1000 \\
    --min-confidence 0.7
```

### Read-Pair Detection
```bash
chimeric_detective_readpair alignment.bam assembly.fasta \\
    --window-size 1000 \\
    --min-mapping-quality 30 \\
    --outlier-method mad \\
    --proper-pair-drop-threshold 0.5
```

### Using Configuration Files
```bash
# Generate configuration templates
chimeric_detective_readpair config sensitive -o sensitive.yaml
chimeric_detective_readpair config specific -o specific.yaml

# Use configuration file
chimeric_detective_readpair -c sensitive.yaml alignment.bam assembly.fasta
```

## Output Examples

### GC-Based Results
```
Contig: scaffold_123
Breakpoint: 15420
GC Left: 0.234
GC Right: 0.678  
GC Difference: 0.444
Confidence: 0.876
```

### Read-Pair Results (JSON)
```json
{
  "breakpoints": [
    {
      "contig": "scaffold_123",
      "position": 15420,
      "confidence": 0.876,
      "anomaly_types": ["proper_pair_drop", "insert_size_shift"],
      "left_proper_ratio": 0.891,
      "right_proper_ratio": 0.234
    }
  ]
}
```

## Advanced Usage

### Python API

```python
# GC-based detection
from chimeric_detective import SimpleChimeraDetector
detector = SimpleChimeraDetector(gc_content_threshold=0.1)
candidates = detector.detect_chimeras("assembly.fasta")

# Read-pair detection  
from chimeric_detective.readpair_based import ReadPairChimeraDetector, DetectorConfig
config = DetectorConfig()
config.detection.proper_pair_drop_threshold = 0.4
detector = ReadPairChimeraDetector(config)
candidates = detector.detect_chimeras("alignment.bam")
```

### Validation and Debugging

```bash
# Validate BAM file for read-pair analysis
chimeric_detective_readpair validate alignment.bam

# Enable debug output
chimeric_detective_readpair alignment.bam assembly.fasta --debug --verbose
```

## Documentation

- **[GC-Based Detection Guide](README_SIMPLE.md)** - Detailed guide for assembly-only detection
- **[Read-Pair Detection Guide](README_READPAIR.md)** - Comprehensive read-pair analysis documentation
- **[Examples](examples/)** - Usage examples and sample configurations

## Performance

| Dataset Size | GC-Based Time | Read-Pair Time | Memory Usage |
|--------------|---------------|----------------|--------------|
| Small (100 contigs) | <1 min | 2-5 min | <1 GB |
| Medium (1K contigs) | 2-5 min | 10-20 min | 2-4 GB |
| Large (10K contigs) | 10-20 min | 1-2 hours | 4-8 GB |

## Citation

If you use Chimeric Detective in your research, please cite:

```
Chimeric Detective: Focused tools for chimeric contig detection in viral metagenomes
Johnson et al. (2024)
```

## Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/megjohnson1999/chimeric-contig-detector/issues)
- **Documentation**: [Full Documentation](https://chimeric-detective.readthedocs.io/)
- **Examples**: [Example Scripts](examples/)