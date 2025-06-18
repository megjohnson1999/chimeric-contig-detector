# Read-Pair Based Chimera Detection

A focused, clean implementation for detecting chimeric contigs in viral metagenomic co-assemblies using only read-pair information. This approach is specifically designed for scenarios where coverage-based methods fail due to varying sample depths.

## Key Features

### 🎯 **Focused Detection Strategy**
- **Single method approach**: Uses only read-pair analysis for clean, debuggable results
- **Statistical rigor**: Pluggable statistical methods (MAD, IQR, z-score, KDE)
- **Baseline-aware**: Establishes proper baselines from well-behaved read pairs

### ⚙️ **Highly Configurable**
- **All parameters configurable**: Window sizes, thresholds, statistical methods
- **Multiple input formats**: Handles various BAM formats gracefully
- **Flexible output**: JSON, TSV, BED, VCF formats with compression support

### 🧩 **Modular Architecture**
- **Clean separation**: Data processing, analysis, and output completely separated
- **Extensible design**: Easy to add new statistical methods or output formats
- **Dependency injection**: Components are independently testable

### 📊 **Comprehensive Output**
- **Machine-readable formats**: JSON and TSV for automated pipelines
- **Genome browser compatible**: BED and VCF formats for manual inspection
- **Rich visualizations**: Interactive HTML reports with plotly

## Quick Start

### Command Line Usage

```bash
# Basic usage with defaults
chimeric_detective_readpair input.bam assembly.fasta

# Specify output directory and formats
chimeric_detective_readpair input.bam assembly.fasta \
    -o results/ \
    --output-formats json tsv bed

# Use custom parameters
chimeric_detective_readpair input.bam assembly.fasta \
    --window-size 2000 \
    --step-size 1000 \
    --min-mapping-quality 30 \
    --outlier-method mad
```

### Python API Usage

```python
from chimeric_detective.readpair_based import ReadPairChimeraDetector, DetectorConfig

# Basic usage
detector = ReadPairChimeraDetector.create_default()
candidates = detector.detect_chimeras("input.bam")
detector.write_results("output_dir/")

# Custom configuration
config = DetectorConfig()
config.detection.proper_pair_drop_threshold = 0.4
config.sliding_window.window_size = 2000
config.output.formats = ["json", "bed"]

detector = ReadPairChimeraDetector(config)
candidates = detector.detect_chimeras("input.bam")
```

## Detection Method

### 1. **Baseline Establishment**
- Analyzes properly paired reads across the entire contig
- Calculates insert size distribution statistics (median, MAD)
- Establishes expected proper pair ratios

### 2. **Sliding Window Analysis**
- Configurable window size and step size
- Calculates metrics for each window:
  - Proper pair ratio
  - Insert size distribution
  - Discordant pair ratio
  - Coverage uniformity

### 3. **Anomaly Detection**
- **Proper pair drops**: Significant decreases in proper pair ratio
- **Insert size shifts**: Statistical changes in insert size distribution
- **Discordant spikes**: Increases in discordant pair rates

### 4. **Breakpoint Identification**
- Clusters nearby anomalies
- Identifies most likely breakpoint positions
- Calculates confidence scores based on multiple evidence types

## Configuration Options

### Read Quality Filtering
```yaml
read_quality:
  min_mapping_quality: 20
  min_base_quality: 20
  require_proper_pairs: true
  max_edit_distance: null
  exclude_duplicates: true
  exclude_secondary: true
  exclude_supplementary: true
```

### Insert Size Analysis
```yaml
insert_size:
  outlier_method: "mad"  # mad, iqr, zscore, kde
  outlier_threshold: 3.0
  min_pairs_for_baseline: 1000
  max_insert_size: 10000
  binning_strategy: "auto"
```

### Sliding Window Parameters
```yaml
sliding_window:
  window_size: 1000
  step_size: 500
  min_pairs_per_window: 50
  window_type: "fixed"
  edge_behavior: "truncate"
```

### Detection Thresholds
```yaml
detection:
  proper_pair_drop_threshold: 0.5
  insert_size_shift_threshold: 2.0
  discordant_pair_threshold: 0.3
  min_anomaly_length: 100
  merge_distance: 500
  statistical_test: "ks"  # ks, anderson, chi2
```

## Output Formats

### JSON Output
```json
{
  "metadata": {
    "input_bam": "sample.bam",
    "analysis_time_seconds": 45.2,
    "total_candidates": 12
  },
  "summary": {
    "high_confidence": 3,
    "medium_confidence": 5,
    "low_confidence": 4
  },
  "breakpoints": [
    {
      "contig": "contig_001",
      "position": 15420,
      "confidence": 0.876,
      "anomaly_types": ["proper_pair_drop", "insert_size_shift"]
    }
  ]
}
```

### TSV Output
```
contig	position	confidence	anomaly_types	left_proper_ratio	right_proper_ratio
contig_001	15420	0.876	proper_pair_drop,insert_size_shift	0.891	0.234
contig_002	8934	0.743	discordant_spike	0.823	0.801
```

### BED Format (for genome browsers)
```
track name="Chimera_Breakpoints" description="Chimeric breakpoints"
contig_001	15320	15520	Breakpoint_1	876	+	15320	15520	255,0,0
```

## Advanced Usage

### Preset Configurations

```python
from chimeric_detective.readpair_based import DetectorFactory

# High sensitivity (more candidates, may include false positives)
detector = DetectorFactory.create_sensitive()

# High specificity (fewer candidates, high confidence)
detector = DetectorFactory.create_specific()

# Default balanced approach
detector = DetectorFactory.create_default()
```

### Configuration Files

```bash
# Generate configuration templates
chimeric_detective_readpair config sensitive -o sensitive_config.yaml
chimeric_detective_readpair config specific -o specific_config.yaml

# Use configuration file
chimeric_detective_readpair -c my_config.yaml input.bam assembly.fasta
```

### Validation and Debugging

```bash
# Validate BAM file compatibility
chimeric_detective_readpair validate input.bam

# Enable debug mode
chimeric_detective_readpair input.bam assembly.fasta --debug
```

## Performance Considerations

### Memory Usage
- Streams BAM data rather than loading everything into memory
- Configurable memory limits for large datasets
- Efficient caching of frequently accessed regions

### Processing Speed
- Parallelizable across contigs
- Configurable threading support
- Optional temporary directory for intermediate files

### Large Dataset Handling
```python
config = DetectorConfig()
config.threads = 8
config.memory_limit_gb = 16.0
config.temp_dir = "/tmp/chimera_detection"

detector = ReadPairChimeraDetector(config)
```

## Integration with Pipelines

### Nextflow Integration
```nextflow
process CHIMERA_DETECTION {
    input:
    tuple val(sample_id), path(bam), path(assembly)
    
    output:
    tuple val(sample_id), path("${sample_id}_chimeras.json")
    
    script:
    """
    chimeric_detective_readpair ${bam} ${assembly} \\
        -o . \\
        --output-formats json \\
        --output-prefix ${sample_id}_chimeras
    """
}
```

### Snakemake Integration
```python
rule detect_chimeras:
    input:
        bam="alignments/{sample}.bam",
        assembly="assemblies/{sample}.fasta"
    output:
        json="chimeras/{sample}.json",
        tsv="chimeras/{sample}.tsv"
    shell:
        """
        chimeric_detective_readpair {input.bam} {input.assembly} \\
            -o chimeras/ \\
            --output-formats json tsv \\
            --output-prefix {wildcards.sample}
        """
```

## Comparison with Other Approaches

| Feature | Read-Pair Based | Multi-Method | GC-Only |
|---------|----------------|--------------|---------|
| **Speed** | Fast | Slow | Very Fast |
| **Accuracy** | High for obvious chimeras | Variable | Moderate |
| **Specificity** | High | Variable | Low |
| **Dependencies** | BAM file only | BAM + assembly + reads | Assembly only |
| **Debugging** | Easy | Difficult | Very Easy |
| **Configuration** | Highly configurable | Complex | Simple |
| **Co-assembly Friendly** | Yes | No (coverage issues) | Yes |

## Troubleshooting

### Common Issues

**No candidates detected:**
- Check BAM file has paired reads: `chimeric_detective_readpair validate input.bam`
- Lower detection thresholds in configuration
- Ensure contigs are long enough (>2x window size)

**Too many false positives:**
- Increase `proper_pair_drop_threshold`
- Use `create_specific()` preset
- Increase `min_pairs_per_window`

**Poor performance:**
- Enable threading: `--threads 4`
- Increase window step size
- Set memory limit appropriately

### Debug Mode

```bash
chimeric_detective_readpair input.bam assembly.fasta --debug --verbose
```

This provides detailed logging of:
- BAM parsing statistics
- Window analysis metrics
- Statistical test results
- Anomaly detection details