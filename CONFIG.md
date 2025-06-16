# Configuration Guide

Chimeric Detective supports flexible configuration through YAML/JSON files, presets, environment variables, and CLI arguments. This guide covers all configuration options and how to use them effectively.

## Quick Start

### Using Presets

The fastest way to get started is with predefined presets:

```bash
# For large assemblies (>500MB, >25K contigs)
chimeric_detective --preset large -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/

# For high sensitivity detection
chimeric_detective --preset sensitive -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/

# For HPC/SLURM environments
chimeric_detective --preset hpc -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

### List Available Presets

```bash
chimeric_detective --list-presets
```

### Generate Configuration File

```bash
# Generate default configuration
chimeric_detective --generate-config my_config.yaml

# Generate configuration based on a preset
chimeric_detective --preset large --generate-config large_assembly_config.yaml
```

## Configuration File Format

Chimeric Detective supports both YAML (recommended) and JSON configuration files.

### YAML Example

```yaml
# Detection algorithm parameters
detection:
  min_contig_length: 1000          # Minimum contig length to analyze (bp)
  min_coverage: 5.0                # Minimum coverage threshold
  coverage_fold_change: 2.0        # Minimum coverage fold change to consider
  gc_content_threshold: 0.1        # Minimum GC content difference threshold
  kmer_distance_threshold: 0.3     # Minimum k-mer distance threshold
  confidence_threshold: 0.5        # Minimum confidence threshold for splitting
  min_split_length: 500            # Minimum length for split contigs (bp)

# Processing and performance settings
processing:
  threads: 1                       # Number of threads to use
  max_workers: 4                   # Maximum parallel workers for multi-sample processing
  parallel: true                   # Enable parallel processing of samples
  keep_intermediates: false        # Keep intermediate files
  batch_size: 5                    # Batch size for batch processing mode

# Output and reporting settings
output:
  generate_report: true            # Generate interactive HTML report
  log_level: "INFO"               # Logging level (DEBUG, INFO, WARNING, ERROR)

# Chimera resolution behavior
behavior:
  split_technical: true            # Split technical artifacts
  split_pcr: true                 # Split PCR chimeras
  preserve_biological: true       # Preserve biological recombination
  sensitivity: "medium"           # Detection sensitivity level (low, medium, high)

# Multi-sample processing settings
multi_sample:
  multi_sample_mode: "separate"   # How to process multiple samples (separate, merged, batch)
  reads_pattern: "*_R{1,2}.fastq.gz"  # Pattern for finding read files
```

### JSON Example

```json
{
  "detection": {
    "min_contig_length": 1000,
    "min_coverage": 5.0,
    "confidence_threshold": 0.5
  },
  "processing": {
    "max_workers": 8,
    "threads": 4
  },
  "behavior": {
    "sensitivity": "high"
  }
}
```

## Using Configuration Files

### Load Configuration File

```bash
chimeric_detective --config my_config.yaml -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

### Configuration File Locations

Chimeric Detective automatically searches for configuration files in these locations (in order):

1. Current directory:
   - `chimeric_detective.yaml`
   - `chimeric_detective.yml`
   - `chimeric_detective.json`

2. User config directory:
   - `~/.config/chimeric_detective/config.yaml`
   - `~/.config/chimeric_detective/config.yml`
   - `~/.config/chimeric_detective/config.json`

If found, the configuration is automatically loaded. You can override this with `--config`.

## Presets

### Available Presets

| Preset | Description | Use Case |
|--------|-------------|----------|
| `small` | Small assemblies (<100MB, <5K contigs) | Quick analysis, development |
| `large` | Large assemblies (>500MB, >25K contigs) | Memory-constrained environments |
| `sensitive` | High sensitivity detection | Comprehensive chimera detection |
| `conservative` | Low false positive rate | High-confidence results only |
| `hpc` | HPC/SLURM optimized | Limited memory environments |
| `development` | Development and debugging | Verbose logging, keep intermediates |

### Preset Details

#### Small Assembly Preset (`--preset small`)
```yaml
detection:
  min_contig_length: 1000
  min_coverage: 3.0
  confidence_threshold: 0.5
processing:
  max_workers: 8
  threads: 4
behavior:
  sensitivity: "medium"
```

#### Large Assembly Preset (`--preset large`)
```yaml
detection:
  min_contig_length: 2000          # Higher threshold for memory efficiency
  min_coverage: 5.0
  confidence_threshold: 0.6        # More conservative
processing:
  max_workers: 4                   # Limit workers
  threads: 2
  batch_size: 3                    # Smaller batches
behavior:
  sensitivity: "medium"
multi_sample:
  multi_sample_mode: "batch"       # Memory-efficient processing
```

#### Sensitive Preset (`--preset sensitive`)
```yaml
detection:
  coverage_fold_change: 1.4        # Lower thresholds = more sensitive
  gc_content_threshold: 0.07
  kmer_distance_threshold: 0.21
  confidence_threshold: 0.4
behavior:
  sensitivity: "high"
output:
  log_level: "DEBUG"              # Verbose logging
```

#### Conservative Preset (`--preset conservative`)
```yaml
detection:
  coverage_fold_change: 3.0        # Higher thresholds = fewer false positives
  gc_content_threshold: 0.15
  kmer_distance_threshold: 0.45
  confidence_threshold: 0.7
behavior:
  sensitivity: "low"
```

#### HPC Preset (`--preset hpc`)
```yaml
detection:
  min_contig_length: 1500
  min_coverage: 4.0
processing:
  max_workers: 2                   # Conservative for memory limits
  threads: 1
  batch_size: 2
  keep_intermediates: false        # Save space
```

## Environment Variables

Override configuration settings using environment variables with the `CHIMERIC_DETECTIVE_` prefix:

```bash
# Set minimum contig length
export CHIMERIC_DETECTIVE_MIN_CONTIG_LENGTH=2000

# Set maximum workers
export CHIMERIC_DETECTIVE_MAX_WORKERS=8

# Set log level
export CHIMERIC_DETECTIVE_LOG_LEVEL=DEBUG

# Set sensitivity
export CHIMERIC_DETECTIVE_SENSITIVITY=high

# Run with environment overrides
chimeric_detective -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

### Supported Environment Variables

| Variable | Type | Description |
|----------|------|-------------|
| `CHIMERIC_DETECTIVE_MIN_CONTIG_LENGTH` | int | Minimum contig length to analyze |
| `CHIMERIC_DETECTIVE_MIN_COVERAGE` | float | Minimum coverage threshold |
| `CHIMERIC_DETECTIVE_CONFIDENCE_THRESHOLD` | float | Confidence threshold for splitting |
| `CHIMERIC_DETECTIVE_MAX_WORKERS` | int | Maximum parallel workers |
| `CHIMERIC_DETECTIVE_THREADS` | int | Number of threads |
| `CHIMERIC_DETECTIVE_LOG_LEVEL` | str | Logging level (DEBUG, INFO, WARNING, ERROR) |
| `CHIMERIC_DETECTIVE_SENSITIVITY` | str | Sensitivity level (low, medium, high) |

## Configuration Precedence

Settings are applied in this order (later overrides earlier):

1. **Default values** - Built-in defaults
2. **Auto-loaded config file** - From default locations
3. **Preset** - Applied with `--preset`
4. **Config file** - Specified with `--config`
5. **Environment variables** - `CHIMERIC_DETECTIVE_*`
6. **CLI arguments** - Command line options (highest precedence)

### Example Precedence

```bash
# This command demonstrates precedence:
export CHIMERIC_DETECTIVE_MAX_WORKERS=16

chimeric_detective \
  --config my_config.yaml \      # Config file: max_workers = 8
  --preset large \               # Preset: max_workers = 4  
  --max-workers 12 \             # CLI: max_workers = 12 (wins!)
  -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/

# Final max_workers value: 12 (from CLI)
```

## Parameter Reference

### Detection Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `min_contig_length` | int | 1000 | ≥100 | Minimum contig length to analyze (bp) |
| `min_coverage` | float | 5.0 | ≥0 | Minimum coverage threshold |
| `coverage_fold_change` | float | 2.0 | >1 | Minimum coverage fold change to consider |
| `gc_content_threshold` | float | 0.1 | 0-1 | Minimum GC content difference threshold |
| `kmer_distance_threshold` | float | 0.3 | 0-1 | Minimum k-mer distance threshold |
| `confidence_threshold` | float | 0.5 | 0-1 | Minimum confidence threshold for splitting |
| `min_split_length` | int | 500 | ≥100 | Minimum length for split contigs (bp) |

### Processing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `threads` | int | 1 | Number of threads to use |
| `max_workers` | int | 4 | Maximum parallel workers for multi-sample processing |
| `parallel` | bool | true | Enable parallel processing of samples |
| `keep_intermediates` | bool | false | Keep intermediate files |
| `batch_size` | int | 5 | Batch size for batch processing mode |

### Output Parameters

| Parameter | Type | Default | Options | Description |
|-----------|------|---------|---------|-------------|
| `generate_report` | bool | true | - | Generate interactive HTML report |
| `log_level` | str | "INFO" | DEBUG, INFO, WARNING, ERROR | Logging level |

### Behavior Parameters

| Parameter | Type | Default | Options | Description |
|-----------|------|---------|---------|-------------|
| `split_technical` | bool | true | - | Split technical artifacts |
| `split_pcr` | bool | true | - | Split PCR chimeras |
| `preserve_biological` | bool | true | - | Preserve biological recombination |
| `sensitivity` | str | "medium" | low, medium, high | Detection sensitivity level |

### Multi-Sample Parameters

| Parameter | Type | Default | Options | Description |
|-----------|------|---------|---------|-------------|
| `multi_sample_mode` | str | "separate" | separate, merged, batch | How to process multiple samples |
| `reads_pattern` | str | "*_R{1,2}.fastq.gz" | - | Pattern for finding read files |

## Advanced Configuration Examples

### Memory-Optimized for Large Assemblies

```yaml
# config/large_assembly.yaml
detection:
  min_contig_length: 2000    # Skip short contigs
  min_coverage: 6.0          # Higher coverage requirement
  confidence_threshold: 0.6  # More conservative splitting

processing:
  max_workers: 2             # Limit parallelism
  threads: 1                 # Single-threaded per worker
  batch_size: 2              # Small batches
  keep_intermediates: false  # Save disk space

multi_sample:
  multi_sample_mode: "batch" # Memory-efficient mode
```

### High-Throughput Processing

```yaml
# config/high_throughput.yaml
detection:
  min_contig_length: 1500    # Moderate filtering
  confidence_threshold: 0.4  # More permissive

processing:
  max_workers: 16            # High parallelism
  threads: 2                 # Multi-threaded workers
  parallel: true
  batch_size: 10             # Large batches

output:
  generate_report: false     # Skip report generation for speed
  log_level: "WARNING"       # Minimal logging
```

### Development and Debugging

```yaml
# config/development.yaml
detection:
  min_contig_length: 500     # Include shorter contigs
  confidence_threshold: 0.3  # Lower threshold for testing

processing:
  max_workers: 1             # Single worker for debugging
  keep_intermediates: true   # Keep all intermediate files

output:
  log_level: "DEBUG"         # Verbose logging
  generate_report: true      # Always generate reports

behavior:
  sensitivity: "high"        # Catch everything for analysis
```

## Best Practices

### 1. Start with Presets

Always start with the most appropriate preset:

```bash
# For most users
chimeric_detective --preset large

# For detailed analysis
chimeric_detective --preset sensitive

# For production pipelines
chimeric_detective --preset conservative
```

### 2. Customize Incrementally

Generate a config file from a preset, then customize:

```bash
# Generate base config
chimeric_detective --preset large --generate-config my_config.yaml

# Edit my_config.yaml as needed
# Then use it
chimeric_detective --config my_config.yaml
```

### 3. Use Environment Variables for Deployment

Set common overrides via environment:

```bash
# In your deployment script
export CHIMERIC_DETECTIVE_MAX_WORKERS=8
export CHIMERIC_DETECTIVE_LOG_LEVEL=INFO

# Run multiple analyses with consistent settings
chimeric_detective --config project_config.yaml -a sample1.fasta ...
chimeric_detective --config project_config.yaml -a sample2.fasta ...
```

### 4. Memory Considerations

For large assemblies or limited memory:

1. **Increase `min_contig_length`** to 2000-3000
2. **Reduce `max_workers`** to 2-4
3. **Use `batch_size` of 2-3**
4. **Set `keep_intermediates: false`**

### 5. Performance Tuning

For faster processing:

1. **Increase `max_workers`** (but watch memory)
2. **Set `generate_report: false`** if not needed
3. **Use `log_level: WARNING`** to reduce I/O
4. **Consider `multi_sample_mode: merged`** for related samples

## Troubleshooting

### Configuration Not Loading

```bash
# Check if config file exists and is valid
chimeric_detective --config my_config.yaml --generate-config test.yaml
# If test.yaml is created successfully, your config is valid
```

### Memory Issues with Large Assemblies

```yaml
# Use these settings for very large assemblies
detection:
  min_contig_length: 3000
processing:
  max_workers: 1
  batch_size: 1
  keep_intermediates: false
```

### Performance Issues

```yaml
# For slow processing
processing:
  max_workers: 1           # Reduce parallelism
  threads: 1               # Single-threaded
output:
  generate_report: false   # Skip report generation
  log_level: "ERROR"       # Minimal logging
```

### Validation Errors

Common validation errors and solutions:

```bash
# Error: confidence_threshold must be between 0 and 1
# Fix: Set value between 0.1 and 0.9
confidence_threshold: 0.5

# Error: min_contig_length must be at least 100
# Fix: Set value ≥ 100
min_contig_length: 1000

# Error: coverage_fold_change must be greater than 1
# Fix: Set value > 1.0
coverage_fold_change: 2.0
```

## Migration from CLI-only Usage

If you're currently using CLI arguments, migrate to configuration files:

### Before (CLI only)
```bash
chimeric_detective \
  --min-contig-length 2000 \
  --max-workers 8 \
  --sensitivity high \
  --confidence-threshold 0.6 \
  -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

### After (Configuration file)
```yaml
# analysis_config.yaml
detection:
  min_contig_length: 2000
  confidence_threshold: 0.6
processing:
  max_workers: 8
behavior:
  sensitivity: "high"
```

```bash
chimeric_detective --config analysis_config.yaml -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
```

This approach is more maintainable, especially for complex workflows with multiple samples.