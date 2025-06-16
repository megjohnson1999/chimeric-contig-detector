# Docker Usage Guide

Chimeric Detective is available as a Docker container with all bioinformatics dependencies pre-installed.

## Quick Start

### Pull and Run
```bash
# Pull the latest image
docker pull ghcr.io/megjohnson1999/chimeric-contig-detector:latest

# Run with help
docker run --rm ghcr.io/megjohnson1999/chimeric-contig-detector:latest --help
```

### Basic Analysis
```bash
# Mount your data directory and run analysis
docker run --rm \
    -v $(pwd)/data:/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/assembly.fasta \
    --reads-dir /data/reads \
    --reads-pattern "*_R{1,2}.fastq.gz" \
    --out /data/results
```

## Available Tags

- `latest` - Latest stable release
- `v1.0.0`, `v1.1.0`, etc. - Specific version releases
- `main` - Latest development version

## Volume Mounts

The container uses `/data` as the working directory. Mount your data:

```bash
# Mount current directory
docker run --rm -v $(pwd):/data ghcr.io/megjohnson1999/chimeric-contig-detector:latest [options]

# Mount specific directories
docker run --rm \
    -v /path/to/assemblies:/assemblies:ro \
    -v /path/to/reads:/reads:ro \
    -v /path/to/output:/output \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /assemblies/my_assembly.fasta \
    --reads-dir /reads \
    --out /output
```

## Resource Limits

For large assemblies, set appropriate resource limits:

```bash
docker run --rm \
    --memory=16g \
    --cpus=8 \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/large_assembly.fasta \
    --reads-dir /data/reads \
    --max-workers 8 \
    --min-contig-length 2000 \
    --out /data/results
```

## Docker Compose

For more complex setups, use docker-compose:

```yaml
version: '3.8'
services:
  chimeric-detective:
    image: ghcr.io/megjohnson1999/chimeric-contig-detector:latest
    volumes:
      - ./data:/data
      - ./results:/results
    command: >
      --assembly /data/assembly.fasta
      --reads-dir /data/reads
      --out /results
    deploy:
      resources:
        limits:
          memory: 16G
          cpus: '8'
```

## Building Locally

To build the image yourself:

```bash
# Clone repository
git clone https://github.com/megjohnson1999/chimeric-contig-detector.git
cd chimeric-contig-detector

# Build image
docker build -t chimeric-detective .

# Run local build
docker run --rm -v $(pwd)/test_data:/data chimeric-detective --help
```

## Troubleshooting

### Permission Issues
If you encounter permission issues with output files:

```bash
# Run with your user ID
docker run --rm \
    --user $(id -u):$(id -g) \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest [options]
```

### Memory Issues
For large assemblies, increase memory limits:

```bash
# Increase Docker memory limit
docker run --rm \
    --memory=32g \
    --memory-swap=64g \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest [options]
```

### Tool Dependencies
All required bioinformatics tools (BWA, minimap2, samtools) are pre-installed in the container.

## Security Notes

- The container runs as a non-root user by default
- Only mount necessary directories
- Use read-only mounts for input data when possible

## Examples

### Single Sample
```bash
docker run --rm \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/sample1.fasta \
    --reads1 /data/sample1_R1.fastq.gz \
    --reads2 /data/sample1_R2.fastq.gz \
    --out /data/sample1_results
```

### Multiple Samples
```bash
docker run --rm \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/assembly.fasta \
    --reads-dir /data/reads \
    --reads-pattern "*_R{1,2}.fastq.gz" \
    --multi-sample-mode separate \
    --out /data/multi_sample_results
```

### Large Assembly Optimization
```bash
docker run --rm \
    --memory=64g \
    --cpus=4 \
    -v $(pwd):/data \
    ghcr.io/megjohnson1999/chimeric-contig-detector:latest \
    --assembly /data/large_assembly.fasta \
    --reads-dir /data/reads \
    --min-contig-length 2000 \
    --max-workers 4 \
    --out /data/results
```