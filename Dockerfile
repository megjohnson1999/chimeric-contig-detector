# Multi-stage Dockerfile for Chimeric Detective
# Optimized for bioinformatics dependencies and small image size

FROM mambaorg/micromamba:1.5-focal as builder

# Install bioinformatics tools via conda/mamba for better dependency management
COPY environment.yml /tmp/environment.yml
USER root
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Production stage
FROM mambaorg/micromamba:1.5-focal

# Copy conda environment from builder
COPY --from=builder /opt/conda /opt/conda

# Install Python package
COPY . /tmp/chimeric-detective
WORKDIR /tmp/chimeric-detective

USER root
RUN /opt/conda/bin/pip install . && \
    rm -rf /tmp/chimeric-detective

# Create non-root user for security
RUN groupadd -r chimeric && useradd -r -g chimeric chimeric
USER chimeric

# Set up working directory
WORKDIR /data

# Expose common ports (if needed for web interface)
EXPOSE 8080

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH"
ENV PYTHONPATH="/opt/conda/lib/python3.9/site-packages"

# Verify installation
RUN chimeric_detective --help

# Default command
ENTRYPOINT ["chimeric_detective"]
CMD ["--help"]

# Labels for metadata
LABEL org.opencontainers.image.title="Chimeric Detective"
LABEL org.opencontainers.image.description="Tool for detecting and resolving chimeric contigs in viral metagenomic assemblies"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.authors="Meg Johnson <meganjohnson1w@gmail.com>"
LABEL org.opencontainers.image.url="https://github.com/megjohnson1999/chimeric-contig-detector"
LABEL org.opencontainers.image.source="https://github.com/megjohnson1999/chimeric-contig-detector"
LABEL org.opencontainers.image.licenses="MIT"