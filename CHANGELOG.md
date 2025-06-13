# Changelog

All notable changes to Chimeric Detective will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of Chimeric Detective

## [1.0.0] - 2024-12-19

### Added
- **Core Detection Engine**
  - Multiple chimera detection methods (coverage, GC content, k-mer composition, read orientation)
  - Confidence scoring system for detection candidates
  - Support for both BAM files and raw FASTQ reads
  - Automatic read alignment using BWA or minimap2

- **Intelligent Classification System**
  - Machine learning-based chimera type classification
  - Support for technical artifacts, PCR chimeras, biological recombination, and more
  - Optional taxonomic classification using BLAST
  - Detailed evidence collection and analysis

- **Automated Resolution**
  - Intelligent splitting of technical artifacts at precise breakpoints
  - Preservation of biological recombination events
  - Configurable minimum split lengths and confidence thresholds
  - Comprehensive decision tracking and documentation

- **Interactive Visualizations**
  - HTML reports with embedded interactive plots
  - Coverage, GC content, and evidence visualizations
  - Individual contig analysis plots
  - Summary statistics and distribution charts
  - Responsive design for all modern browsers

- **Multi-Sample Processing**
  - Three processing modes: separate, merged, and batch
  - Parallel processing with configurable workers
  - Support for multiple file naming patterns
  - Cross-sample comparison reports
  - Memory-efficient batch processing for large studies

- **Command-Line Interface**
  - Intuitive Unix-style command structure
  - Comprehensive parameter configuration
  - Multiple input format support
  - Progress reporting and detailed logging
  - Flexible sensitivity settings

- **Documentation and Examples**
  - Comprehensive tutorial with real-world examples
  - API documentation with detailed docstrings
  - Example scripts for common use cases
  - Multi-sample processing examples
  - Troubleshooting guides

### Features
- **Detection Methods**
  - Coverage discontinuity analysis
  - GC content shift detection
  - K-mer composition change analysis
  - Read pair orientation inconsistency detection
  - Sequence homology analysis around breakpoints

- **Classification Types**
  - Technical assembly artifacts
  - PCR chimeras
  - Biological recombination events
  - Provirus integration sites
  - Horizontal gene transfer events

- **Output Formats**
  - Cleaned FASTA assemblies
  - Interactive HTML reports
  - Machine-readable JSON results
  - Tab-separated decision tables
  - Static PNG/SVG visualizations

- **Input Support**
  - FASTA assembly files
  - BAM/SAM alignment files
  - FASTQ/FASTA read files (single or paired-end)
  - Compressed files (gzip, bzip2)
  - Multiple file naming conventions
  - Directory-based multi-sample input

### Performance
- Efficient algorithms for large datasets
- Parallel processing support
- Memory-optimized data structures
- Typical processing time: <30 minutes for 5K contigs
- Memory usage: <8GB for most datasets

### Dependencies
- **Core**: Python 3.8+, Biopython, NumPy, pandas, scikit-learn
- **Visualization**: Plotly, matplotlib, seaborn
- **Interface**: Click, tqdm, Jinja2
- **External**: BWA/minimap2, samtools, BLAST (optional)

### Supported Platforms
- Linux (primary)
- macOS
- Windows (via WSL recommended)

## [0.9.0] - 2024-12-15 (Beta)

### Added
- Initial beta release for testing
- Core detection and analysis functionality
- Basic command-line interface
- Preliminary documentation

### Known Issues
- Limited multi-sample support
- Basic visualization capabilities
- Performance optimization needed

## Project Milestones

### Version 1.0.0 Goals ✅
- [x] Complete chimera detection pipeline
- [x] Intelligent classification system
- [x] Interactive HTML reports
- [x] Multi-sample processing
- [x] Comprehensive documentation
- [x] Example data and tutorials

### Future Releases (Planned)

#### Version 1.1.0
- [ ] Machine learning model improvements
- [ ] Additional visualization options
- [ ] Performance optimizations
- [ ] Extended file format support
- [ ] Integration with popular assembly tools

#### Version 1.2.0
- [ ] Advanced statistical analysis
- [ ] Comparative genomics features
- [ ] Cloud computing support
- [ ] Web interface option
- [ ] Extended database support

#### Version 2.0.0
- [ ] Real-time analysis capabilities
- [ ] Advanced machine learning models
- [ ] Integration with sequencing platforms
- [ ] Collaborative features
- [ ] Extended API for integration

## Notes

### Breaking Changes
- None in version 1.0.0 (initial release)

### Deprecations
- None in version 1.0.0

### Security
- All external tool calls are validated
- Input sanitization implemented
- No known security vulnerabilities

### Contributors
- Chimeric Detective Development Team
- Community contributors (see CONTRIBUTORS.md)

---

For detailed information about any release, please see the corresponding GitHub release notes and documentation.