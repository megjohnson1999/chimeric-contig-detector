# Getting Started with Chimeric Detective Development

This guide will help you set up Chimeric Detective for development and contributing.

## Quick Setup for Contributors

### 1. Fork and Clone
```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOURUSERNAME/chimeric-contig-detector.git
cd chimeric-contig-detector
```

### 2. Set up Development Environment
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e .[dev,test,docs]

# Install external bioinformatics tools
# Ubuntu/Debian:
sudo apt-get install bwa minimap2 samtools ncbi-blast+

# macOS:
brew install bwa minimap2 samtools blast

# Or via conda:
conda install -c bioconda bwa minimap2 samtools blast
```

### 3. Verify Installation
```bash
# Run tests
pytest

# Check code style
black --check src/
flake8 src/

# Test CLI
chimeric_detective --help
```

### 4. Run Examples
```bash
# Basic example
cd examples
python example_usage.py

# Multi-sample example
python multi_sample_example.py
```

## Development Workflow

### Making Changes
1. Create a new branch: `git checkout -b feature/your-feature-name`
2. Make your changes
3. Run tests: `pytest`
4. Check code style: `black src/ && flake8 src/`
5. Commit changes: `git commit -m "Description of changes"`
6. Push branch: `git push origin feature/your-feature-name`
7. Create pull request on GitHub

### Code Style
- Use Black for formatting: `black src/`
- Follow PEP 8 with 100 character line limit
- Add type hints where appropriate
- Write docstrings for public functions

### Testing
- Write tests for new features
- Run `pytest` before committing
- Aim for good test coverage
- Include integration tests for CLI features

## Project Structure

```
chimeric-contig-detector/
├── src/chimeric_detective/     # Main package
│   ├── __init__.py
│   ├── cli.py                  # Command-line interface
│   ├── detector.py             # Chimera detection
│   ├── analyzer.py             # Classification and analysis
│   ├── resolver.py             # Assembly cleaning
│   ├── visualizer.py           # Reporting and visualization
│   ├── multi_sample.py         # Multi-sample processing
│   ├── utils.py                # Utility functions
│   └── templates/              # HTML templates
├── tests/                      # Test suite
├── examples/                   # Example scripts
├── docs/                       # Documentation
├── .github/                    # GitHub templates and workflows
└── setup files                # Configuration files
```

## Common Development Tasks

### Adding a New Detection Method
1. Add method to `detector.py`
2. Add tests in `tests/test_detector.py`
3. Update documentation
4. Add example usage

### Adding a New Chimera Type
1. Update classification in `analyzer.py`
2. Add explanation templates
3. Update visualization code
4. Add tests

### Improving Visualizations
1. Update `visualizer.py`
2. Test with various data types
3. Update templates if needed
4. Test HTML output

## Release Process

### Version Numbering
- Follow semantic versioning (MAJOR.MINOR.PATCH)
- Update version in `pyproject.toml`
- Update `CHANGELOG.md`

### Creating a Release
1. Update version numbers
2. Update CHANGELOG.md
3. Run full test suite
4. Create release on GitHub
5. Publish to PyPI (maintainers only)

## Getting Help

- **Documentation**: Check the tutorial and API docs
- **Issues**: Search existing issues or create new ones
- **Discussions**: Use GitHub Discussions for questions
- **Contributing**: See CONTRIBUTING.md for detailed guidelines

## External Dependencies

The tool requires several external bioinformatics tools:

### Required
- **BWA** or **minimap2**: For read alignment
- **samtools**: For BAM file processing

### Optional
- **BLAST**: For taxonomic classification
- **Python packages**: Installed automatically via pip

### Installation Commands

#### Conda (recommended)
```bash
conda install -c bioconda bwa minimap2 samtools blast
```

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install bwa minimap2 samtools ncbi-blast+
```

#### macOS (Homebrew)
```bash
brew install bwa minimap2 samtools blast
```

## Performance Tips

- Use SSD storage for large datasets
- Allocate sufficient RAM (8GB+ recommended)
- Use parallel processing with `--max-workers`
- Consider batch processing for many samples
- Pre-align reads to BAM format when possible

## Troubleshooting

### Common Issues
1. **Missing external tools**: Install BWA, minimap2, samtools
2. **Permission errors**: Check file permissions
3. **Memory issues**: Reduce batch size or parallel workers
4. **Python version**: Ensure Python 3.8+

### Getting Debug Information
```bash
# Enable debug logging
chimeric_detective --log-level DEBUG ...

# Check tool availability
which bwa minimap2 samtools

# Verify installation
python -c "import chimeric_detective; print(chimeric_detective.__version__)"
```

Happy coding! 🔬