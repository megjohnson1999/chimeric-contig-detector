# Contributing to Chimeric Detective

We welcome contributions to Chimeric Detective! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your feature or bugfix
4. Make your changes
5. Submit a pull request

## Development Setup

### Prerequisites

- Python 3.8 or later
- Git
- External bioinformatics tools (BWA, minimap2, samtools, BLAST)

### Installation for Development

```bash
# Clone your fork
git clone https://github.com/YOURUSERNAME/chimeric-contig-detector.git
cd chimeric-contig-detector

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .

# Install development dependencies
pip install pytest pytest-cov black flake8 mypy
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=chimeric_detective

# Run specific test file
pytest tests/test_detector.py
```

## Code Style

We follow Python PEP 8 style guidelines with some modifications:

- Line length: 100 characters
- Use Black for code formatting
- Use flake8 for linting
- Use type hints where appropriate

### Formatting Code

```bash
# Format code with Black
black src/chimeric_detective/

# Check with flake8
flake8 src/chimeric_detective/

# Type checking with mypy
mypy src/chimeric_detective/
```

## Contribution Guidelines

### Bug Reports

When reporting bugs, please include:

1. **Description**: Clear description of the bug
2. **Steps to reproduce**: Minimal example to reproduce the issue
3. **Expected behavior**: What should happen
4. **Actual behavior**: What actually happens
5. **Environment**: OS, Python version, package versions
6. **Log files**: Relevant log output or error messages

### Feature Requests

When requesting features, please include:

1. **Use case**: Why is this feature needed?
2. **Description**: Detailed description of the proposed feature
3. **Examples**: Example usage or interface
4. **Alternatives**: Alternative solutions you've considered

### Pull Requests

#### Before Submitting

- [ ] Code follows the project style guidelines
- [ ] Tests pass locally
- [ ] New features have appropriate tests
- [ ] Documentation is updated if needed
- [ ] Commit messages are clear and descriptive

#### Pull Request Process

1. **Create a branch**: Use a descriptive name like `feature/new-detection-method` or `bugfix/coverage-calculation`
2. **Make changes**: Implement your feature or fix
3. **Add tests**: Ensure your changes are tested
4. **Update documentation**: Update relevant docs and docstrings
5. **Commit changes**: Use clear, descriptive commit messages
6. **Push branch**: Push to your fork
7. **Create PR**: Submit pull request with detailed description

#### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Tests pass locally
- [ ] New tests added for new functionality
- [ ] Manual testing performed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No new warnings introduced
```

## Development Workflow

### Adding New Detection Methods

1. Create new method in `detector.py`
2. Add corresponding tests in `tests/test_detector.py`
3. Update documentation
4. Add example usage

### Adding New Chimera Types

1. Update `analyzer.py` classification logic
2. Add explanation templates
3. Update visualization code
4. Add tests and examples

### Improving Visualizations

1. Update `visualizer.py`
2. Test with various data types
3. Ensure accessibility compliance
4. Update documentation

## Testing Guidelines

### Test Structure

```python
# tests/test_module.py
import unittest
from chimeric_detective.module import ClassName

class TestClassName(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        pass
    
    def test_specific_functionality(self):
        """Test specific functionality with descriptive name."""
        # Arrange
        # Act
        # Assert
        pass
```

### Test Data

- Use small, synthetic test data when possible
- Include edge cases and error conditions
- Mock external dependencies when appropriate
- Test with realistic but minimal examples

### Performance Tests

For performance-critical code:
- Include benchmarks
- Test with realistic data sizes
- Monitor memory usage
- Document expected performance characteristics

## Documentation

### Docstring Style

We use NumPy-style docstrings:

```python
def function_name(param1: str, param2: int = 10) -> bool:
    """
    Brief description of function.
    
    Longer description if needed, explaining the purpose,
    algorithm, or important details.
    
    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int, optional
        Description of param2 (default is 10)
    
    Returns
    -------
    bool
        Description of return value
    
    Raises
    ------
    ValueError
        Description of when this exception is raised
    
    Examples
    --------
    >>> function_name("example", 5)
    True
    """
```

### README Updates

When adding features, update:
- Installation instructions (if needed)
- Usage examples
- Feature list
- Requirements

### Tutorial Updates

For user-facing features:
- Add examples to tutorial
- Include common use cases
- Provide troubleshooting tips
- Update command reference

## Release Process

### Version Numbering

We use semantic versioning (SemVer):
- `MAJOR.MINOR.PATCH`
- `MAJOR`: Breaking changes
- `MINOR`: New features (backward compatible)
- `PATCH`: Bug fixes (backward compatible)

### Pre-release Checklist

- [ ] All tests pass
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] Version numbers updated
- [ ] Example data tested
- [ ] Performance benchmarks run

## Community Guidelines

### Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Help newcomers
- Focus on the project goals
- Follow GitHub community guidelines

### Communication

- Use GitHub issues for bug reports and feature requests
- Use GitHub Discussions for questions and general discussion
- Be patient with responses
- Provide context and details in communications

## Getting Help

### Resources

- **Documentation**: Check the tutorial and API docs first
- **Examples**: Look at the examples directory
- **Issues**: Search existing issues for similar problems
- **Tests**: Look at test files for usage examples

### Where to Ask

1. **GitHub Issues**: Bug reports, feature requests
2. **GitHub Discussions**: Questions, ideas, general help
3. **Documentation**: Check docs first for common questions

## Recognition

Contributors will be acknowledged in:
- CONTRIBUTORS.md file
- Release notes for significant contributions
- Documentation credits

Thank you for contributing to Chimeric Detective! 🔬