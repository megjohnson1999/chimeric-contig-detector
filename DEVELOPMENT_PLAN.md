# Chimeric Detective Development Plan

**Date Created:** June 16, 2025  
**Status:** Planning Phase  
**Goal:** Transform Chimeric Detective into a robust, widely-usable bioinformatics tool

## Executive Summary

Based on codebase analysis, Chimeric Detective has strong core functionality but needs improvements in three critical areas:
1. **Testing Infrastructure** - Currently minimal test coverage (~20 test cases)
2. **Distribution & Installation** - Not published to PyPI, complex dependency management
3. **Configuration Management** - Hardcoded parameters, no config file support

## Development Phases

### Phase 1: Testing Foundation (Weeks 1-2)
**Goal**: Establish robust testing infrastructure for reliability

#### 1.1 Testing Infrastructure Setup
- [ ] Implement pytest framework with proper fixtures
- [ ] Create test data generators and mock objects  
- [ ] Set up test configuration and utilities
- [ ] Configure test discovery and execution
- [ ] Add coverage reporting with pytest-cov

#### 1.2 Core Algorithm Testing
- [ ] Unit tests for detection methods:
  - [ ] Coverage discontinuity detection (`detector.py:analyze_coverage_patterns`)
  - [ ] GC content shift analysis (`detector.py:analyze_gc_content`)
  - [ ] K-mer composition analysis (`detector.py:analyze_kmer_composition`)
  - [ ] Read orientation pattern analysis
- [ ] Edge case testing:
  - [ ] Empty assemblies and single-contig cases
  - [ ] Very short contigs (< min-length threshold)
  - [ ] Single-read coverage scenarios
  - [ ] Missing or corrupted input files
- [ ] Performance regression tests for large assemblies

#### 1.3 Integration Testing
- [ ] End-to-end workflow validation:
  - [ ] Single-sample analysis pipeline
  - [ ] Multi-sample processing modes
  - [ ] BAM vs FASTQ input workflows
- [ ] Tool integration tests:
  - [ ] minimap2/BWA fallback mechanism
  - [ ] samtools BAM processing
  - [ ] External tool error handling
- [ ] Output validation tests:
  - [ ] File format correctness
  - [ ] Results consistency across runs
  - [ ] HTML report generation

#### 1.4 CI/CD Integration
- [ ] Set up automated testing in GitHub Actions
- [ ] Add test matrix for different Python versions
- [ ] Configure test failure notifications
- [ ] Add performance benchmarking automation

**Deliverables:**
- Comprehensive test suite with >80% code coverage
- Automated CI/CD pipeline
- Test documentation and contribution guidelines

---

### Phase 2: Distribution & Accessibility (Weeks 3-4)
**Goal**: Make tool easily installable and deployable

#### 2.1 PyPI Publication Preparation
- [ ] Validate `pyproject.toml` configuration
- [ ] Ensure all dependencies are properly specified
- [ ] Create proper package metadata and descriptions
- [ ] Test installation in clean environments
- [ ] Set up automated PyPI releases via GitHub Actions

#### 2.2 PyPI Publication
- [ ] Register package name on PyPI
- [ ] Create initial release (v1.0.0)
- [ ] Publish to PyPI with proper versioning
- [ ] Validate pip installation works correctly
- [ ] Update README with correct installation instructions

#### 2.3 Containerization
- [ ] Create multi-stage Dockerfile:
  - [ ] Base image with system dependencies
  - [ ] Python environment with package installation
  - [ ] External bioinformatics tools (minimap2, BWA, samtools)
- [ ] Optimize container size and build time
- [ ] Create container for different architectures (x86_64, ARM64)
- [ ] Publish to Docker Hub and GitHub Container Registry

#### 2.4 Alternative Distribution Methods
- [ ] Create conda-forge recipe for easier bioinformatics integration
- [ ] Generate pre-built binaries using PyInstaller (optional)
- [ ] Create Singularity container for HPC environments

**Deliverables:**
- Published PyPI package with automated releases
- Multi-platform Docker containers
- Simplified installation documentation

---

### Phase 3: Configuration & Usability (Weeks 5-6)
**Goal**: Improve user experience and flexibility

#### 3.1 Configuration File Support
- [ ] Design configuration schema (YAML/JSON)
- [ ] Implement configuration file parsing
- [ ] Add support for environment variable overrides
- [ ] Create parameter validation and type checking
- [ ] Add configuration file generation command

#### 3.2 Parameter Presets
- [ ] Create preset configurations for common scenarios:
  - [ ] `--preset small` - Small assemblies (<100MB)
  - [ ] `--preset large` - Large assemblies (>500MB)  
  - [ ] `--preset sensitive` - High sensitivity detection
  - [ ] `--preset conservative` - Low false-positive rate
- [ ] Allow preset customization and user-defined presets
- [ ] Add preset documentation with use case recommendations

#### 3.3 User Experience Improvements
- [ ] Enhanced error messages with specific solutions
- [ ] Progress indicators for long-running operations
- [ ] Better logging with structured output
- [ ] Input validation with helpful error messages
- [ ] Auto-detection of optimal parameters based on assembly characteristics

#### 3.4 Documentation Enhancements
- [ ] Create user tutorial with real-world examples
- [ ] Add troubleshooting flowchart
- [ ] Create developer documentation for contributions
- [ ] Add performance tuning guide
- [ ] Create comparison guide with other tools (CheckV, etc.)

**Deliverables:**
- Flexible configuration system
- User-friendly presets and defaults
- Comprehensive documentation suite

---

## Success Metrics

### Phase 1 Success Criteria
- [ ] Test coverage >80% for core functionality
- [ ] All CI/CD tests passing
- [ ] Integration tests cover major workflows
- [ ] Performance regression tests prevent slowdowns

### Phase 2 Success Criteria  
- [ ] Package installable via `pip install chimeric-detective`
- [ ] Docker container runs successfully on major platforms
- [ ] Installation time <5 minutes including dependencies
- [ ] Clean installation in fresh environments

### Phase 3 Success Criteria
- [ ] Configuration files reduce command complexity
- [ ] Presets work out-of-box for 80% of use cases
- [ ] User documentation reduces support requests
- [ ] New user onboarding time <30 minutes

## Risk Assessment & Mitigation

### High Risk Items
1. **Dependency Hell** - External bioinformatics tools may conflict
   - *Mitigation*: Comprehensive testing across environments, Docker containers
   
2. **Performance Regression** - New features may slow core algorithms
   - *Mitigation*: Continuous performance monitoring, regression tests

3. **Breaking Changes** - Configuration changes may break existing workflows
   - *Mitigation*: Backward compatibility, migration guides, deprecation warnings

### Medium Risk Items
1. **PyPI Publication Issues** - Package conflicts or metadata problems
   - *Mitigation*: Thorough testing in isolated environments
   
2. **Container Size** - Docker images may become too large
   - *Mitigation*: Multi-stage builds, dependency optimization

## Timeline & Milestones

| Week | Milestone | Deliverable |
|------|-----------|-------------|
| 1 | Testing Infrastructure Complete | pytest framework, test fixtures |
| 2 | Core Algorithm Tests Complete | >80% coverage, CI/CD running |
| 3 | PyPI Package Published | `pip install chimeric-detective` works |
| 4 | Containerization Complete | Docker images published |
| 5 | Configuration System Complete | YAML/JSON configs, presets |
| 6 | Documentation & Polish Complete | User guides, tutorials |

## Next Steps

1. **Start with Phase 1.1** - Set up pytest testing infrastructure
2. **Create test data generators** - Build on existing synthetic data tools
3. **Establish CI/CD pipeline** - Ensure quality gates from day 1
4. **Regular progress reviews** - Weekly check-ins on milestone progress

---

**Last Updated:** June 16, 2025  
**Next Review:** Weekly milestone check-ins  
**Contact:** Development team lead