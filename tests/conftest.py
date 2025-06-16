"""
Pytest configuration and shared fixtures.
"""

import pytest
import tempfile
import os
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from typing import List, Dict


@pytest.fixture
def temp_dir():
    """Create temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield tmp_dir


@pytest.fixture
def sample_sequences():
    """Sample sequences for testing."""
    return {
        'normal_sequence': "ATCGATCGATCG" * 100,  # Normal 1200bp sequence
        'chimeric_sequence': "A" * 600 + "G" * 600,  # Clear GC shift chimera
        'coverage_chimera': "ATCGATCG" * 150 + "GCTAGCTA" * 150,  # Different k-mer composition
        'short_sequence': "ATCGATCG",  # Short sequence (24bp)
        'with_ns': "ATCGATCGNNNNGATCGATC",  # Sequence with Ns
    }


@pytest.fixture
def sample_fasta_file(temp_dir, sample_sequences):
    """Create a sample FASTA file for testing."""
    fasta_path = Path(temp_dir) / "test_assembly.fasta"
    
    records = []
    for seq_id, seq_str in sample_sequences.items():
        record = SeqRecord(Seq(seq_str), id=seq_id, description="")
        records.append(record)
    
    SeqIO.write(records, fasta_path, "fasta")
    return str(fasta_path)


@pytest.fixture
def coverage_data():
    """Sample coverage data for testing."""
    return {
        'uniform_coverage': np.array([10.0] * 1000),  # Uniform coverage
        'step_coverage': np.array([10.0] * 500 + [30.0] * 500),  # Step change
        'gradual_change': np.array([10.0 + i * 0.02 for i in range(1000)]),  # Gradual change
        'noisy_coverage': np.random.normal(10.0, 2.0, 1000),  # Noisy but uniform
        'low_coverage': np.array([2.0] * 1000),  # Below min threshold
    }


@pytest.fixture
def gc_profile_data():
    """Sample GC content profiles for testing."""
    return {
        'uniform_gc': [0.5] * 10,  # Uniform 50% GC
        'gc_shift': [0.3] * 5 + [0.7] * 5,  # Clear GC shift
        'gradual_gc': [0.3 + i * 0.04 for i in range(10)],  # Gradual change
        'short_profile': [0.5, 0.6],  # Too short for analysis
    }


@pytest.fixture
def kmer_profile_data():
    """Sample k-mer profiles for testing."""
    return {
        'uniform_kmers': [{'AAAA': 10, 'TTTT': 10, 'GGGG': 10, 'CCCC': 10}] * 5,
        'shifted_kmers': [
            {'AAAA': 20, 'TTTT': 5, 'GGGG': 5, 'CCCC': 5},  # A-rich
            {'AAAA': 20, 'TTTT': 5, 'GGGG': 5, 'CCCC': 5},
            {'AAAA': 5, 'TTTT': 5, 'GGGG': 20, 'CCCC': 5},   # G-rich
            {'AAAA': 5, 'TTTT': 5, 'GGGG': 20, 'CCCC': 5},
        ],
        'empty_kmers': [{}, {}, {}],  # Empty k-mer profiles
    }


@pytest.fixture
def mock_bam_file(temp_dir):
    """Create a mock BAM file for testing (empty but valid)."""
    bam_path = Path(temp_dir) / "test_reads.bam"
    # Create empty BAM file - in real tests we'd need pysam to create proper BAM
    bam_path.touch()
    return str(bam_path)


@pytest.fixture
def detector_config():
    """Standard detector configuration for testing."""
    return {
        'min_contig_length': 500,
        'min_coverage': 3.0,
        'coverage_fold_change': 2.0,
        'gc_content_threshold': 0.1,
        'kmer_distance_threshold': 0.3,
        'window_size': 100,  # Smaller window for testing
        'step_size': 50,
        'log_level': 'ERROR'  # Suppress logging during tests
    }