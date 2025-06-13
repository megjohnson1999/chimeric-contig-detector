"""
Utility functions for chimeric detective.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam


def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """Set up logging configuration."""
    logger = logging.getLogger("chimeric_detective")
    logger.setLevel(getattr(logging, log_level.upper()))
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    
    return logger


def validate_file_exists(filepath: str, description: str = "File") -> None:
    """Validate that a file exists and is readable."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{description} not found: {filepath}")
    if not os.access(filepath, os.R_OK):
        raise PermissionError(f"{description} is not readable: {filepath}")


def create_output_directory(output_dir: str) -> None:
    """Create output directory structure."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    Path(output_dir, "chimeric_contigs").mkdir(exist_ok=True)
    Path(output_dir, "figures").mkdir(exist_ok=True)


def calculate_gc_content(sequence: str, window_size: int = 1000) -> List[float]:
    """Calculate GC content in sliding windows."""
    gc_contents = []
    seq_len = len(sequence)
    
    for i in range(0, seq_len, window_size // 2):
        window_end = min(i + window_size, seq_len)
        window_seq = sequence[i:window_end]
        
        if len(window_seq) > 0:
            gc_count = window_seq.count('G') + window_seq.count('C')
            gc_content = gc_count / len(window_seq)
            gc_contents.append(gc_content)
        else:
            gc_contents.append(0.0)
    
    return gc_contents


def calculate_kmer_frequencies(sequence: str, k: int = 4) -> Dict[str, int]:
    """Calculate k-mer frequencies in a sequence."""
    kmers = {}
    seq_len = len(sequence)
    
    for i in range(seq_len - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:
            kmers[kmer] = kmers.get(kmer, 0) + 1
    
    return kmers


def calculate_kmer_distance(kmers1: Dict[str, int], kmers2: Dict[str, int]) -> float:
    """Calculate Jensen-Shannon distance between k-mer frequency distributions."""
    all_kmers = set(kmers1.keys()) | set(kmers2.keys())
    
    if not all_kmers:
        return 0.0
    
    total1 = sum(kmers1.values()) or 1
    total2 = sum(kmers2.values()) or 1
    
    freq1 = np.array([kmers1.get(kmer, 0) / total1 for kmer in all_kmers])
    freq2 = np.array([kmers2.get(kmer, 0) / total2 for kmer in all_kmers])
    
    # Add small pseudocount to avoid log(0)
    freq1 += 1e-10
    freq2 += 1e-10
    
    # Jensen-Shannon distance
    m = 0.5 * (freq1 + freq2)
    js_div = 0.5 * np.sum(freq1 * np.log(freq1 / m)) + 0.5 * np.sum(freq2 * np.log(freq2 / m))
    js_distance = np.sqrt(js_div)
    
    return js_distance


def run_command(command: List[str], cwd: Optional[str] = None, 
                capture_output: bool = True) -> subprocess.CompletedProcess:
    """Run a system command and return the result."""
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            capture_output=capture_output,
            text=True,
            check=True
        )
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {' '.join(command)}")
        logging.error(f"Error: {e.stderr}")
        raise


def check_external_tools() -> Dict[str, bool]:
    """Check if required external tools are available."""
    tools = {
        'bwa': False,
        'minimap2': False,
        'samtools': False
    }
    
    for tool in tools:
        try:
            subprocess.run([tool, '--version'], 
                         capture_output=True, check=True)
            tools[tool] = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass
    
    return tools


def parse_reads_pattern(reads_dir: str, pattern: str) -> List[Tuple[str, Optional[str]]]:
    """Parse reads directory with pattern to find read pairs."""
    import glob
    
    if '{1,2}' in pattern:
        # Paired-end pattern
        pattern1 = pattern.replace('{1,2}', '1')
        pattern2 = pattern.replace('{1,2}', '2')
        
        files1 = sorted(glob.glob(os.path.join(reads_dir, pattern1)))
        files2 = sorted(glob.glob(os.path.join(reads_dir, pattern2)))
        
        if len(files1) != len(files2):
            raise ValueError("Unequal number of R1 and R2 files found")
        
        return list(zip(files1, files2))
    else:
        # Single-end pattern
        files = sorted(glob.glob(os.path.join(reads_dir, pattern)))
        return [(f, None) for f in files]


def merge_overlapping_intervals(intervals: List[Tuple[int, int]], 
                              min_gap: int = 100) -> List[Tuple[int, int]]:
    """Merge overlapping or nearby intervals."""
    if not intervals:
        return []
    
    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]
    
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        
        if start <= last_end + min_gap:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    
    return merged


def calculate_n50(lengths: List[int]) -> int:
    """Calculate N50 statistic."""
    if not lengths:
        return 0
    
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length / 2
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= target_length:
            return length
    
    return sorted_lengths[-1]


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(sequence.upper()))