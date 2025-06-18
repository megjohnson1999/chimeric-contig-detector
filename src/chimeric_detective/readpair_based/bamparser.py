"""
BAM file parser for extracting read-pair information.

Handles various BAM formats and provides efficient streaming access to read pairs.
"""

import pysam
import logging
from typing import Iterator, Dict, List, Tuple, Optional, NamedTuple
from dataclasses import dataclass
import numpy as np
from collections import defaultdict

from .config import ReadQualityConfig


class ReadPair(NamedTuple):
    """Container for read pair information."""
    read1_start: int
    read1_end: int
    read2_start: int
    read2_end: int
    insert_size: int
    is_proper_pair: bool
    mapping_quality: int
    is_reverse: Tuple[bool, bool]
    template_length: int
    read_names: Tuple[str, str]


@dataclass
class ContigStats:
    """Statistics for a contig from BAM file."""
    name: str
    length: int
    total_pairs: int
    proper_pairs: int
    coverage_depth: float
    insert_size_mean: float
    insert_size_std: float


class BamParser:
    """Parser for extracting read-pair information from BAM files."""
    
    def __init__(self, bam_path: str, config: ReadQualityConfig):
        """
        Initialize BAM parser.
        
        Args:
            bam_path: Path to indexed BAM file
            config: Read quality configuration
        """
        self.bam_path = bam_path
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Open BAM file
        try:
            self.bamfile = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            raise ValueError(f"Failed to open BAM file: {e}")
        
        # Check if index exists
        if not self.bamfile.check_index():
            self.logger.warning("BAM index not found, creating one...")
            pysam.index(bam_path)
            self.bamfile = pysam.AlignmentFile(bam_path, "rb")
        
        # Cache contig information
        self.contigs = {ref: length for ref, length in 
                       zip(self.bamfile.references, self.bamfile.lengths)}
        
        self.logger.info(f"Loaded BAM file with {len(self.contigs)} contigs")
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
    
    def close(self):
        """Close BAM file."""
        if hasattr(self, 'bamfile'):
            self.bamfile.close()
    
    def get_contig_stats(self) -> Dict[str, ContigStats]:
        """Calculate basic statistics for all contigs."""
        stats = {}
        
        for contig, length in self.contigs.items():
            pairs = list(self.get_read_pairs(contig))
            if not pairs:
                continue
            
            proper_pairs = [p for p in pairs if p.is_proper_pair]
            insert_sizes = [p.insert_size for p in proper_pairs if p.insert_size > 0]
            
            stats[contig] = ContigStats(
                name=contig,
                length=length,
                total_pairs=len(pairs),
                proper_pairs=len(proper_pairs),
                coverage_depth=len(pairs) * 200 / length,  # Rough estimate
                insert_size_mean=np.mean(insert_sizes) if insert_sizes else 0,
                insert_size_std=np.std(insert_sizes) if insert_sizes else 0
            )
        
        return stats
    
    def get_read_pairs(self, contig: str, 
                      start: Optional[int] = None,
                      end: Optional[int] = None) -> Iterator[ReadPair]:
        """
        Extract read pairs for a contig or region.
        
        Args:
            contig: Contig name
            start: Start position (0-based, inclusive)
            end: End position (0-based, exclusive)
            
        Yields:
            ReadPair objects passing quality filters
        """
        if contig not in self.contigs:
            raise ValueError(f"Contig {contig} not found in BAM file")
        
        # Track reads by name to match pairs
        read_cache = {}
        
        # Fetch reads from region
        for read in self.bamfile.fetch(contig, start, end):
            # Apply quality filters
            if not self._passes_quality_filter(read):
                continue
            
            # Skip if not paired
            if not read.is_paired:
                continue
            
            qname = read.query_name
            
            # Check if we've seen the mate
            if qname in read_cache:
                mate = read_cache.pop(qname)
                
                # Determine which is read1 and read2
                if read.is_read1:
                    read1, read2 = read, mate
                else:
                    read1, read2 = mate, read
                
                # Create ReadPair
                pair = self._create_read_pair(read1, read2)
                if pair:
                    yield pair
            else:
                # Cache this read
                read_cache[qname] = read
        
        # Log any orphaned reads
        if read_cache:
            self.logger.debug(f"Found {len(read_cache)} orphaned reads in region")
    
    def _passes_quality_filter(self, read: pysam.AlignedSegment) -> bool:
        """Check if read passes quality filters."""
        # Check mapping quality
        if read.mapping_quality < self.config.min_mapping_quality:
            return False
        
        # Check flags
        if self.config.exclude_duplicates and read.is_duplicate:
            return False
        if self.config.exclude_secondary and read.is_secondary:
            return False
        if self.config.exclude_supplementary and read.is_supplementary:
            return False
        
        # Check if properly paired (if required)
        if self.config.require_proper_pairs and not read.is_proper_pair:
            return False
        
        # Check edit distance if specified
        if self.config.max_edit_distance is not None:
            nm_tag = read.get_tag('NM') if read.has_tag('NM') else 0
            if nm_tag > self.config.max_edit_distance:
                return False
        
        return True
    
    def _create_read_pair(self, read1: pysam.AlignedSegment, 
                         read2: pysam.AlignedSegment) -> Optional[ReadPair]:
        """Create ReadPair from two aligned segments."""
        # Both reads must be mapped
        if read1.is_unmapped or read2.is_unmapped:
            return None
        
        # Both must be on same contig
        if read1.reference_name != read2.reference_name:
            return None
        
        # Calculate insert size (outer distance)
        if read1.reference_start < read2.reference_start:
            insert_size = read2.reference_end - read1.reference_start
        else:
            insert_size = read1.reference_end - read2.reference_start
        
        return ReadPair(
            read1_start=read1.reference_start,
            read1_end=read1.reference_end,
            read2_start=read2.reference_start,
            read2_end=read2.reference_end,
            insert_size=insert_size,
            is_proper_pair=read1.is_proper_pair,
            mapping_quality=min(read1.mapping_quality, read2.mapping_quality),
            is_reverse=(read1.is_reverse, read2.is_reverse),
            template_length=abs(read1.template_length),
            read_names=(read1.query_name, read2.query_name)
        )
    
    def get_insert_size_distribution(self, contig: str, 
                                   sample_size: Optional[int] = None) -> np.ndarray:
        """
        Get insert size distribution for properly paired reads.
        
        Args:
            contig: Contig name
            sample_size: Maximum number of pairs to sample
            
        Returns:
            Array of insert sizes
        """
        insert_sizes = []
        
        for i, pair in enumerate(self.get_read_pairs(contig)):
            if pair.is_proper_pair and 0 < pair.insert_size < 10000:
                insert_sizes.append(pair.insert_size)
            
            if sample_size and len(insert_sizes) >= sample_size:
                break
        
        return np.array(insert_sizes)
    
    def stream_windows(self, contig: str, window_size: int, 
                      step_size: int) -> Iterator[Tuple[int, int, List[ReadPair]]]:
        """
        Stream read pairs in sliding windows.
        
        Args:
            contig: Contig name
            window_size: Size of sliding window
            step_size: Step between windows
            
        Yields:
            Tuples of (start, end, read_pairs)
        """
        contig_length = self.contigs[contig]
        
        for start in range(0, contig_length - window_size + 1, step_size):
            end = start + window_size
            pairs = list(self.get_read_pairs(contig, start, end))
            yield start, end, pairs


class BamValidator:
    """Validator for BAM file compatibility."""
    
    @staticmethod
    def validate_bam_file(bam_path: str) -> List[str]:
        """
        Validate BAM file for compatibility.
        
        Returns:
            List of validation issues (empty if valid)
        """
        issues = []
        
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bamfile:
                # Check if file is empty
                try:
                    next(bamfile.fetch())
                except StopIteration:
                    issues.append("BAM file appears to be empty")
                
                # Check for index
                if not bamfile.check_index():
                    issues.append("BAM file is not indexed (run samtools index)")
                
                # Check for paired reads
                paired_count = 0
                for i, read in enumerate(bamfile.fetch()):
                    if read.is_paired:
                        paired_count += 1
                    if i >= 1000:  # Sample first 1000 reads
                        break
                
                if paired_count == 0:
                    issues.append("No paired reads found in first 1000 reads")
                
                # Check for required tags
                recommended_tags = ['NM', 'AS', 'XS']
                found_tags = set()
                for read in bamfile.fetch():
                    found_tags.update(read.tags)
                    if len(found_tags) > 10:
                        break
                
                missing_tags = set(recommended_tags) - {tag[0] for tag in found_tags}
                if missing_tags:
                    issues.append(f"Recommended tags missing: {missing_tags}")
        
        except Exception as e:
            issues.append(f"Failed to open BAM file: {e}")
        
        return issues