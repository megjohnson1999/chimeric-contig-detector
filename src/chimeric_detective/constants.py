"""
Biological and algorithmic constants for chimeric contig detection.

This module centralizes all magic numbers to make biological parameters
clear and easy to adjust for different use cases.
"""

# =============================================================================
# SEQUENCE ANALYSIS PARAMETERS
# =============================================================================

class SequenceConstants:
    """Constants related to sequence analysis and processing."""
    
    # Minimum sequence lengths for meaningful analysis
    MIN_CONTIG_LENGTH_FOR_ANALYSIS = 100      # bp - Skip very short contigs
    MIN_SEQUENCE_FOR_GC_ANALYSIS = 20         # bp - Minimum for GC content
    MIN_SEQUENCE_FOR_KMER_ANALYSIS = 4        # bp - Minimum for k-mer analysis
    MIN_SEQUENCE_FOR_WINDOW_ANALYSIS = 20     # bp - Minimum for windowed analysis
    
    # Default contig length thresholds
    DEFAULT_MIN_CONTIG_LENGTH = 1000          # bp - Default minimum for processing
    MIN_SPLIT_CONTIG_LENGTH = 500             # bp - Minimum length after splitting
    
    # K-mer analysis parameters
    DEFAULT_KMER_SIZE = 4                     # Standard k-mer size for viral analysis
    

# =============================================================================
# WINDOW AND STEP SIZE PARAMETERS
# =============================================================================

class WindowConstants:
    """Constants for sliding window analysis."""
    
    # Adaptive window sizing
    MIN_ADAPTIVE_WINDOW = 200                 # bp - Minimum adaptive window size
    MAX_ADAPTIVE_WINDOW = 2000                # bp - Maximum adaptive window size
    WINDOW_TO_CONTIG_RATIO = 20               # Contig length / ratio = window size
    
    # Step sizes for sliding windows
    MIN_ADAPTIVE_STEP = 25                    # bp - Minimum step size
    MAX_ADAPTIVE_STEP = 50                    # bp - Maximum step size
    STEP_TO_WINDOW_RATIO = 8                  # Window size / ratio = step size
    
    # Standard window sizes for specific analyses
    DEFAULT_GC_WINDOW = 1000                  # bp - Default GC content window
    DEFAULT_COVERAGE_WINDOW = 500             # bp - Default coverage analysis window
    BREAKPOINT_ANALYSIS_WINDOW = 500          # bp - Window for breakpoint analysis
    
    # Window adjustment parameters
    WINDOW_TO_SEQUENCE_MAX_RATIO = 2          # Max window size = length / ratio
    FALLBACK_WINDOW_DIVISOR = 4               # Fallback window = length / divisor


# =============================================================================
# DETECTION THRESHOLDS
# =============================================================================

class DetectionThresholds:
    """Thresholds for chimera detection algorithms."""
    
    # Coverage analysis thresholds
    MIN_COVERAGE_DEPTH = 5.0                 # Minimum coverage for analysis
    COVERAGE_FOLD_CHANGE_THRESHOLD = 2.0     # Fold change for coverage breakpoints
    
    # GC content analysis thresholds
    GC_CONTENT_DIFFERENCE_THRESHOLD = 0.1    # Minimum GC difference for breakpoint
    
    # K-mer composition thresholds
    KMER_DISTANCE_THRESHOLD = 0.3            # Jensen-Shannon distance threshold
    
    # Confidence and classification thresholds
    CONFIDENCE_THRESHOLD = 0.5               # Minimum confidence for splitting
    HIGH_CONFIDENCE_THRESHOLD = 0.8          # High confidence classification
    MEDIUM_CONFIDENCE_THRESHOLD = 0.5        # Medium confidence classification
    
    # Signal analysis thresholds
    PEAK_DETECTION_MIN_HEIGHT = 0.1          # Minimum peak height for detection
    SIGNAL_SMOOTHING_WINDOW = 5              # Window for signal smoothing


# =============================================================================
# ALGORITHM TUNING PARAMETERS
# =============================================================================

class AlgorithmConstants:
    """Constants for algorithm behavior and performance tuning."""
    
    # Sensitivity adjustments for different presets
    SENSITIVE_PRESET_MULTIPLIER = 0.7        # Multiply thresholds for sensitive mode
    CONSERVATIVE_PRESET_MULTIPLIER = 1.3     # Multiply thresholds for conservative mode
    
    # Refinement parameters
    REFINEMENT_WINDOW_DIVISOR = 4             # Refinement window = window / divisor
    BREAKPOINT_REFINEMENT_MARGIN = 50        # bp - Safety margin for breakpoint refinement
    
    # Sub-grid detection parameters
    SUB_GRID_SCAN_STEP = 25                  # bp - Step size for sub-grid scanning
    
    # Adaptive analysis parameters
    ADAPTIVE_THRESHOLD_FACTOR = 0.8          # Factor for adaptive threshold adjustment
    

# =============================================================================
# BIOLOGICAL VALIDATION PARAMETERS
# =============================================================================

class BiologicalConstants:
    """Constants based on biological knowledge of viral sequences."""
    
    # Typical viral genome characteristics
    TYPICAL_VIRAL_GENOME_SIZE_MIN = 1000     # bp - Minimum meaningful viral sequence
    TYPICAL_VIRAL_GENOME_SIZE_MAX = 1000000  # bp - Maximum expected viral genome
    
    # Recombination and chimera characteristics
    MIN_MEANINGFUL_SEGMENT_LENGTH = 200      # bp - Minimum segment to be biologically meaningful
    TYPICAL_RECOMBINATION_BREAKPOINT_MARGIN = 100  # bp - Expected precision of recombination
    
    # Assembly artifact characteristics
    MIN_TECHNICAL_CHIMERA_SIZE = 50          # bp - Minimum size for technical artifacts
    MAX_OVERLAPPING_REGION = 20              # bp - Maximum expected overlap in assemblies


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_adaptive_window_size(contig_length: int) -> int:
    """Calculate adaptive window size based on contig length."""
    return max(
        WindowConstants.MIN_ADAPTIVE_WINDOW,
        min(WindowConstants.MAX_ADAPTIVE_WINDOW, contig_length // WindowConstants.WINDOW_TO_CONTIG_RATIO)
    )


def get_adaptive_step_size(window_size: int) -> int:
    """Calculate adaptive step size based on window size."""
    return max(
        WindowConstants.MIN_ADAPTIVE_STEP,
        min(WindowConstants.MAX_ADAPTIVE_STEP, window_size // WindowConstants.STEP_TO_WINDOW_RATIO)
    )


def adjust_window_for_sequence_length(window_size: int, sequence_length: int) -> int:
    """Adjust window size if it's too large for the sequence."""
    if window_size > sequence_length // WindowConstants.WINDOW_TO_SEQUENCE_MAX_RATIO:
        return max(
            WindowConstants.MIN_ADAPTIVE_WINDOW,
            sequence_length // WindowConstants.FALLBACK_WINDOW_DIVISOR
        )
    return window_size


def get_sensitivity_multiplier(sensitivity: str) -> float:
    """Get threshold multiplier for different sensitivity settings."""
    multipliers = {
        'sensitive': AlgorithmConstants.SENSITIVE_PRESET_MULTIPLIER,
        'conservative': AlgorithmConstants.CONSERVATIVE_PRESET_MULTIPLIER,
        'default': 1.0
    }
    return multipliers.get(sensitivity.lower(), 1.0)