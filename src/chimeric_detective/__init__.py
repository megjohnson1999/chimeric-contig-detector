"""
Chimeric Detective: A tool for detecting and resolving chimeric contigs in viral metagenomic assemblies.
"""

__version__ = "1.0.2"
__author__ = "Chimeric Detective Team"

# Simplified GC-only modules
from .detector_simple import SimpleChimeraDetector
from .analyzer_simple import SimpleChimeraAnalyzer
from .resolver_simple import SimpleChimeraResolver
from .visualizer_simple import SimpleChimeraVisualizer

# Read-pair based modules (new focused approach)
try:
    from .readpair_based import (
        ReadPairChimeraDetector,
        DetectorConfig,
        BamParser,
        ReadPairAnalyzer,
        OutputFormatter
    )
    READPAIR_AVAILABLE = True
except ImportError:
    READPAIR_AVAILABLE = False

__all__ = [
    # Simplified modules
    "SimpleChimeraDetector",
    "SimpleChimeraAnalyzer",
    "SimpleChimeraResolver",
    "SimpleChimeraVisualizer"
]

# Add read-pair modules if available
if READPAIR_AVAILABLE:
    __all__.extend([
        "ReadPairChimeraDetector",
        "DetectorConfig",
        "BamParser",
        "ReadPairAnalyzer",
        "OutputFormatter"
    ])