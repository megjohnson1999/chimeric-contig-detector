"""
Chimeric Detective: A tool for detecting and resolving chimeric contigs in viral metagenomic assemblies.
"""

__version__ = "1.0.2"
__author__ = "Chimeric Detective Team"

from .detector import ChimeraDetector
from .analyzer import ChimeraAnalyzer
from .resolver import ChimeraResolver
from .visualizer import ChimeraVisualizer
from .multi_sample import MultiSampleProcessor

# Simplified GC-only modules
from .detector_simple import SimpleChimeraDetector
from .analyzer_simple import SimpleChimeraAnalyzer
from .resolver_simple import SimpleChimeraResolver
from .visualizer_simple import SimpleChimeraVisualizer

__all__ = [
    "ChimeraDetector",
    "ChimeraAnalyzer", 
    "ChimeraResolver",
    "ChimeraVisualizer",
    "MultiSampleProcessor",
    # Simplified modules
    "SimpleChimeraDetector",
    "SimpleChimeraAnalyzer",
    "SimpleChimeraResolver",
    "SimpleChimeraVisualizer"
]