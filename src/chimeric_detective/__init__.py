"""
Chimeric Detective: A tool for detecting and resolving chimeric contigs in viral metagenomic assemblies.
"""

__version__ = "1.0.0"
__author__ = "Chimeric Detective Team"

from .detector import ChimeraDetector
from .analyzer import ChimeraAnalyzer
from .resolver import ChimeraResolver
from .visualizer import ChimeraVisualizer
from .multi_sample import MultiSampleProcessor

__all__ = [
    "ChimeraDetector",
    "ChimeraAnalyzer", 
    "ChimeraResolver",
    "ChimeraVisualizer",
    "MultiSampleProcessor"
]