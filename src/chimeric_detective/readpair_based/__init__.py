"""
Read-pair based chimera detection module.

A focused approach using only read-pair information to detect chimeric contigs
in viral metagenomic co-assemblies.
"""

from .core import ReadPairChimeraDetector
from .config import DetectorConfig
from .bamparser import BamParser
from .analyzer import ReadPairAnalyzer
from .output import OutputFormatter

__all__ = [
    "ReadPairChimeraDetector",
    "DetectorConfig",
    "BamParser",
    "ReadPairAnalyzer",
    "OutputFormatter"
]

__version__ = "2.0.0"