"""
Simplified chimera analysis module for classifying GC-based chimera candidates.
"""

import logging
from typing import List, Dict, Optional
from dataclasses import dataclass
import numpy as np

from .detector_simple import ChimeraCandidate
from .utils import setup_logging


@dataclass
class ChimeraClassification:
    """Classification result for a chimera candidate."""
    candidate: ChimeraCandidate
    is_likely_chimera: bool
    classification: str  # 'technical_artifact', 'biological_recombination', 'uncertain'
    confidence: float
    reason: str


class SimpleChimeraAnalyzer:
    """Simplified analyzer for classifying chimera candidates based on GC evidence."""
    
    def __init__(self,
                 gc_difference_threshold: float = 0.15,
                 high_confidence_gc_threshold: float = 0.25,
                 log_level: str = "INFO"):
        """
        Initialize SimpleChimeraAnalyzer.
        
        Args:
            gc_difference_threshold: Minimum GC difference to consider as chimera
            high_confidence_gc_threshold: GC difference for high confidence classification
            log_level: Logging level
        """
        self.gc_difference_threshold = gc_difference_threshold
        self.high_confidence_gc_threshold = high_confidence_gc_threshold
        self.logger = setup_logging(log_level)
    
    def analyze_candidates(self, candidates: List[ChimeraCandidate]) -> List[ChimeraClassification]:
        """
        Analyze and classify chimera candidates.
        
        Args:
            candidates: List of chimera candidates to analyze
            
        Returns:
            List of classification results
        """
        self.logger.info(f"Analyzing {len(candidates)} chimera candidates")
        
        classifications = []
        for candidate in candidates:
            classification = self._classify_candidate(candidate)
            classifications.append(classification)
        
        # Log summary
        likely_chimeras = [c for c in classifications if c.is_likely_chimera]
        self.logger.info(f"Classified {len(likely_chimeras)} candidates as likely chimeras")
        
        return classifications
    
    def _classify_candidate(self, candidate: ChimeraCandidate) -> ChimeraClassification:
        """Classify a single chimera candidate based on GC evidence."""
        gc_diff = candidate.gc_difference
        
        # Determine if it's a likely chimera based on GC difference
        if gc_diff >= self.high_confidence_gc_threshold:
            is_likely_chimera = True
            classification = 'technical_artifact'
            confidence = min(1.0, gc_diff / 0.4)
            reason = f"Very high GC content difference ({gc_diff:.2f})"
        elif gc_diff >= self.gc_difference_threshold:
            is_likely_chimera = True
            classification = 'uncertain'
            confidence = 0.5 + (gc_diff - self.gc_difference_threshold) / (self.high_confidence_gc_threshold - self.gc_difference_threshold) * 0.3
            reason = f"Significant GC content difference ({gc_diff:.2f})"
        else:
            is_likely_chimera = False
            classification = 'uncertain'
            confidence = gc_diff / self.gc_difference_threshold * 0.3
            reason = f"Low GC content difference ({gc_diff:.2f})"
        
        return ChimeraClassification(
            candidate=candidate,
            is_likely_chimera=is_likely_chimera,
            classification=classification,
            confidence=confidence,
            reason=reason
        )
    
    def filter_high_confidence(self, classifications: List[ChimeraClassification],
                             min_confidence: float = 0.7) -> List[ChimeraClassification]:
        """
        Filter classifications to only include high confidence chimeras.
        
        Args:
            classifications: List of all classifications
            min_confidence: Minimum confidence threshold
            
        Returns:
            Filtered list of high confidence classifications
        """
        return [c for c in classifications if c.is_likely_chimera and c.confidence >= min_confidence]