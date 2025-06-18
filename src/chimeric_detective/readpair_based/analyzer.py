"""
Read-pair analysis module for chimera detection.

Implements sliding window analysis and statistical anomaly detection.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Callable
from dataclasses import dataclass
import logging
from scipy import stats
from collections import defaultdict

from .bamparser import ReadPair, BamParser
from .config import DetectorConfig


@dataclass
class WindowMetrics:
    """Metrics calculated for a sliding window."""
    start: int
    end: int
    total_pairs: int
    proper_pairs: int
    proper_pair_ratio: float
    insert_sizes: List[int]
    insert_size_median: float
    insert_size_mad: float
    discordant_pairs: int
    discordant_ratio: float
    coverage_uniformity: float


@dataclass
class AnomalyRegion:
    """Detected anomaly region with evidence."""
    contig: str
    start: int
    end: int
    anomaly_type: str  # proper_pair_drop, insert_size_shift, discordant_spike
    score: float
    metrics: WindowMetrics
    baseline_metrics: Optional[WindowMetrics] = None


@dataclass
class BreakpointCandidate:
    """Candidate chimeric breakpoint."""
    contig: str
    position: int
    confidence: float
    supporting_anomalies: List[AnomalyRegion]
    left_metrics: WindowMetrics
    right_metrics: WindowMetrics


class ReadPairAnalyzer:
    """Analyzer for detecting anomalies in read-pair patterns."""
    
    def __init__(self, config: DetectorConfig):
        """Initialize analyzer with configuration."""
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Statistical method registry
        self.outlier_methods = {
            'mad': self._detect_outliers_mad,
            'iqr': self._detect_outliers_iqr,
            'zscore': self._detect_outliers_zscore,
            'kde': self._detect_outliers_kde
        }
        
        # Statistical test registry
        self.statistical_tests = {
            'ks': stats.ks_2samp,
            'anderson': self._anderson_darling_test,
            'chi2': self._chi_squared_test
        }
    
    def analyze_contig(self, parser: BamParser, contig: str) -> List[BreakpointCandidate]:
        """
        Analyze a contig for chimeric breakpoints.
        
        Args:
            parser: BAM parser instance
            contig: Contig name
            
        Returns:
            List of breakpoint candidates
        """
        self.logger.info(f"Analyzing contig {contig}")
        
        # Calculate baseline metrics
        baseline = self._calculate_baseline(parser, contig)
        if not baseline:
            self.logger.warning(f"Insufficient data for baseline on {contig}")
            return []
        
        # Perform sliding window analysis
        window_metrics = self._sliding_window_analysis(parser, contig)
        
        # Detect anomalies
        anomalies = self._detect_anomalies(window_metrics, baseline)
        
        # Identify breakpoints from anomalies
        breakpoints = self._identify_breakpoints(anomalies, window_metrics)
        
        return breakpoints
    
    def _calculate_baseline(self, parser: BamParser, contig: str) -> Optional[WindowMetrics]:
        """Calculate baseline metrics from properly paired reads."""
        all_pairs = []
        proper_pairs = []
        insert_sizes = []
        
        # Collect pairs for baseline
        for pair in parser.get_read_pairs(contig):
            all_pairs.append(pair)
            if pair.is_proper_pair:
                proper_pairs.append(pair)
                if 0 < pair.insert_size < self.config.insert_size.max_insert_size:
                    insert_sizes.append(pair.insert_size)
            
            # Stop if we have enough pairs
            if len(proper_pairs) >= self.config.insert_size.min_pairs_for_baseline * 2:
                break
        
        if len(proper_pairs) < self.config.insert_size.min_pairs_for_baseline:
            return None
        
        # Calculate metrics
        insert_sizes_array = np.array(insert_sizes)
        
        return WindowMetrics(
            start=0,
            end=parser.contigs[contig],
            total_pairs=len(all_pairs),
            proper_pairs=len(proper_pairs),
            proper_pair_ratio=len(proper_pairs) / len(all_pairs) if all_pairs else 0,
            insert_sizes=insert_sizes,
            insert_size_median=np.median(insert_sizes_array),
            insert_size_mad=self._calculate_mad(insert_sizes_array),
            discordant_pairs=len(all_pairs) - len(proper_pairs),
            discordant_ratio=1 - (len(proper_pairs) / len(all_pairs)) if all_pairs else 0,
            coverage_uniformity=1.0  # Baseline assumption
        )
    
    def _sliding_window_analysis(self, parser: BamParser, contig: str) -> List[WindowMetrics]:
        """Perform sliding window analysis on contig."""
        metrics_list = []
        
        window_size = self.config.sliding_window.window_size
        step_size = self.config.sliding_window.step_size
        
        for start, end, pairs in parser.stream_windows(contig, window_size, step_size):
            # Skip windows with insufficient data
            if len(pairs) < self.config.sliding_window.min_pairs_per_window:
                continue
            
            # Calculate window metrics
            metrics = self._calculate_window_metrics(start, end, pairs)
            metrics_list.append(metrics)
        
        return metrics_list
    
    def _calculate_window_metrics(self, start: int, end: int, 
                                 pairs: List[ReadPair]) -> WindowMetrics:
        """Calculate metrics for a window of read pairs."""
        proper_pairs = [p for p in pairs if p.is_proper_pair]
        insert_sizes = [p.insert_size for p in proper_pairs 
                       if 0 < p.insert_size < self.config.insert_size.max_insert_size]
        
        # Calculate coverage uniformity (simplified)
        coverage_positions = defaultdict(int)
        for pair in pairs:
            for pos in range(pair.read1_start, pair.read1_end):
                coverage_positions[pos] += 1
            for pos in range(pair.read2_start, pair.read2_end):
                coverage_positions[pos] += 1
        
        if coverage_positions:
            coverage_values = list(coverage_positions.values())
            coverage_uniformity = 1.0 - (np.std(coverage_values) / (np.mean(coverage_values) + 1))
        else:
            coverage_uniformity = 0.0
        
        return WindowMetrics(
            start=start,
            end=end,
            total_pairs=len(pairs),
            proper_pairs=len(proper_pairs),
            proper_pair_ratio=len(proper_pairs) / len(pairs) if pairs else 0,
            insert_sizes=insert_sizes,
            insert_size_median=np.median(insert_sizes) if insert_sizes else 0,
            insert_size_mad=self._calculate_mad(np.array(insert_sizes)) if insert_sizes else 0,
            discordant_pairs=len(pairs) - len(proper_pairs),
            discordant_ratio=1 - (len(proper_pairs) / len(pairs)) if pairs else 0,
            coverage_uniformity=coverage_uniformity
        )
    
    def _detect_anomalies(self, window_metrics: List[WindowMetrics], 
                         baseline: WindowMetrics) -> List[AnomalyRegion]:
        """Detect anomalous windows compared to baseline."""
        anomalies = []
        
        for metrics in window_metrics:
            # Check proper pair ratio drop
            if baseline.proper_pair_ratio > 0:
                ratio_drop = 1 - (metrics.proper_pair_ratio / baseline.proper_pair_ratio)
                if ratio_drop > self.config.detection.proper_pair_drop_threshold:
                    anomalies.append(AnomalyRegion(
                        contig="",  # Will be filled by caller
                        start=metrics.start,
                        end=metrics.end,
                        anomaly_type="proper_pair_drop",
                        score=ratio_drop,
                        metrics=metrics,
                        baseline_metrics=baseline
                    ))
            
            # Check insert size distribution shift
            if metrics.insert_sizes and baseline.insert_sizes:
                # Use configured statistical test
                test_func = self.statistical_tests[self.config.detection.statistical_test]
                statistic, p_value = test_func(metrics.insert_sizes, baseline.insert_sizes)
                
                if p_value < 0.05:  # Significant difference
                    # Calculate effect size
                    median_shift = abs(metrics.insert_size_median - baseline.insert_size_median)
                    normalized_shift = median_shift / (baseline.insert_size_mad + 1)
                    
                    if normalized_shift > self.config.detection.insert_size_shift_threshold:
                        anomalies.append(AnomalyRegion(
                            contig="",
                            start=metrics.start,
                            end=metrics.end,
                            anomaly_type="insert_size_shift",
                            score=normalized_shift,
                            metrics=metrics,
                            baseline_metrics=baseline
                        ))
            
            # Check discordant pair ratio
            if metrics.discordant_ratio > self.config.detection.discordant_pair_threshold:
                anomalies.append(AnomalyRegion(
                    contig="",
                    start=metrics.start,
                    end=metrics.end,
                    anomaly_type="discordant_spike",
                    score=metrics.discordant_ratio,
                    metrics=metrics,
                    baseline_metrics=baseline
                ))
        
        return anomalies
    
    def _identify_breakpoints(self, anomalies: List[AnomalyRegion], 
                            all_metrics: List[WindowMetrics]) -> List[BreakpointCandidate]:
        """Identify breakpoint positions from anomaly regions."""
        if not anomalies:
            return []
        
        # Group nearby anomalies
        anomaly_clusters = self._cluster_anomalies(anomalies)
        
        breakpoints = []
        for cluster in anomaly_clusters:
            # Find most likely breakpoint position in cluster
            position = self._find_breakpoint_position(cluster, all_metrics)
            
            # Get metrics flanking the breakpoint
            left_metrics, right_metrics = self._get_flanking_metrics(position, all_metrics)
            
            if left_metrics and right_metrics:
                # Calculate confidence based on multiple factors
                confidence = self._calculate_confidence(cluster, left_metrics, right_metrics)
                
                breakpoints.append(BreakpointCandidate(
                    contig="",  # Will be filled by caller
                    position=position,
                    confidence=confidence,
                    supporting_anomalies=cluster,
                    left_metrics=left_metrics,
                    right_metrics=right_metrics
                ))
        
        return breakpoints
    
    def _cluster_anomalies(self, anomalies: List[AnomalyRegion]) -> List[List[AnomalyRegion]]:
        """Group nearby anomalies into clusters."""
        if not anomalies:
            return []
        
        # Sort by start position
        sorted_anomalies = sorted(anomalies, key=lambda a: a.start)
        
        clusters = []
        current_cluster = [sorted_anomalies[0]]
        
        for anomaly in sorted_anomalies[1:]:
            # Check if anomaly is close to current cluster
            last_end = max(a.end for a in current_cluster)
            if anomaly.start - last_end <= self.config.detection.merge_distance:
                current_cluster.append(anomaly)
            else:
                clusters.append(current_cluster)
                current_cluster = [anomaly]
        
        clusters.append(current_cluster)
        return clusters
    
    def _find_breakpoint_position(self, cluster: List[AnomalyRegion], 
                                all_metrics: List[WindowMetrics]) -> int:
        """Find most likely breakpoint position within anomaly cluster."""
        # Use the position with highest combined anomaly score
        position_scores = defaultdict(float)
        
        for anomaly in cluster:
            # Add score to midpoint
            midpoint = (anomaly.start + anomaly.end) // 2
            position_scores[midpoint] += anomaly.score
        
        # Return position with highest score
        return max(position_scores.items(), key=lambda x: x[1])[0]
    
    def _get_flanking_metrics(self, position: int, 
                            all_metrics: List[WindowMetrics]) -> Tuple[Optional[WindowMetrics], 
                                                                       Optional[WindowMetrics]]:
        """Get metrics for windows flanking a position."""
        left_metrics = None
        right_metrics = None
        
        for metrics in all_metrics:
            if metrics.end <= position:
                left_metrics = metrics
            elif metrics.start >= position and right_metrics is None:
                right_metrics = metrics
        
        return left_metrics, right_metrics
    
    def _calculate_confidence(self, anomalies: List[AnomalyRegion],
                            left_metrics: WindowMetrics,
                            right_metrics: WindowMetrics) -> float:
        """Calculate confidence score for a breakpoint."""
        # Combine multiple evidence sources
        scores = []
        
        # Average anomaly scores
        if anomalies:
            avg_anomaly_score = np.mean([a.score for a in anomalies])
            scores.append(min(avg_anomaly_score, 1.0))
        
        # Difference in proper pair ratios
        if left_metrics and right_metrics:
            pair_ratio_diff = abs(left_metrics.proper_pair_ratio - right_metrics.proper_pair_ratio)
            scores.append(pair_ratio_diff)
        
        # Insert size distribution difference
        if left_metrics.insert_sizes and right_metrics.insert_sizes:
            ks_stat, p_value = stats.ks_2samp(left_metrics.insert_sizes, 
                                             right_metrics.insert_sizes)
            scores.append(1 - p_value)
        
        # Return weighted average
        return np.mean(scores) if scores else 0.0
    
    # Statistical helper methods
    
    def _calculate_mad(self, data: np.ndarray) -> float:
        """Calculate Median Absolute Deviation."""
        if len(data) == 0:
            return 0.0
        median = np.median(data)
        return np.median(np.abs(data - median))
    
    def _detect_outliers_mad(self, data: np.ndarray, threshold: float = 3.0) -> np.ndarray:
        """Detect outliers using MAD method."""
        median = np.median(data)
        mad = self._calculate_mad(data)
        if mad == 0:
            return np.zeros(len(data), dtype=bool)
        return np.abs(data - median) > threshold * mad
    
    def _detect_outliers_iqr(self, data: np.ndarray, threshold: float = 1.5) -> np.ndarray:
        """Detect outliers using IQR method."""
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        return (data < q1 - threshold * iqr) | (data > q3 + threshold * iqr)
    
    def _detect_outliers_zscore(self, data: np.ndarray, threshold: float = 3.0) -> np.ndarray:
        """Detect outliers using z-score method."""
        z_scores = np.abs(stats.zscore(data))
        return z_scores > threshold
    
    def _detect_outliers_kde(self, data: np.ndarray, threshold: float = 0.01) -> np.ndarray:
        """Detect outliers using kernel density estimation."""
        # Simplified KDE outlier detection
        kde = stats.gaussian_kde(data)
        densities = kde(data)
        return densities < np.percentile(densities, threshold * 100)
    
    def _anderson_darling_test(self, sample1: List[float], 
                              sample2: List[float]) -> Tuple[float, float]:
        """Perform Anderson-Darling test for two samples."""
        # Simplified implementation
        combined = np.concatenate([sample1, sample2])
        result = stats.anderson_ksamp([sample1, sample2])
        return result.statistic, result.significance_level
    
    def _chi_squared_test(self, sample1: List[float], 
                         sample2: List[float]) -> Tuple[float, float]:
        """Perform chi-squared test for two samples."""
        # Create bins
        combined = np.concatenate([sample1, sample2])
        bins = np.histogram_bin_edges(combined, bins='auto')
        
        # Create histograms
        hist1, _ = np.histogram(sample1, bins=bins)
        hist2, _ = np.histogram(sample2, bins=bins)
        
        # Perform test
        statistic, p_value = stats.chisquare(hist1 + 1, hist2 + 1)  # Add 1 to avoid zeros
        return statistic, p_value