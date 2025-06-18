"""
Configuration management for read-pair based chimera detection.

Provides flexible configuration with sensible defaults and validation.
"""

import json
import yaml
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Any
from pathlib import Path
import logging


@dataclass
class ReadQualityConfig:
    """Configuration for read quality filtering."""
    min_mapping_quality: int = 20
    min_base_quality: int = 20
    require_proper_pairs: bool = True
    max_edit_distance: Optional[int] = None
    exclude_duplicates: bool = True
    exclude_secondary: bool = True
    exclude_supplementary: bool = True


@dataclass
class InsertSizeConfig:
    """Configuration for insert size analysis."""
    outlier_method: str = "mad"  # mad, iqr, zscore, kde
    outlier_threshold: float = 3.0
    min_pairs_for_baseline: int = 1000
    max_insert_size: int = 10000
    binning_strategy: str = "auto"  # auto, fixed, adaptive


@dataclass
class SlidingWindowConfig:
    """Configuration for sliding window analysis."""
    window_size: int = 1000
    step_size: int = 500
    min_pairs_per_window: int = 50
    window_type: str = "fixed"  # fixed, adaptive, coverage-based
    edge_behavior: str = "truncate"  # truncate, pad, skip


@dataclass
class DetectionConfig:
    """Configuration for anomaly detection thresholds."""
    proper_pair_drop_threshold: float = 0.5
    insert_size_shift_threshold: float = 2.0
    discordant_pair_threshold: float = 0.3
    min_anomaly_length: int = 100
    merge_distance: int = 500
    statistical_test: str = "ks"  # ks, anderson, chi2


@dataclass
class OutputConfig:
    """Configuration for output formatting."""
    formats: List[str] = field(default_factory=lambda: ["json", "tsv"])
    include_debug_info: bool = False
    visualization: bool = True
    decimal_places: int = 4
    compress_output: bool = False
    output_prefix: str = "chimera_detection"


@dataclass
class LoggingConfig:
    """Configuration for logging behavior."""
    level: str = "INFO"
    file: Optional[str] = None
    format: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    debug_bam_parsing: bool = False
    profile_performance: bool = False


@dataclass
class DetectorConfig:
    """Main configuration container for the detector."""
    read_quality: ReadQualityConfig = field(default_factory=ReadQualityConfig)
    insert_size: InsertSizeConfig = field(default_factory=InsertSizeConfig)
    sliding_window: SlidingWindowConfig = field(default_factory=SlidingWindowConfig)
    detection: DetectionConfig = field(default_factory=DetectionConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    
    # Runtime options
    threads: int = 1
    memory_limit_gb: Optional[float] = None
    temp_dir: Optional[str] = None
    
    @classmethod
    def from_file(cls, config_path: str) -> "DetectorConfig":
        """Load configuration from YAML or JSON file."""
        path = Path(config_path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        with open(path) as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)
            elif path.suffix == '.json':
                data = json.load(f)
            else:
                raise ValueError(f"Unsupported config format: {path.suffix}")
        
        return cls.from_dict(data)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DetectorConfig":
        """Create configuration from dictionary."""
        config = cls()
        
        # Update nested configs
        if 'read_quality' in data:
            config.read_quality = ReadQualityConfig(**data['read_quality'])
        if 'insert_size' in data:
            config.insert_size = InsertSizeConfig(**data['insert_size'])
        if 'sliding_window' in data:
            config.sliding_window = SlidingWindowConfig(**data['sliding_window'])
        if 'detection' in data:
            config.detection = DetectionConfig(**data['detection'])
        if 'output' in data:
            config.output = OutputConfig(**data['output'])
        if 'logging' in data:
            config.logging = LoggingConfig(**data['logging'])
        
        # Update runtime options
        for key in ['threads', 'memory_limit_gb', 'temp_dir']:
            if key in data:
                setattr(config, key, data[key])
        
        return config
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return asdict(self)
    
    def save(self, path: str, format: str = "yaml") -> None:
        """Save configuration to file."""
        data = self.to_dict()
        
        with open(path, 'w') as f:
            if format == "yaml":
                yaml.dump(data, f, default_flow_style=False, sort_keys=False)
            elif format == "json":
                json.dump(data, f, indent=2)
            else:
                raise ValueError(f"Unsupported format: {format}")
    
    def validate(self) -> List[str]:
        """Validate configuration and return list of issues."""
        issues = []
        
        # Validate read quality
        if self.read_quality.min_mapping_quality < 0:
            issues.append("min_mapping_quality must be >= 0")
        
        # Validate insert size
        if self.insert_size.outlier_method not in ["mad", "iqr", "zscore", "kde"]:
            issues.append(f"Unknown outlier method: {self.insert_size.outlier_method}")
        if self.insert_size.outlier_threshold <= 0:
            issues.append("outlier_threshold must be > 0")
        
        # Validate sliding window
        if self.sliding_window.window_size <= 0:
            issues.append("window_size must be > 0")
        if self.sliding_window.step_size <= 0:
            issues.append("step_size must be > 0")
        if self.sliding_window.step_size > self.sliding_window.window_size:
            issues.append("step_size should not exceed window_size")
        
        # Validate detection thresholds
        for attr in ['proper_pair_drop_threshold', 'insert_size_shift_threshold', 
                     'discordant_pair_threshold']:
            value = getattr(self.detection, attr)
            if not 0 <= value <= 1:
                issues.append(f"{attr} must be between 0 and 1")
        
        # Validate output formats
        valid_formats = ["json", "tsv", "csv", "bed", "vcf"]
        for fmt in self.output.formats:
            if fmt not in valid_formats:
                issues.append(f"Unknown output format: {fmt}")
        
        # Validate logging level
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if self.logging.level.upper() not in valid_levels:
            issues.append(f"Invalid logging level: {self.logging.level}")
        
        return issues


def create_argument_parser():
    """Create argument parser that can override config file settings."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Read-pair based chimera detection for viral metagenomes"
    )
    
    # Required arguments
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("assembly", help="Assembly FASTA file")
    
    # Config file
    parser.add_argument("-c", "--config", help="Configuration file (YAML/JSON)")
    
    # Override specific parameters
    parser.add_argument("--min-mapping-quality", type=int,
                       help="Override minimum mapping quality")
    parser.add_argument("--window-size", type=int,
                       help="Override sliding window size")
    parser.add_argument("--step-size", type=int,
                       help="Override sliding window step size")
    parser.add_argument("--outlier-method", choices=["mad", "iqr", "zscore", "kde"],
                       help="Override outlier detection method")
    parser.add_argument("--threads", type=int,
                       help="Number of threads to use")
    
    # Output options
    parser.add_argument("-o", "--output-dir", default="chimera_output",
                       help="Output directory")
    parser.add_argument("--output-formats", nargs="+",
                       choices=["json", "tsv", "csv", "bed", "vcf"],
                       help="Output formats")
    
    # Logging
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose logging")
    parser.add_argument("--debug", action="store_true",
                       help="Enable debug mode")
    
    return parser