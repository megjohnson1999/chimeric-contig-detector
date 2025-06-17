"""
Configuration management for Chimeric Detective.

This module handles configuration file parsing, validation, and preset management.
"""

import os
import yaml
import json
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
from dataclasses import dataclass, field, asdict
import logging

logger = logging.getLogger(__name__)


@dataclass
class DetectionConfig:
    """Configuration for chimera detection algorithms."""
    min_contig_length: int = 1000
    min_coverage: float = 5.0
    coverage_fold_change: float = 2.0
    gc_content_threshold: float = 0.1
    kmer_distance_threshold: float = 0.3
    confidence_threshold: float = 0.5
    min_split_length: int = 500


@dataclass
class ProcessingConfig:
    """Configuration for processing and performance."""
    threads: int = 1
    max_workers: int = 4
    parallel: bool = True
    keep_intermediates: bool = False
    batch_size: int = 5


@dataclass
class OutputConfig:
    """Configuration for output and reporting."""
    generate_report: bool = True
    log_level: str = "INFO"
    

@dataclass
class BehaviorConfig:
    """Configuration for chimera resolution behavior."""
    split_technical: bool = True
    split_pcr: bool = True
    preserve_biological: bool = True
    sensitivity: str = "conservative"  # conservative, balanced, sensitive, very_sensitive


@dataclass
class MultiSampleConfig:
    """Configuration for multi-sample processing."""
    multi_sample_mode: str = "separate"  # separate, merged, batch
    reads_pattern: str = "*_R{1,2}.fastq.gz"


@dataclass
class ChimericDetectiveConfig:
    """Complete configuration for Chimeric Detective."""
    detection: DetectionConfig = field(default_factory=DetectionConfig)
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    behavior: BehaviorConfig = field(default_factory=BehaviorConfig)
    multi_sample: MultiSampleConfig = field(default_factory=MultiSampleConfig)


# Predefined presets
PRESETS = {
    "small": {
        "description": "Optimized for small assemblies (<100MB, <5K contigs)",
        "config": {
            "detection": {
                "min_contig_length": 1000,
                "min_coverage": 3.0,
                "confidence_threshold": 0.5
            },
            "processing": {
                "max_workers": 8,
                "threads": 4
            },
            "behavior": {
                "sensitivity": "balanced"
            }
        }
    },
    
    "large": {
        "description": "Optimized for large assemblies (>500MB, >25K contigs)",
        "config": {
            "detection": {
                "min_contig_length": 2000,
                "min_coverage": 5.0,
                "confidence_threshold": 0.6
            },
            "processing": {
                "max_workers": 4,
                "threads": 2,
                "batch_size": 3
            },
            "behavior": {
                "sensitivity": "balanced"
            }
        }
    },
    
    "sensitive": {
        "description": "High sensitivity detection (more chimeras detected, potential false positives)",
        "config": {
            "detection": {
                "coverage_fold_change": 1.4,
                "gc_content_threshold": 0.07,
                "kmer_distance_threshold": 0.21,
                "confidence_threshold": 0.4
            },
            "behavior": {
                "sensitivity": "sensitive"
            }
        }
    },
    
    "conservative": {
        "description": "Conservative detection (fewer false positives, may miss some chimeras)",
        "config": {
            "detection": {
                "coverage_fold_change": 3.0,
                "gc_content_threshold": 0.15,
                "kmer_distance_threshold": 0.45,
                "confidence_threshold": 0.7
            },
            "behavior": {
                "sensitivity": "conservative"
            }
        }
    },
    
    "hpc": {
        "description": "Optimized for HPC/SLURM environments with limited memory",
        "config": {
            "detection": {
                "min_contig_length": 1500,
                "min_coverage": 4.0
            },
            "processing": {
                "max_workers": 2,
                "threads": 1,
                "parallel": True,
                "batch_size": 2
            },
            "output": {
                "keep_intermediates": False
            }
        }
    },
    
    "development": {
        "description": "Development and debugging settings",
        "config": {
            "detection": {
                "min_contig_length": 500,
                "confidence_threshold": 0.3
            },
            "processing": {
                "keep_intermediates": True
            },
            "output": {
                "log_level": "DEBUG"
            }
        }
    }
}


class ConfigManager:
    """Manages configuration loading, validation, and preset handling."""
    
    def __init__(self):
        self.config = ChimericDetectiveConfig()
        self._loaded_from_file = False
        self._config_file_path = None
    
    def load_config(self, config_path: Union[str, Path]) -> None:
        """Load configuration from YAML or JSON file."""
        config_path = Path(config_path)
        
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, 'r') as f:
                if config_path.suffix.lower() in ['.yaml', '.yml']:
                    data = yaml.safe_load(f)
                elif config_path.suffix.lower() == '.json':
                    data = json.load(f)
                else:
                    raise ValueError(f"Unsupported config file format: {config_path.suffix}")
            
            self._apply_config_data(data)
            self._loaded_from_file = True
            self._config_file_path = config_path
            logger.info(f"Configuration loaded from {config_path}")
            
        except Exception as e:
            raise ValueError(f"Failed to load configuration from {config_path}: {e}")
    
    def load_preset(self, preset_name: str) -> None:
        """Load a predefined preset configuration."""
        if preset_name not in PRESETS:
            available = ", ".join(PRESETS.keys())
            raise ValueError(f"Unknown preset '{preset_name}'. Available presets: {available}")
        
        preset_data = PRESETS[preset_name]["config"]
        self._apply_config_data(preset_data)
        logger.info(f"Loaded preset configuration: {preset_name}")
        logger.info(f"Preset description: {PRESETS[preset_name]['description']}")
    
    def save_config(self, output_path: Union[str, Path], format: str = "yaml") -> None:
        """Save current configuration to file."""
        output_path = Path(output_path)
        
        config_dict = self.to_dict()
        
        try:
            with open(output_path, 'w') as f:
                if format.lower() in ['yaml', 'yml']:
                    yaml.dump(config_dict, f, default_flow_style=False, indent=2)
                elif format.lower() == 'json':
                    json.dump(config_dict, f, indent=2)
                else:
                    raise ValueError(f"Unsupported format: {format}")
            
            logger.info(f"Configuration saved to {output_path}")
            
        except Exception as e:
            raise ValueError(f"Failed to save configuration to {output_path}: {e}")
    
    def apply_env_overrides(self) -> None:
        """Apply environment variable overrides."""
        env_prefix = "CHIMERIC_DETECTIVE_"
        
        # Map environment variables to config paths
        env_mappings = {
            f"{env_prefix}MIN_CONTIG_LENGTH": ("detection", "min_contig_length", int),
            f"{env_prefix}MIN_COVERAGE": ("detection", "min_coverage", float),
            f"{env_prefix}CONFIDENCE_THRESHOLD": ("detection", "confidence_threshold", float),
            f"{env_prefix}MAX_WORKERS": ("processing", "max_workers", int),
            f"{env_prefix}THREADS": ("processing", "threads", int),
            f"{env_prefix}LOG_LEVEL": ("output", "log_level", str),
            f"{env_prefix}SENSITIVITY": ("behavior", "sensitivity", str),
        }
        
        overrides_applied = []
        
        for env_var, (section, key, type_func) in env_mappings.items():
            if env_var in os.environ:
                try:
                    value = type_func(os.environ[env_var])
                    setattr(getattr(self.config, section), key, value)
                    overrides_applied.append(f"{env_var}={value}")
                except (ValueError, TypeError) as e:
                    logger.warning(f"Invalid environment variable {env_var}: {e}")
        
        if overrides_applied:
            logger.info(f"Applied environment overrides: {', '.join(overrides_applied)}")
    
    def validate_config(self) -> None:
        """Validate the current configuration."""
        errors = []
        
        # Validate detection parameters
        if not 0 < self.config.detection.confidence_threshold <= 1:
            errors.append("confidence_threshold must be between 0 and 1")
        
        if self.config.detection.min_contig_length < 100:
            errors.append("min_contig_length must be at least 100")
        
        if self.config.detection.min_coverage < 0:
            errors.append("min_coverage must be non-negative")
        
        if self.config.detection.coverage_fold_change <= 1:
            errors.append("coverage_fold_change must be greater than 1")
        
        # Validate processing parameters
        if self.config.processing.threads < 1:
            errors.append("threads must be at least 1")
        
        if self.config.processing.max_workers < 1:
            errors.append("max_workers must be at least 1")
        
        # Validate behavior parameters
        valid_sensitivities = ["conservative", "balanced", "sensitive", "very_sensitive"]
        if self.config.behavior.sensitivity not in valid_sensitivities:
            errors.append(f"sensitivity must be one of: {', '.join(valid_sensitivities)}")
        
        # Validate output parameters
        valid_log_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
        if self.config.output.log_level not in valid_log_levels:
            errors.append(f"log_level must be one of: {', '.join(valid_log_levels)}")
        
        if errors:
            raise ValueError(f"Configuration validation failed: {'; '.join(errors)}")
        
        logger.debug("Configuration validation passed")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            "detection": asdict(self.config.detection),
            "processing": asdict(self.config.processing),
            "output": asdict(self.config.output),
            "behavior": asdict(self.config.behavior),
            "multi_sample": asdict(self.config.multi_sample)
        }
    
    def to_cli_args(self) -> Dict[str, Any]:
        """Convert configuration to CLI argument format."""
        cli_args = {}
        
        # Map config to CLI parameter names
        cli_args.update({
            'min_contig_length': self.config.detection.min_contig_length,
            'min_coverage': self.config.detection.min_coverage,
            'coverage_fold_change': self.config.detection.coverage_fold_change,
            'gc_content_threshold': self.config.detection.gc_content_threshold,
            'kmer_distance_threshold': self.config.detection.kmer_distance_threshold,
            'confidence_threshold': self.config.detection.confidence_threshold,
            'min_split_length': self.config.detection.min_split_length,
            'threads': self.config.processing.threads,
            'max_workers': self.config.processing.max_workers,
            'parallel': self.config.processing.parallel,
            'keep_intermediates': self.config.processing.keep_intermediates,
            'batch_size': self.config.processing.batch_size,
            'generate_report': self.config.output.generate_report,
            'log_level': self.config.output.log_level,
            'split_technical': self.config.behavior.split_technical,
            'split_pcr': self.config.behavior.split_pcr,
            'preserve_biological': self.config.behavior.preserve_biological,
            'sensitivity': self.config.behavior.sensitivity,
            'multi_sample_mode': self.config.multi_sample.multi_sample_mode,
            'reads_pattern': self.config.multi_sample.reads_pattern,
        })
        
        return cli_args
    
    def _apply_config_data(self, data: Dict[str, Any]) -> None:
        """Apply configuration data to the current config."""
        for section_name, section_data in data.items():
            if hasattr(self.config, section_name):
                section = getattr(self.config, section_name)
                for key, value in section_data.items():
                    if hasattr(section, key):
                        setattr(section, key, value)
                    else:
                        logger.warning(f"Unknown configuration key: {section_name}.{key}")
            else:
                logger.warning(f"Unknown configuration section: {section_name}")
    
    @staticmethod
    def list_presets() -> Dict[str, str]:
        """List available presets with descriptions."""
        return {name: preset["description"] for name, preset in PRESETS.items()}
    
    @staticmethod
    def get_default_config_paths() -> List[Path]:
        """Get default configuration file search paths."""
        paths = []
        
        # Current directory
        paths.extend([
            Path.cwd() / "chimeric_detective.yaml",
            Path.cwd() / "chimeric_detective.yml",
            Path.cwd() / "chimeric_detective.json",
        ])
        
        # User config directory
        if config_dir := os.environ.get("XDG_CONFIG_HOME"):
            config_dir = Path(config_dir)
        else:
            config_dir = Path.home() / ".config"
        
        paths.extend([
            config_dir / "chimeric_detective" / "config.yaml",
            config_dir / "chimeric_detective" / "config.yml",
            config_dir / "chimeric_detective" / "config.json",
        ])
        
        return paths
    
    def auto_load_config(self) -> Optional[Path]:
        """Automatically load configuration from default locations."""
        for config_path in self.get_default_config_paths():
            if config_path.exists():
                try:
                    self.load_config(config_path)
                    return config_path
                except Exception as e:
                    logger.warning(f"Failed to load config from {config_path}: {e}")
        
        return None