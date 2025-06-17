"""
Tests for configuration management system.
"""

import pytest
import tempfile
import os
import json
import yaml
from pathlib import Path
from unittest.mock import patch

from chimeric_detective.config import (
    ConfigManager, ChimericDetectiveConfig, DetectionConfig, 
    ProcessingConfig, OutputConfig, BehaviorConfig, MultiSampleConfig,
    PRESETS
)


class TestConfigDataClasses:
    """Test configuration data classes."""
    
    def test_detection_config_defaults(self):
        """Test DetectionConfig default values."""
        config = DetectionConfig()
        assert config.min_contig_length == 1000
        assert config.min_coverage == 5.0
        assert config.coverage_fold_change == 2.0
        assert config.gc_content_threshold == 0.1
        assert config.kmer_distance_threshold == 0.3
        assert config.confidence_threshold == 0.5
        assert config.min_split_length == 500
    
    def test_processing_config_defaults(self):
        """Test ProcessingConfig default values."""
        config = ProcessingConfig()
        assert config.threads == 1
        assert config.max_workers == 4
        assert config.parallel is True
        assert config.keep_intermediates is False
        assert config.batch_size == 5
    
    def test_output_config_defaults(self):
        """Test OutputConfig default values."""
        config = OutputConfig()
        assert config.generate_report is True
        assert config.log_level == "INFO"
    
    def test_behavior_config_defaults(self):
        """Test BehaviorConfig default values."""
        config = BehaviorConfig()
        assert config.split_technical is True
        assert config.split_pcr is True
        assert config.preserve_biological is True
        assert config.sensitivity == "conservative"
    
    def test_multi_sample_config_defaults(self):
        """Test MultiSampleConfig default values."""
        config = MultiSampleConfig()
        assert config.multi_sample_mode == "separate"
        assert config.reads_pattern == "*_R{1,2}.fastq.gz"
    
    def test_main_config_composition(self):
        """Test that main config properly composes all sub-configs."""
        config = ChimericDetectiveConfig()
        assert isinstance(config.detection, DetectionConfig)
        assert isinstance(config.processing, ProcessingConfig)
        assert isinstance(config.output, OutputConfig)
        assert isinstance(config.behavior, BehaviorConfig)
        assert isinstance(config.multi_sample, MultiSampleConfig)


class TestConfigManager:
    """Test ConfigManager functionality."""
    
    def test_config_manager_initialization(self):
        """Test ConfigManager initialization."""
        manager = ConfigManager()
        assert isinstance(manager.config, ChimericDetectiveConfig)
        assert manager._loaded_from_file is False
        assert manager._config_file_path is None
    
    def test_to_dict_conversion(self):
        """Test configuration to dictionary conversion."""
        manager = ConfigManager()
        config_dict = manager.to_dict()
        
        assert "detection" in config_dict
        assert "processing" in config_dict
        assert "output" in config_dict
        assert "behavior" in config_dict
        assert "multi_sample" in config_dict
        
        # Check that values are present
        assert config_dict["detection"]["min_contig_length"] == 1000
        assert config_dict["processing"]["max_workers"] == 4
        assert config_dict["output"]["log_level"] == "INFO"
    
    def test_to_cli_args_conversion(self):
        """Test configuration to CLI args conversion."""
        manager = ConfigManager()
        cli_args = manager.to_cli_args()
        
        # Check that CLI parameter names are correctly mapped
        assert "min_contig_length" in cli_args
        assert "max_workers" in cli_args
        assert "log_level" in cli_args
        assert "sensitivity" in cli_args
        
        # Check values
        assert cli_args["min_contig_length"] == 1000
        assert cli_args["max_workers"] == 4
        assert cli_args["log_level"] == "INFO"
    
    def test_yaml_config_loading(self):
        """Test loading configuration from YAML file."""
        config_data = {
            "detection": {
                "min_contig_length": 2000,
                "min_coverage": 10.0
            },
            "processing": {
                "max_workers": 8,
                "threads": 4
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_file = f.name
        
        try:
            manager = ConfigManager()
            manager.load_config(config_file)
            
            assert manager.config.detection.min_contig_length == 2000
            assert manager.config.detection.min_coverage == 10.0
            assert manager.config.processing.max_workers == 8
            assert manager.config.processing.threads == 4
            assert manager._loaded_from_file is True
            
        finally:
            os.unlink(config_file)
    
    def test_json_config_loading(self):
        """Test loading configuration from JSON file."""
        config_data = {
            "behavior": {
                "sensitivity": "sensitive",
                "split_technical": False
            },
            "output": {
                "log_level": "DEBUG"
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config_data, f)
            config_file = f.name
        
        try:
            manager = ConfigManager()
            manager.load_config(config_file)
            
            assert manager.config.behavior.sensitivity == "sensitive"
            assert manager.config.behavior.split_technical is False
            assert manager.config.output.log_level == "DEBUG"
            
        finally:
            os.unlink(config_file)
    
    def test_config_file_not_found(self):
        """Test error handling for missing config file."""
        manager = ConfigManager()
        
        with pytest.raises(FileNotFoundError):
            manager.load_config("nonexistent_config.yaml")
    
    def test_invalid_config_format(self):
        """Test error handling for invalid config file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("invalid config content")
            config_file = f.name
        
        try:
            manager = ConfigManager()
            with pytest.raises(ValueError, match="Unsupported config file format"):
                manager.load_config(config_file)
        finally:
            os.unlink(config_file)
    
    def test_config_saving_yaml(self):
        """Test saving configuration to YAML file."""
        manager = ConfigManager()
        manager.config.detection.min_contig_length = 3000
        manager.config.processing.max_workers = 16
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            output_file = f.name
        
        try:
            manager.save_config(output_file, 'yaml')
            
            # Load back and verify
            with open(output_file, 'r') as f:
                loaded_data = yaml.safe_load(f)
            
            assert loaded_data["detection"]["min_contig_length"] == 3000
            assert loaded_data["processing"]["max_workers"] == 16
            
        finally:
            os.unlink(output_file)
    
    def test_config_saving_json(self):
        """Test saving configuration to JSON file."""
        manager = ConfigManager()
        manager.config.behavior.sensitivity = "conservative"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            output_file = f.name
        
        try:
            manager.save_config(output_file, 'json')
            
            # Load back and verify
            with open(output_file, 'r') as f:
                loaded_data = json.load(f)
            
            assert loaded_data["behavior"]["sensitivity"] == "conservative"
            
        finally:
            os.unlink(output_file)


class TestConfigPresets:
    """Test configuration preset functionality."""
    
    def test_list_presets(self):
        """Test listing available presets."""
        presets = ConfigManager.list_presets()
        
        assert isinstance(presets, dict)
        assert "small" in presets
        assert "large" in presets
        assert "sensitive" in presets
        assert "conservative" in presets
        assert "hpc" in presets
        assert "development" in presets
        
        # Check that descriptions are strings
        for name, description in presets.items():
            assert isinstance(description, str)
            assert len(description) > 0
    
    def test_load_small_preset(self):
        """Test loading the 'small' preset."""
        manager = ConfigManager()
        manager.load_preset("small")
        
        expected = PRESETS["small"]["config"]
        assert manager.config.detection.min_contig_length == expected["detection"]["min_contig_length"]
        assert manager.config.processing.max_workers == expected["processing"]["max_workers"]
        assert manager.config.behavior.sensitivity == expected["behavior"]["sensitivity"]
    
    def test_load_large_preset(self):
        """Test loading the 'large' preset."""
        manager = ConfigManager()
        manager.load_preset("large")
        
        expected = PRESETS["large"]["config"]
        assert manager.config.detection.min_contig_length == expected["detection"]["min_contig_length"]
        assert manager.config.processing.max_workers == expected["processing"]["max_workers"]
        assert manager.config.detection.confidence_threshold == expected["detection"]["confidence_threshold"]
    
    def test_load_sensitive_preset(self):
        """Test loading the 'sensitive' preset."""
        manager = ConfigManager()
        manager.load_preset("sensitive")
        
        expected = PRESETS["sensitive"]["config"]
        assert manager.config.detection.coverage_fold_change == expected["detection"]["coverage_fold_change"]
        assert manager.config.behavior.sensitivity == expected["behavior"]["sensitivity"]
    
    def test_load_conservative_preset(self):
        """Test loading the 'conservative' preset."""
        manager = ConfigManager()
        manager.load_preset("conservative")
        
        expected = PRESETS["conservative"]["config"]
        assert manager.config.detection.coverage_fold_change == expected["detection"]["coverage_fold_change"]
        assert manager.config.behavior.sensitivity == expected["behavior"]["sensitivity"]
    
    def test_load_invalid_preset(self):
        """Test error handling for invalid preset name."""
        manager = ConfigManager()
        
        with pytest.raises(ValueError, match="Unknown preset 'invalid_preset'"):
            manager.load_preset("invalid_preset")


class TestConfigValidation:
    """Test configuration validation."""
    
    def test_valid_config_passes_validation(self):
        """Test that a valid configuration passes validation."""
        manager = ConfigManager()
        # Should not raise an exception
        manager.validate_config()
    
    def test_invalid_confidence_threshold(self):
        """Test validation of confidence threshold."""
        manager = ConfigManager()
        
        # Test too low
        manager.config.detection.confidence_threshold = 0.0
        with pytest.raises(ValueError, match="confidence_threshold must be between 0 and 1"):
            manager.validate_config()
        
        # Test too high
        manager.config.detection.confidence_threshold = 1.5
        with pytest.raises(ValueError, match="confidence_threshold must be between 0 and 1"):
            manager.validate_config()
    
    def test_invalid_min_contig_length(self):
        """Test validation of minimum contig length."""
        manager = ConfigManager()
        manager.config.detection.min_contig_length = 50
        
        with pytest.raises(ValueError, match="min_contig_length must be at least 100"):
            manager.validate_config()
    
    def test_invalid_coverage_fold_change(self):
        """Test validation of coverage fold change."""
        manager = ConfigManager()
        manager.config.detection.coverage_fold_change = 0.5
        
        with pytest.raises(ValueError, match="coverage_fold_change must be greater than 1"):
            manager.validate_config()
    
    def test_invalid_sensitivity(self):
        """Test validation of sensitivity level."""
        manager = ConfigManager()
        manager.config.behavior.sensitivity = "invalid"
        
        with pytest.raises(ValueError, match="sensitivity must be one of: conservative, balanced, sensitive, very_sensitive"):
            manager.validate_config()
    
    def test_invalid_log_level(self):
        """Test validation of log level."""
        manager = ConfigManager()
        manager.config.output.log_level = "INVALID"
        
        with pytest.raises(ValueError, match="log_level must be one of: DEBUG, INFO, WARNING, ERROR"):
            manager.validate_config()


class TestEnvironmentOverrides:
    """Test environment variable override functionality."""
    
    def test_env_override_min_contig_length(self):
        """Test environment override for min_contig_length."""
        with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_MIN_CONTIG_LENGTH": "2500"}):
            manager = ConfigManager()
            manager.apply_env_overrides()
            
            assert manager.config.detection.min_contig_length == 2500
    
    def test_env_override_max_workers(self):
        """Test environment override for max_workers."""
        with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_MAX_WORKERS": "12"}):
            manager = ConfigManager()
            manager.apply_env_overrides()
            
            assert manager.config.processing.max_workers == 12
    
    def test_env_override_log_level(self):
        """Test environment override for log_level."""
        with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_LOG_LEVEL": "DEBUG"}):
            manager = ConfigManager()
            manager.apply_env_overrides()
            
            assert manager.config.output.log_level == "DEBUG"
    
    def test_env_override_sensitivity(self):
        """Test environment override for sensitivity."""
        with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_SENSITIVITY": "sensitive"}):
            manager = ConfigManager()
            manager.apply_env_overrides()
            
            assert manager.config.behavior.sensitivity == "sensitive"
    
    def test_invalid_env_override(self):
        """Test handling of invalid environment override values."""
        with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_MAX_WORKERS": "invalid"}):
            manager = ConfigManager()
            # Should not crash, just log warning
            manager.apply_env_overrides()
            
            # Value should remain default
            assert manager.config.processing.max_workers == 4


class TestConfigIntegration:
    """Test configuration system integration."""
    
    def test_config_precedence_order(self):
        """Test that configuration precedence order works correctly."""
        # Create a config file
        config_data = {
            "detection": {"min_contig_length": 1500},
            "processing": {"max_workers": 6}
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_file = f.name
        
        try:
            manager = ConfigManager()
            
            # 1. Load from file
            manager.load_config(config_file)
            assert manager.config.detection.min_contig_length == 1500
            assert manager.config.processing.max_workers == 6
            
            # 2. Apply preset (should override file)
            manager.load_preset("large")
            assert manager.config.detection.min_contig_length == 2000  # From preset
            assert manager.config.processing.max_workers == 4  # From preset
            
            # 3. Apply environment override (should override preset)
            with patch.dict(os.environ, {"CHIMERIC_DETECTIVE_MAX_WORKERS": "8"}):
                manager.apply_env_overrides()
                assert manager.config.processing.max_workers == 8
            
        finally:
            os.unlink(config_file)
    
    def test_partial_config_application(self):
        """Test that partial configurations don't override unspecified values."""
        manager = ConfigManager()
        original_min_coverage = manager.config.detection.min_coverage
        
        # Apply config that only specifies some values
        config_data = {
            "detection": {"min_contig_length": 2000}
            # Note: min_coverage not specified
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_file = f.name
        
        try:
            manager.load_config(config_file)
            
            # Should have updated the specified value
            assert manager.config.detection.min_contig_length == 2000
            
            # Should have kept the original value for unspecified
            assert manager.config.detection.min_coverage == original_min_coverage
            
        finally:
            os.unlink(config_file)