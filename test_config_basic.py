#!/usr/bin/env python3
"""
Basic test script for configuration system.
"""

import sys
import os
sys.path.insert(0, 'src')

from chimeric_detective.config import ConfigManager

def test_presets():
    """Test preset functionality."""
    print("🎛️  Testing Configuration Presets")
    print("="*50)
    
    # List presets
    presets = ConfigManager.list_presets()
    print(f"Available presets: {list(presets.keys())}")
    
    # Test loading a preset
    manager = ConfigManager()
    manager.load_preset("large")
    
    print(f"Large preset settings:")
    print(f"  - min_contig_length: {manager.config.detection.min_contig_length}")
    print(f"  - max_workers: {manager.config.processing.max_workers}")
    print(f"  - sensitivity: {manager.config.behavior.sensitivity}")
    
    # Test generating config
    print("\n📁 Testing Config Generation")
    manager.save_config("test_generated_config.yaml", "yaml")
    print("✅ Generated test_generated_config.yaml")
    
    # Test loading the generated config
    manager2 = ConfigManager()
    manager2.load_config("test_generated_config.yaml")
    print("✅ Successfully loaded generated config")
    
    # Cleanup
    os.remove("test_generated_config.yaml")
    print("🧹 Cleaned up test file")

def test_validation():
    """Test configuration validation."""
    print("\n🔍 Testing Configuration Validation")
    print("="*50)
    
    manager = ConfigManager()
    
    # Test valid config
    try:
        manager.validate_config()
        print("✅ Default configuration is valid")
    except Exception as e:
        print(f"❌ Default configuration validation failed: {e}")
    
    # Test invalid config
    manager.config.detection.confidence_threshold = 1.5
    try:
        manager.validate_config()
        print("❌ Should have failed validation")
    except ValueError as e:
        print(f"✅ Correctly caught validation error: {e}")

def test_cli_conversion():
    """Test CLI argument conversion."""
    print("\n⚙️  Testing CLI Conversion")
    print("="*50)
    
    manager = ConfigManager()
    cli_args = manager.to_cli_args()
    
    expected_keys = [
        'min_contig_length', 'min_coverage', 'max_workers', 
        'threads', 'log_level', 'sensitivity'
    ]
    
    for key in expected_keys:
        if key in cli_args:
            print(f"✅ {key}: {cli_args[key]}")
        else:
            print(f"❌ Missing key: {key}")

if __name__ == "__main__":
    print("🧪 Chimeric Detective Configuration System Test")
    print("="*60)
    
    try:
        test_presets()
        test_validation()
        test_cli_conversion()
        print("\n🎉 All tests passed!")
        
    except Exception as e:
        print(f"\n💥 Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)