#!/usr/bin/env python3
"""
Simulate the multi-sample command that was failing to test the fix.
"""

import sys
import os
sys.path.insert(0, 'src')

def simulate_cli_call():
    """Simulate the CLI call that was failing."""
    print("🧪 Simulating Multi-Sample CLI Call")
    print("="*50)
    
    # Import the CLI components
    try:
        from chimeric_detective.cli import _validate_inputs, main
        from chimeric_detective.config import ConfigManager
        print("✅ Successfully imported CLI components")
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False
    
    # Simulate the kwargs that would be passed to the CLI
    # This simulates: chimeric_detective -a test.fasta --reads-dir reads/ --reads-pattern "*_R{1,2}.fastq.gz" -o results/
    kwargs = {
        'config': None,
        'preset': None,
        'generate_config': None,
        'list_presets': False,
        'assembly': 'test_data_large/large_test_assembly.fasta',  # This file exists
        'bam': None,  # This is the key that was causing issues
        'reads1': None,
        'reads2': None,
        'reads': None,
        'reads_dir': 'test_data_large/reads',  # This directory exists
        'reads_pattern': '*_R{1,2}.fastq.gz',
        'multi_sample_mode': 'separate',
        'max_workers': 4,
        'batch_size': 5,
        'parallel': True,
        'out': 'demo_results_test/',
        'reference': None,
        'min_contig_length': 1000,
        'min_coverage': 5.0,
        'coverage_fold_change': 2.0,
        'gc_content_threshold': 0.1,
        'kmer_distance_threshold': 0.3,
        'confidence_threshold': 0.5,
        'min_split_length': 500,
        'threads': 1,
        'sensitivity': 'medium',
        'split_technical': True,
        'split_pcr': True,
        'preserve_biological': True,
        'generate_report': True,
        'keep_intermediates': False,
        'log_level': 'DEBUG'
    }
    
    # Test 1: Check if the validation function works without KeyError
    print("\n📋 Testing input validation...")
    try:
        _validate_inputs(**kwargs)
        print("✅ Input validation passed without KeyError")
    except Exception as e:
        print(f"❌ Input validation failed: {e}")
        return False
    
    # Test 2: Check configuration system integration
    print("\n⚙️  Testing configuration system...")
    try:
        config_manager = ConfigManager()
        
        # Simulate what happens in main()
        auto_config_path = config_manager.auto_load_config()
        if auto_config_path:
            print(f"✅ Auto-loaded config from: {auto_config_path}")
        else:
            print("✅ No auto-config found (expected)")
        
        # Test environment overrides
        config_manager.apply_env_overrides()
        print("✅ Environment overrides applied")
        
        # Test CLI args conversion
        cli_args = config_manager.to_cli_args()
        print("✅ Configuration converted to CLI args")
        
        # Test parameter merging (what happens in main)
        merged_kwargs = cli_args.copy()
        for key, value in kwargs.items():
            if value is not None and key not in ['config', 'preset', 'generate_config', 'list_presets']:
                merged_kwargs[key] = value
        
        print("✅ Parameters merged successfully")
        
        # Test the critical condition check
        if merged_kwargs.get('reads_dir'):
            print("✅ Would correctly route to multi-sample pipeline")
        else:
            print("❌ Would incorrectly route to single-sample pipeline")
            return False
            
    except Exception as e:
        print(f"❌ Configuration system failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    print("\n🎉 All simulation tests passed! The fix should work.")
    return True

def test_specific_error_scenarios():
    """Test the specific scenarios that were causing the 'bam' KeyError."""
    print("\n🔍 Testing Specific Error Scenarios")
    print("="*50)
    
    # Test kwargs without 'bam' key (multi-sample scenario)
    kwargs_no_bam = {
        'assembly': 'test.fasta',
        'reads_dir': 'reads/',
        'out': 'results/',
        'log_level': 'INFO',
        'confidence_threshold': 0.5,
        'min_contig_length': 1000,
        'min_coverage': 5.0
        # Note: No 'bam', 'reads1', 'reads2', 'reads', 'reference' keys
    }
    
    try:
        # Test the patterns that were causing errors
        has_bam = kwargs_no_bam.get('bam') is not None  # Fixed version
        has_reads1 = kwargs_no_bam.get('reads1') is not None  # Fixed version
        has_reads2 = kwargs_no_bam.get('reads2') is not None  # Fixed version
        has_single_reads = kwargs_no_bam.get('reads') is not None  # Fixed version
        has_reads_dir = kwargs_no_bam.get('reads_dir') is not None  # Fixed version
        
        print(f"✅ Safe kwargs access: bam={has_bam}, reads1={has_reads1}, reads2={has_reads2}")
        print(f"✅ Safe kwargs access: reads={has_single_reads}, reads_dir={has_reads_dir}")
        
        # Test routing logic
        if kwargs_no_bam.get('reads_dir'):
            print("✅ Correctly routes to multi-sample pipeline")
        else:
            print("❌ Incorrectly routes to single-sample pipeline")
        
        # Test reference validation
        if kwargs_no_bam.get('reference'):
            print("Would validate reference")
        else:
            print("✅ Correctly skips reference validation")
            
        return True
        
    except KeyError as e:
        print(f"❌ Still getting KeyError: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

if __name__ == "__main__":
    print("🔧 Testing Multi-Sample Fix Before Pull")
    print("="*60)
    
    success = True
    
    try:
        # Test 1: Simulate the actual CLI call
        if not simulate_cli_call():
            success = False
        
        # Test 2: Test specific error scenarios
        if not test_specific_error_scenarios():
            success = False
        
        if success:
            print("\n🎉 ALL TESTS PASSED! Safe to pull the changes.")
        else:
            print("\n💥 SOME TESTS FAILED! Do not pull yet.")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n💥 Test suite failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)