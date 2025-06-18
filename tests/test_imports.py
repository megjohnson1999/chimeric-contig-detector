"""
Test basic imports to ensure the package structure is correct.
"""

import unittest


class TestImports(unittest.TestCase):
    """Test that all modules can be imported without errors."""
    
    def test_import_main_package(self):
        """Test importing the main package."""
        try:
            import chimeric_detective
            self.assertTrue(hasattr(chimeric_detective, '__version__'))
        except ImportError as e:
            self.fail(f"Failed to import chimeric_detective: {e}")
    
    def test_import_simple_detector(self):
        """Test importing simple detector module."""
        try:
            from chimeric_detective.detector_simple import SimpleChimeraDetector
            self.assertTrue(SimpleChimeraDetector is not None)
        except ImportError as e:
            self.fail(f"Failed to import simple detector: {e}")
    
    def test_import_readpair_detector(self):
        """Test importing read-pair detector module."""
        try:
            from chimeric_detective.readpair_based import ReadPairChimeraDetector
            self.assertTrue(ReadPairChimeraDetector is not None)
        except ImportError as e:
            self.fail(f"Failed to import read-pair detector: {e}")
    
    def test_import_utils(self):
        """Test importing utils module."""
        try:
            from chimeric_detective.utils import calculate_gc_content, setup_logging
            self.assertTrue(calculate_gc_content is not None)
            self.assertTrue(setup_logging is not None)
        except ImportError as e:
            self.fail(f"Failed to import utils: {e}")
    
    def test_import_simple_cli(self):
        """Test importing simple CLI module."""
        try:
            from chimeric_detective.cli_simple import main
            self.assertTrue(main is not None)
        except ImportError as e:
            self.fail(f"Failed to import simple CLI: {e}")
    
    def test_import_readpair_cli(self):
        """Test importing read-pair CLI module."""
        try:
            from chimeric_detective.readpair_based.cli import main
            self.assertTrue(main is not None)
        except ImportError as e:
            self.fail(f"Failed to import read-pair CLI: {e}")


if __name__ == '__main__':
    unittest.main()