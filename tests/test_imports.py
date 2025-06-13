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
    
    def test_import_detector(self):
        """Test importing detector module."""
        try:
            from chimeric_detective.detector import ChimeraDetector, ChimeraCandidate
            self.assertTrue(ChimeraDetector is not None)
            self.assertTrue(ChimeraCandidate is not None)
        except ImportError as e:
            self.fail(f"Failed to import detector: {e}")
    
    def test_import_utils(self):
        """Test importing utils module."""
        try:
            from chimeric_detective.utils import calculate_gc_content, setup_logging
            self.assertTrue(calculate_gc_content is not None)
            self.assertTrue(setup_logging is not None)
        except ImportError as e:
            self.fail(f"Failed to import utils: {e}")
    
    def test_import_cli(self):
        """Test importing CLI module."""
        try:
            from chimeric_detective.cli import main
            self.assertTrue(main is not None)
        except ImportError as e:
            self.fail(f"Failed to import CLI: {e}")


if __name__ == '__main__':
    unittest.main()