import unittest
import os

from workflow.runner import run_parse
from workflow.models import WorkflowResult
from mime_setup import register_mimetypes

register_mimetypes()


class TestRunParseFunction(unittest.TestCase):
    """Test run_parse function interface."""

    def test_run_parse_valid_fasta_returns_workflow_result(self):
        """Test that sending a valid FASTA returns a valid WorkflowResult."""
        valid_fasta_path = os.path.join(os.path.dirname(__file__), "..", "fastas", "valid.fasta")
        
        result = run_parse(valid_fasta_path)
        
        self.assertIsInstance(result, WorkflowResult)
        self.assertIsNone(result.error)
        self.assertGreater(len(result.seq_dict), 0)


if __name__ == "__main__":
    unittest.main()