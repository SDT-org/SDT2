import unittest
from unittest.mock import Mock, patch
import subprocess

from workflow.analyze.lzani import run
from mime_setup import register_mimetypes

register_mimetypes()


class TestLzaniAnalyze(unittest.TestCase):
    """Test lzani.py analyze module."""

    def setUp(self):
        self.mock_result = Mock()
        self.mock_result.seq_dict = {"seq1": "ATCG", "seq2": "GCTA"}
        self.mock_result.is_aa = False
        self.mock_result._replace = Mock(return_value=self.mock_result)

        self.mock_settings = Mock()
        self.mock_settings.fasta_path = "/path/to/test.fasta"
        self.mock_settings.export_alignments = False
        self.mock_settings.alignment_export_path = None

        self.mock_settings.lzani = Mock()
        self.mock_settings.lzani.exec_path = "/usr/bin/lzani"
        self.mock_settings.lzani.aw = None
        self.mock_settings.lzani.am = None
        self.mock_settings.lzani.mal = None
        self.mock_settings.lzani.msl = None
        self.mock_settings.lzani.mrd = None
        self.mock_settings.lzani.mqd = None
        self.mock_settings.lzani.reg = None
        self.mock_settings.lzani.ar = None

        self.mock_settings.doc_paths = Mock()
        self.mock_settings.doc_paths.lzani_results = "/tmp/lzani_results.tsv"
        self.mock_settings.doc_paths.lzani_results_ids = "/tmp/lzani_ids.txt"

        self.mock_set_progress = Mock()
        self.mock_canceled = Mock()

    @patch("workflow.analyze.lzani.lzani_tsv_to_distance_matrix")
    @patch("subprocess.Popen")
    def test_run_success(self, mock_popen, mock_transform):
        """Test successful lzani run."""
        # Mock subprocess
        mock_process = Mock()
        mock_process.communicate.return_value = ("success output", "")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        # Mock transformation
        mock_matrix = [[0, 10], [10, 0]]
        mock_ids = ["seq1", "seq2"]
        mock_transform.return_value = (mock_matrix, mock_ids)

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        # Check that _replace was called with correct arguments
        self.mock_result._replace.assert_called_with(
            distance_matrix=mock_matrix, ordered_ids=mock_ids
        )

    @patch("subprocess.Popen")
    def test_run_subprocess_timeout(self, mock_popen):
        """Test lzani timeout."""
        mock_process = Mock()
        # First call raises timeout, second call (after kill) returns empty strings
        mock_process.communicate.side_effect = [
            subprocess.TimeoutExpired("lzani", 300),
            ("", ""),
        ]
        mock_process.kill.return_value = None
        mock_popen.return_value = mock_process

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        self.mock_result._replace.assert_called_with(
            error="LZ-ANI timed out after 5 minutes"
        )

    @patch("subprocess.Popen")
    def test_run_subprocess_error(self, mock_popen):
        """Test lzani subprocess error."""
        mock_process = Mock()
        mock_process.communicate.return_value = ("", "Error message")
        mock_process.returncode = 1
        mock_popen.return_value = mock_process

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        self.mock_result._replace.assert_called_with(
            error="Error running LZ-ANI: Error message"
        )

    @patch("workflow.analyze.lzani.lzani_tsv_to_distance_matrix")
    @patch("subprocess.Popen")
    @patch("os.makedirs")
    def test_run_with_alignment_export(self, mock_makedirs, mock_popen, mock_transform):
        """Test lzani run with alignment export."""
        self.mock_settings.export_alignments = True
        self.mock_settings.alignment_export_path = "/tmp/alignments"

        mock_process = Mock()
        mock_process.communicate.return_value = ("success", "")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        mock_transform.return_value = ([[0, 10], [10, 0]], ["seq1", "seq2"])

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        mock_makedirs.assert_called_once_with("/tmp/alignments", exist_ok=True)
        # Check that --out-alignment was added to command
        call_args = mock_popen.call_args[0][0]
        self.assertIn("--out-alignment", call_args)

    @patch("workflow.analyze.lzani.lzani_tsv_to_distance_matrix")
    @patch("subprocess.Popen")
    def test_run_with_all_parameters(self, mock_popen, mock_transform):
        """Test lzani run with all optional parameters set."""
        self.mock_settings.lzani.aw = 5
        self.mock_settings.lzani.am = 2
        self.mock_settings.lzani.mal = 50
        self.mock_settings.lzani.msl = 50
        self.mock_settings.lzani.mrd = 0.1
        self.mock_settings.lzani.mqd = 0.01
        self.mock_settings.lzani.reg = 0
        self.mock_settings.lzani.ar = 0.95

        mock_process = Mock()
        mock_process.communicate.return_value = ("success", "")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        mock_transform.return_value = ([[0, 10], [10, 0]], ["seq1", "seq2"])

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        # Check that all parameters were added to command
        call_args = mock_popen.call_args[0][0]
        self.assertIn("--aw", call_args)
        self.assertIn("5", call_args)
        self.assertIn("--am", call_args)
        self.assertIn("2", call_args)
        self.assertIn("--mal", call_args)
        self.assertIn("50", call_args)


if __name__ == "__main__":
    unittest.main()
