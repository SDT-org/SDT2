import unittest
import os
from unittest.mock import Mock, patch, MagicMock
from datetime import datetime

from workflow.runner import run_parse, output_summary, run_process
from workflow.models import (
    WorkflowResult,
)
from mime_setup import register_mimetypes

register_mimetypes()


class TestRunParseFunction(unittest.TestCase):
    """Test run_parse function interface."""

    def test_run_parse_valid_fasta_returns_workflow_result(self):
        """Test that sending a valid FASTA returns a valid WorkflowResult."""
        valid_fasta_path = os.path.join(
            os.path.dirname(__file__), "..", "fastas", "valid.fasta"
        )

        result = run_parse(valid_fasta_path)

        self.assertIsInstance(result, WorkflowResult)
        self.assertIsNone(result.error)
        self.assertGreater(len(result.seq_dict), 0)


class TestRunProcess(unittest.TestCase):
    """Test run_process function."""

    def setUp(self):
        self.mock_result = Mock()
        self.mock_result.error = None

        self.mock_settings = Mock()
        self.mock_settings.analysis_method = "parasail"
        self.mock_settings.cluster_method = "single"
        self.mock_settings.fasta_path = "test.fasta"
        self.mock_settings.doc_paths = Mock()
        self.mock_settings.doc_paths.run_settings = "run_settings.json"
        self.mock_settings.doc_paths.summary = "summary.txt"
        self.mock_settings.lzani = Mock()
        self.mock_settings.lzani.score_type = "ani"

        self.mock_workflow_run = Mock()
        self.mock_workflow_run.result = self.mock_result
        self.mock_workflow_run.settings = self.mock_settings
        self.mock_workflow_run.analyze_start_time = None
        self.mock_workflow_run.analyze_start_counter = None

        self.mock_cancel_event = Mock()

    @patch("workflow.analyze.jobs")
    @patch("workflow.cluster.run")
    @patch("workflow.postprocess.run")
    @patch("workflow.runner.save_run_settings_to_json")
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("builtins.print")
    def test_run_process_success(
        self,
        mock_print,
        mock_open,
        mock_save_json,
        mock_postprocess,
        mock_cluster,
        mock_analyze_jobs,
    ):
        """Test successful run_process execution."""
        mock_analyze_job = Mock()
        mock_analyze_job.run.return_value = self.mock_result
        mock_analyze_jobs.__getitem__.return_value = mock_analyze_job

        mock_cluster.return_value = self.mock_result
        mock_postprocess.return_value = self.mock_result

        result = run_process(self.mock_workflow_run, self.mock_cancel_event)

        self.assertEqual(result, self.mock_result)
        mock_analyze_job.run.assert_called_once()
        mock_cluster.assert_called_once()
        mock_postprocess.assert_called_once()
        mock_save_json.assert_called_once()
        mock_open.assert_called_once()

    @patch("workflow.analyze.jobs")
    def test_run_process_analyze_error(self, mock_analyze_jobs):
        """Test run_process when analyze returns error."""
        error_result = Mock()
        error_result.error = "ANALYZE_ERROR"

        mock_analyze_job = Mock()
        mock_analyze_job.run.return_value = error_result
        mock_analyze_jobs.__getitem__.return_value = mock_analyze_job

        result = run_process(self.mock_workflow_run, self.mock_cancel_event)

        self.assertEqual(result, error_result)

    @patch("workflow.analyze.jobs")
    @patch("workflow.cluster.run")
    def test_run_process_cluster_error(self, mock_cluster, mock_analyze_jobs):
        """Test run_process when cluster returns error."""
        mock_analyze_job = Mock()
        mock_analyze_job.run.return_value = self.mock_result
        mock_analyze_jobs.__getitem__.return_value = mock_analyze_job

        error_result = Mock()
        error_result.error = "CLUSTER_ERROR"
        mock_cluster.return_value = error_result

        result = run_process(self.mock_workflow_run, self.mock_cancel_event)

        self.assertEqual(result, error_result)

    @patch("workflow.analyze.jobs")
    @patch("workflow.cluster.run")
    @patch("workflow.postprocess.run")
    def test_run_process_postprocess_error(
        self, mock_postprocess, mock_cluster, mock_analyze_jobs
    ):
        """Test run_process when postprocess returns error."""
        mock_analyze_job = Mock()
        mock_analyze_job.run.return_value = self.mock_result
        mock_analyze_jobs.__getitem__.return_value = mock_analyze_job

        mock_cluster.return_value = self.mock_result

        error_result = Mock()
        error_result.error = "POSTPROCESS_ERROR"
        mock_postprocess.return_value = error_result

        result = run_process(self.mock_workflow_run, self.mock_cancel_event)

        self.assertEqual(result, error_result)

    @patch("workflow.analyze.jobs")
    @patch("workflow.postprocess.run")
    @patch("workflow.runner.save_run_settings_to_json")
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    @patch("builtins.print")
    def test_run_process_no_clustering(
        self, mock_print, mock_open, mock_save_json, mock_postprocess, mock_analyze_jobs
    ):
        """Test run_process with no clustering (None or 'None')."""
        mock_analyze_job = Mock()
        mock_analyze_job.run.return_value = self.mock_result
        mock_analyze_jobs.__getitem__.return_value = mock_analyze_job

        mock_postprocess.return_value = self.mock_result

        self.mock_settings.cluster_method = None
        result = run_process(self.mock_workflow_run, self.mock_cancel_event)
        self.assertEqual(result, self.mock_result)

        self.mock_settings.cluster_method = "None"
        result = run_process(self.mock_workflow_run, self.mock_cancel_event)
        self.assertEqual(result, self.mock_result)


class TestOutputSummary(unittest.TestCase):
    """Test output_summary function."""

    def setUp(self):
        self.mock_settings = Mock()
        self.mock_settings.analysis_method = "parasail"
        self.mock_settings.parasail = Mock()
        self.mock_settings.parasail.scoring_matrix = "BLOSUM62"
        self.mock_settings.parasail.open_penalty = 10
        self.mock_settings.parasail.extend_penalty = 1

        self.mock_result = Mock()
        self.mock_result.is_aa = True

    def test_output_summary_parasail(self):
        """Test output_summary with parasail method."""
        start_time = datetime.now()
        end_time = datetime.now()

        summary = output_summary(
            "test.fasta",
            start_time,
            end_time,
            0,
            1,
            self.mock_settings,
            self.mock_result,
        )

        self.assertIn("PARASAIL", summary)
        self.assertIn("BLOSUM62", summary)

    def test_output_summary_lzani(self):
        """Test output_summary with lzani method."""
        self.mock_settings.analysis_method = "lzani"
        self.mock_settings.lzani = Mock()
        self.mock_settings.lzani.score_type = "ani"
        self.mock_settings.lzani.aw = 5
        self.mock_settings.lzani.am = 2
        self.mock_settings.lzani.mal = 50
        self.mock_settings.lzani.msl = 50
        self.mock_settings.lzani.mrd = 0.1
        self.mock_settings.lzani.mqd = 0.01
        self.mock_settings.lzani.reg = 0
        self.mock_settings.lzani.ar = 0.95

        start_time = datetime.now()
        end_time = datetime.now()

        summary = output_summary(
            "test.fasta",
            start_time,
            end_time,
            0,
            1,
            self.mock_settings,
            self.mock_result,
        )

        self.assertIn("LZ-ANI", summary)
        self.assertIn("ANI", summary)

    def test_output_summary_parasail_no_scoring_matrix_dna(self):
        """Test output_summary parasail with no scoring matrix and DNA."""
        self.mock_settings.parasail.scoring_matrix = None
        self.mock_result.is_aa = False

        start_time = datetime.now()
        end_time = datetime.now()

        summary = output_summary(
            "test.fasta",
            start_time,
            end_time,
            0,
            1,
            self.mock_settings,
            self.mock_result,
        )

        self.assertIn("1-1", summary)

    def test_output_summary_unknown_method(self):
        """Test output_summary with unknown analysis method."""
        self.mock_settings.analysis_method = "unknown_method"

        start_time = datetime.now()
        end_time = datetime.now()

        summary = output_summary(
            "test.fasta",
            start_time,
            end_time,
            0,
            1,
            self.mock_settings,
            self.mock_result,
        )

        self.assertIn("UNKNOWN", summary)
        self.assertIn("unknown_method", summary)


if __name__ == "__main__":
    unittest.main()
