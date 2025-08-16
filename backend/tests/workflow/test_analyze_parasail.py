import unittest
from unittest.mock import Mock, patch

from workflow.analyze.parasail import (
    run,
    supports_striped_32,
    process_pair,
    get_stats_score,
    get_similarity,
    get_traceback_score,
)
from mime_setup import register_mimetypes

register_mimetypes()


class TestParasailAnalyze(unittest.TestCase):
    """Test parasail.py analyze module."""

    def setUp(self):
        self.mock_result = Mock()
        self.mock_result.seq_dict = {"seq1": "ATCG", "seq2": "GCTA"}
        self.mock_result.is_aa = False
        self.mock_result._replace = Mock(return_value=self.mock_result)

        self.mock_settings = Mock()
        self.mock_settings.export_alignments = False
        self.mock_settings.alignment_export_path = None

        self.mock_settings.parasail = Mock()
        self.mock_settings.parasail.scoring_matrix = None
        self.mock_settings.parasail.open_penalty = None
        self.mock_settings.parasail.extend_penalty = None
        self.mock_settings.parasail.process_count = 1

        self.mock_set_progress = Mock()
        self.mock_canceled = Mock()
        self.mock_canceled.value = False

    @patch("workflow.analyze.parasail.Pool")
    def test_run_success(self, mock_pool_class):
        """Test successful parasail run."""
        mock_pool = Mock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Mock pool results
        mock_pool.imap.return_value = [
            (["seq1", "seq1"], 0.0),
            (["seq1", "seq2"], 25.0),
            (["seq2", "seq2"], 0.0),
        ]

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        # Check that _replace was called
        self.mock_result._replace.assert_called()

    @patch("workflow.analyze.parasail.Pool")
    def test_run_canceled(self, mock_pool_class):
        """Test parasail run canceled."""
        mock_pool = Mock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Set canceled after first iteration
        def side_effect():
            self.mock_canceled.value = True
            return [(["seq1", "seq1"], 0.0)]

        mock_pool.imap.return_value = side_effect()

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        self.mock_result._replace.assert_called_with(
            distance_matrix=None, error="PROCESS_CANCELED"
        )

    @patch("workflow.analyze.parasail.Pool")
    @patch("os.makedirs")
    def test_run_with_alignment_export(self, mock_makedirs, mock_pool_class):
        """Test parasail run with alignment export."""
        self.mock_settings.export_alignments = True
        self.mock_settings.alignment_export_path = "/tmp/alignments"

        mock_pool = Mock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool
        mock_pool.imap.return_value = [(["seq1", "seq2"], 25.0)]

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        mock_makedirs.assert_called_once_with("/tmp/alignments", exist_ok=True)

    @patch("workflow.analyze.parasail.Pool")
    def test_run_with_negative_score_min_score_update(self, mock_pool_class):
        """Test parasail run with negative score that updates min_score."""
        mock_pool = Mock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Return results with a negative score to trigger min_score update at line 76
        mock_pool.imap.return_value = [
            (["seq1", "seq1"], 0.0),
            (["seq1", "seq2"], -5.0),  # Negative score should update min_score
            (["seq2", "seq2"], 0.0),
        ]

        result = run(
            self.mock_result,
            self.mock_settings,
            self.mock_set_progress,
            self.mock_canceled,
        )

        self.mock_result._replace.assert_called()

    @patch("parasail.nw_trace")
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_process_pair_with_alignment_export_success(self, mock_open, mock_nw_trace):
        """Test process_pair with successful alignment export."""
        # Mock parasail traceback result
        mock_traceback_result = Mock()
        mock_traceback_result.traceback.query = "ATCG"
        mock_traceback_result.traceback.ref = "ATCC"
        mock_nw_trace.return_value = mock_traceback_result

        id_sequence_pair = [["seq1", "seq2"], ["ATCG", "ATCC"]]

        with patch("workflow.analyze.parasail.get_traceback_score", return_value=25.0):
            result = process_pair(
                id_sequence_pair,
                matrix_id="blosum62",
                open_penalty=10,
                extend_penalty=1,
                use_traceback=True,
                export_alignments=True,
                alignment_export_path="/tmp/alignments",
            )

        self.assertEqual(result, (["seq1", "seq2"], 25.0))
        mock_open.assert_called_once()

    @patch("parasail.nw_trace")
    @patch("builtins.open", side_effect=Exception("File write error"))
    @patch("builtins.print")
    def test_process_pair_with_alignment_export_failure(
        self, mock_print, mock_open, mock_nw_trace
    ):
        """Test process_pair with alignment export failure."""
        # Mock parasail traceback result
        mock_traceback_result = Mock()
        mock_traceback_result.traceback.query = "ATCG"
        mock_traceback_result.traceback.ref = "ATCC"
        mock_nw_trace.return_value = mock_traceback_result

        id_sequence_pair = [["seq1", "seq2"], ["ATCG", "ATCC"]]

        with patch("workflow.analyze.parasail.get_traceback_score", return_value=25.0):
            result = process_pair(
                id_sequence_pair,
                matrix_id="blosum62",
                open_penalty=10,
                extend_penalty=1,
                use_traceback=True,
                export_alignments=True,
                alignment_export_path="/tmp/alignments",
            )

        self.assertEqual(result, (["seq1", "seq2"], 25.0))
        mock_print.assert_called_once()

    def test_supports_striped_32(self):
        """Test supports_striped_32 function."""
        with patch("parasail.can_use_sse41", return_value=True):
            self.assertTrue(supports_striped_32())

        with patch("parasail.can_use_sse41", return_value=False), patch(
            "parasail.can_use_avx2", return_value=True
        ):
            self.assertTrue(supports_striped_32())

        with patch("parasail.can_use_sse41", return_value=False), patch(
            "parasail.can_use_avx2", return_value=False
        ), patch("parasail.can_use_neon", return_value=True):
            self.assertTrue(supports_striped_32())

        with patch("parasail.can_use_sse41", return_value=False), patch(
            "parasail.can_use_avx2", return_value=False
        ), patch("parasail.can_use_neon", return_value=False):
            self.assertFalse(supports_striped_32())

    def test_process_pair_identical_sequences(self):
        """Test process_pair with identical sequences."""
        id_sequence_pair = [["seq1", "seq1"], ["ATCG", "ATCG"]]

        result = process_pair(
            id_sequence_pair, matrix_id="simple_2_-1", open_penalty=8, extend_penalty=1
        )

        self.assertEqual(result, (["seq1", "seq1"], 0.0))

    @patch("workflow.analyze.parasail.get_stats_score")
    def test_process_pair_different_sequences_stats(self, mock_get_stats):
        """Test process_pair with different sequences using stats."""
        mock_get_stats.return_value = 25.0

        id_sequence_pair = [["seq1", "seq2"], ["ATCG", "GCTA"]]

        result = process_pair(
            id_sequence_pair,
            matrix_id="simple_2_-1",
            open_penalty=8,
            extend_penalty=1,
            use_traceback=False,
        )

        self.assertEqual(result, (["seq1", "seq2"], 25.0))
        mock_get_stats.assert_called_once()

    @patch("workflow.analyze.parasail.get_traceback_score")
    def test_process_pair_different_sequences_traceback(self, mock_get_traceback):
        """Test process_pair with different sequences using traceback."""
        mock_get_traceback.return_value = 25.0

        id_sequence_pair = [["seq1", "seq2"], ["ATCG", "GCTA"]]

        result = process_pair(
            id_sequence_pair,
            matrix_id="blosum62",
            open_penalty=10,
            extend_penalty=1,
            use_traceback=True,
        )

        self.assertEqual(result, (["seq1", "seq2"], 25.0))
        mock_get_traceback.assert_called_once()

    @patch("parasail.nw_stats_striped_32")
    def test_get_stats_score(self, mock_nw_stats):
        """Test get_stats_score function."""
        mock_result = Mock()
        mock_result.matches = 3
        mock_result.length = 5
        mock_nw_stats.return_value = mock_result

        mock_matrix = Mock()

        score = get_stats_score("ATCG", "ATCG", 8, 1, mock_matrix)

        # num_ungapped_cols = 4 + 4 - 5 = 3
        # identity_score_percent = (3 / 3) * 100 = 100
        # distance_score = 100 - 100 = 0
        self.assertEqual(score, 0.0)

    def test_get_similarity(self):
        """Test get_similarity function."""
        # Perfect match
        similarity = get_similarity("ATCG", "ATCG")
        self.assertEqual(similarity, 0.0)

        # One mismatch
        similarity = get_similarity("ATCG", "ATCC")
        self.assertEqual(similarity, 0.25)

        # With gaps
        similarity = get_similarity("AT-G", "AT-C")
        self.assertEqual(similarity, 1.0 / 3.0)  # 1 mismatch out of 3 total positions

    def test_get_similarity_unequal_length(self):
        """Test get_similarity with unequal length sequences."""
        with self.assertRaises(ValueError):
            get_similarity("ATCG", "ATCGG")

    @patch("parasail.nw_trace")
    @patch("workflow.analyze.parasail.get_similarity")
    def test_get_traceback_score(self, mock_get_similarity, mock_nw_trace):
        """Test get_traceback_score function."""
        mock_result = Mock()
        mock_result.traceback.query = "ATCG"
        mock_result.traceback.ref = "ATCC"
        mock_nw_trace.return_value = mock_result

        mock_get_similarity.return_value = 0.25

        mock_matrix = Mock()

        score = get_traceback_score("ATCG", "ATCC", 10, 1, mock_matrix)

        self.assertEqual(score, 25.0)
        mock_get_similarity.assert_called_once_with("ATCG", "ATCC")

    @patch("parasail.nw_trace")
    def test_get_traceback_score_exception(self, mock_nw_trace):
        """Test get_traceback_score with parasail exception."""
        mock_nw_trace.side_effect = Exception("Parasail error")

        mock_matrix = Mock()

        with self.assertRaises(Exception) as context:
            get_traceback_score("ATCG", "ATCC", 10, 1, mock_matrix)

        self.assertEqual(str(context.exception), "PARASAIL_TRACEBACK")


if __name__ == "__main__":
    unittest.main()
