import unittest
from unittest.mock import Mock, patch
import numpy as np

from workflow.postprocess import run, get_seq_stats
from mime_setup import register_mimetypes

register_mimetypes()


class TestPostprocess(unittest.TestCase):
    """Test postprocess.py module."""

    def setUp(self):
        self.mock_result = Mock()
        self.mock_result.seq_dict = {
            "seq1": "ATCGATCG",
            "seq2": "GCTAGCTA",
            "seq3": "AAAATTTT",
        }
        self.mock_result.is_aa = False
        self.mock_result.ordered_ids = ["seq1", "seq2", "seq3"]
        self.mock_result.reordered_ids = []
        self.mock_result.distance_matrix = np.array(
            [[0.0, 25.0, 50.0], [25.0, 0.0, 30.0], [50.0, 30.0, 0.0]]
        )

        self.mock_doc_paths = Mock()
        self.mock_doc_paths.columns = "/tmp/columns.csv"
        self.mock_doc_paths.stats = "/tmp/stats.csv"
        self.mock_doc_paths.matrix = "/tmp/matrix.csv"
        self.mock_doc_paths.triangle = "/tmp/triangle.csv"
        self.mock_doc_paths.seq_dict = "/tmp/seq_dict.json"

        self.mock_settings = Mock()
        self.mock_settings.doc_paths = self.mock_doc_paths

    @patch("workflow.postprocess.save_seq_dict_to_json")
    @patch("workflow.postprocess.save_matrix_to_csv")
    @patch("workflow.postprocess.save_stats_to_csv")
    @patch("workflow.postprocess.save_cols_to_csv")
    def test_run_with_ordered_ids(
        self, mock_save_cols, mock_save_stats, mock_save_matrix, mock_save_seq_dict
    ):
        """Test postprocess run with ordered_ids (no clustering)."""
        result = run(self.mock_result, self.mock_settings)

        # Verify all save functions were called
        mock_save_cols.assert_called_once()
        mock_save_stats.assert_called_once()
        mock_save_matrix.assert_called_once()
        mock_save_seq_dict.assert_called_once()

        # Verify the result is returned unchanged
        self.assertEqual(result, self.mock_result)

    @patch("workflow.postprocess.save_seq_dict_to_json")
    @patch("workflow.postprocess.save_matrix_to_csv")
    @patch("workflow.postprocess.save_stats_to_csv")
    @patch("workflow.postprocess.save_cols_to_csv")
    def test_run_with_reordered_ids(
        self, mock_save_cols, mock_save_stats, mock_save_matrix, mock_save_seq_dict
    ):
        """Test postprocess run with reordered_ids (clustering was performed)."""
        self.mock_result.reordered_ids = [
            "seq3",
            "seq1",
            "seq2",
        ]  # Different order from clustering

        result = run(self.mock_result, self.mock_settings)

        # Verify all save functions were called
        mock_save_cols.assert_called_once()
        mock_save_stats.assert_called_once()
        mock_save_matrix.assert_called_once()
        mock_save_seq_dict.assert_called_once()

        # Verify the result is returned unchanged
        self.assertEqual(result, self.mock_result)

    @patch("workflow.postprocess.get_seq_stats")
    @patch("workflow.postprocess.save_seq_dict_to_json")
    @patch("workflow.postprocess.save_matrix_to_csv")
    @patch("workflow.postprocess.save_stats_to_csv")
    @patch("workflow.postprocess.save_cols_to_csv")
    def test_run_calls_get_seq_stats(
        self,
        mock_save_cols,
        mock_save_stats,
        mock_save_matrix,
        mock_save_seq_dict,
        mock_get_seq_stats,
    ):
        """Test that run calls get_seq_stats with correct parameters."""
        mock_stats = {"seq1": [50.0, 8], "seq2": [62.5, 8], "seq3": [25.0, 8]}
        mock_get_seq_stats.return_value = mock_stats

        result = run(self.mock_result, self.mock_settings)

        # Verify get_seq_stats was called with correct parameters
        mock_get_seq_stats.assert_called_once_with(
            self.mock_result.seq_dict, self.mock_result.is_aa
        )

        # Verify save_stats_to_csv was called with the stats from get_seq_stats
        mock_save_stats.assert_called_once_with(mock_stats, self.mock_doc_paths.stats)

    def test_get_seq_stats_dna_sequences(self):
        """Test get_seq_stats with DNA sequences."""
        seq_dict = {
            "seq1": "ATCG",  # 50% GC
            "seq2": "AAAA",  # 0% GC
            "seq3": "GCGC",  # 100% GC
            "seq4": "ATCGATCG",  # 50% GC
        }
        is_aa = False

        with patch("workflow.postprocess.gc_fraction") as mock_gc_fraction:
            # Mock gc_fraction to return known values
            mock_gc_fraction.side_effect = [0.5, 0.0, 1.0, 0.5]

            stats = get_seq_stats(seq_dict, is_aa)

        expected_stats = {
            "seq1": [50.0, 4],  # 50% GC, length 4
            "seq2": [0.0, 4],  # 0% GC, length 4
            "seq3": [100.0, 4],  # 100% GC, length 4
            "seq4": [50.0, 8],  # 50% GC, length 8
        }

        self.assertEqual(stats, expected_stats)

        # Verify gc_fraction was called for each sequence
        self.assertEqual(mock_gc_fraction.call_count, 4)

    def test_get_seq_stats_amino_acid_sequences(self):
        """Test get_seq_stats with amino acid sequences."""
        seq_dict = {
            "prot1": "MKTFFLLLLFTI",
            "prot2": "ARNDCQEGHILKMFPSTWYV",
            "prot3": "AAA",
        }
        is_aa = True

        stats = get_seq_stats(seq_dict, is_aa)

        expected_stats = {
            "prot1": [0.0, 12],  # 0% GC (AA sequences), length 12
            "prot2": [0.0, 20],  # 0% GC (AA sequences), length 20
            "prot3": [0.0, 3],  # 0% GC (AA sequences), length 3
        }

        self.assertEqual(stats, expected_stats)

    def test_get_seq_stats_empty_dict(self):
        """Test get_seq_stats with empty sequence dictionary."""
        seq_dict = {}
        is_aa = False

        stats = get_seq_stats(seq_dict, is_aa)

        self.assertEqual(stats, {})

    def test_get_seq_stats_single_sequence(self):
        """Test get_seq_stats with single sequence."""
        seq_dict = {"only_seq": "ATCGATCG"}
        is_aa = False

        with patch("workflow.postprocess.gc_fraction", return_value=0.5):
            stats = get_seq_stats(seq_dict, is_aa)

        expected_stats = {"only_seq": [50.0, 8]}  # 50% GC, length 8

        self.assertEqual(stats, expected_stats)

    @patch("workflow.postprocess.gc_fraction")
    def test_get_seq_stats_gc_fraction_rounding(self, mock_gc_fraction):
        """Test that GC fraction is properly rounded to 2 decimal places."""
        seq_dict = {"test_seq": "ATCGATCGATCG"}
        is_aa = False

        # Return a value that needs rounding
        mock_gc_fraction.return_value = 0.583333333

        stats = get_seq_stats(seq_dict, is_aa)

        # The actual calculation is: round(0.583333333, 2) * 100 = 0.58 * 100 = 58.0
        # But due to floating point precision, we need to check the rounded result
        expected_gc_percent = round(0.583333333, 2) * 100
        expected_stats = {
            "test_seq": [expected_gc_percent, 12]  # GC (rounded), length 12
        }

        self.assertEqual(stats, expected_stats)


if __name__ == "__main__":
    unittest.main()
