import unittest
import tempfile
import os
import json
import shutil
from unittest.mock import patch
import pandas as pd
from pandas import DataFrame

from export_utils import (
    save_matrix_to_csv,
    save_cols_to_csv,
    save_stats_to_csv,
    save_seq_dict_to_json,
    save_run_settings_to_json,
)


class TestExportUtils(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        # Sample data for testing
        self.sample_df = DataFrame(
            [[0, 10, 25], [10, 0, 15], [25, 15, 0]],
            index=["SeqA", "SeqB", "SeqC"],
            columns=["SeqA", "SeqB", "SeqC"],
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch("export_utils.to_triangle")
    def test_save_matrix_to_csv(self, mock_to_triangle):
        matrix_path = os.path.join(self.temp_dir, "matrix.csv")
        triangle_path = os.path.join(self.temp_dir, "triangle.csv")

        # Mock triangle result
        mock_triangle = DataFrame(
            [[0, None, None], [10, 0, None], [25, 15, 0]],
            index=["SeqA", "SeqB", "SeqC"],
            columns=["SeqA", "SeqB", "SeqC"],
        )
        mock_to_triangle.return_value = mock_triangle

        save_matrix_to_csv(self.sample_df, matrix_path, triangle_path)

        # Verify both files were created
        self.assertTrue(os.path.exists(matrix_path))
        self.assertTrue(os.path.exists(triangle_path))

        # Verify matrix file content (read without header since export writes without header)
        matrix_result = pd.read_csv(matrix_path, index_col=0, header=None)
        matrix_result.columns = self.sample_df.columns
        pd.testing.assert_frame_equal(matrix_result, self.sample_df, check_names=False)

        # Verify triangle function was called
        mock_to_triangle.assert_called_once_with(self.sample_df)

    def test_save_cols_to_csv(self):
        cols_path = os.path.join(self.temp_dir, "columns.csv")

        save_cols_to_csv(self.sample_df, cols_path)

        self.assertTrue(os.path.exists(cols_path))

        # Read and verify the columnar format
        result_df = pd.read_csv(cols_path)
        expected_columns = ["First Sequence", "Second Sequence", "Identity Score"]
        self.assertEqual(list(result_df.columns), expected_columns)

        # Should only have lower triangle comparisons
        self.assertEqual(len(result_df), 3)  # 3 choose 2 = 3 pairs

        # Check specific values (convert distance back to similarity)
        first_row = result_df.iloc[0]
        self.assertEqual(first_row["First Sequence"], "SeqB")
        self.assertEqual(first_row["Second Sequence"], "SeqA")
        self.assertEqual(first_row["Identity Score"], 90.0)  # 100 - 10

    def test_save_stats_to_csv(self):
        stats_path = os.path.join(self.temp_dir, "stats.csv")

        seq_stats = {"SeqA": [45.2, 1000], "SeqB": [52.1, 1200], "SeqC": [38.9, 800]}

        save_stats_to_csv(seq_stats, stats_path)

        self.assertTrue(os.path.exists(stats_path))

        # Read and verify content
        result_df = pd.read_csv(stats_path)
        expected_columns = ["Sequence", "GC %", "Sequence Length"]
        self.assertEqual(list(result_df.columns), expected_columns)
        self.assertEqual(len(result_df), 3)

        # Check specific values
        seq_a_row = result_df[result_df["Sequence"] == "SeqA"].iloc[0]
        self.assertEqual(seq_a_row["GC %"], 45.2)
        self.assertEqual(seq_a_row["Sequence Length"], 1000)

    def test_save_seq_dict_to_json(self):
        json_path = os.path.join(self.temp_dir, "sequences.json")

        seq_dict = {
            "SeqA": "ATCGATCGATCG",
            "SeqB": "GCTAGCTAGCTA",
            "SeqC": "TTAACCGGTTAA",
        }

        save_seq_dict_to_json(seq_dict, json_path)

        self.assertTrue(os.path.exists(json_path))

        # Read and verify content
        with open(json_path, "r") as f:
            result = json.load(f)

        self.assertEqual(result, seq_dict)

        # Verify JSON formatting (indented)
        with open(json_path, "r") as f:
            content = f.read()

        self.assertIn("    ", content)  # Should have indentation

    def test_save_run_settings_to_json(self):
        settings_path = os.path.join(self.temp_dir, "settings.json")

        run_settings = {
            "algorithm": "lzani",
            "scoring_matrix": "BLOSUM62",
            "gap_open": -10,
            "gap_extend": -1,
            "threads": 4,
        }

        save_run_settings_to_json(run_settings, settings_path)

        self.assertTrue(os.path.exists(settings_path))

        # Read and verify content
        with open(settings_path, "r") as f:
            result = json.load(f)

        self.assertEqual(result, run_settings)

        # Verify JSON formatting (indented)
        with open(settings_path, "r") as f:
            content = f.read()

        self.assertIn("    ", content)  # Should have indentation

    def test_save_cols_to_csv_single_sequence(self):
        # Edge case: single sequence (no pairs)
        single_df = DataFrame([[0]], index=["SeqA"], columns=["SeqA"])
        cols_path = os.path.join(self.temp_dir, "single_columns.csv")

        save_cols_to_csv(single_df, cols_path)

        self.assertTrue(os.path.exists(cols_path))

        result_df = pd.read_csv(cols_path)
        self.assertEqual(len(result_df), 0)  # No pairs for single sequence

    def test_save_cols_to_csv_large_matrix(self):
        # Test with larger matrix
        n = 5
        large_df = DataFrame(
            [[abs(i - j) * 10 for j in range(n)] for i in range(n)],
            index=[f"Seq{i}" for i in range(n)],
            columns=[f"Seq{i}" for i in range(n)],
        )
        cols_path = os.path.join(self.temp_dir, "large_columns.csv")

        save_cols_to_csv(large_df, cols_path)

        self.assertTrue(os.path.exists(cols_path))

        result_df = pd.read_csv(cols_path)
        # Should have n choose 2 pairs = 5*4/2 = 10 pairs
        self.assertEqual(len(result_df), 10)

    def test_save_stats_to_csv_empty_stats(self):
        stats_path = os.path.join(self.temp_dir, "empty_stats.csv")

        save_stats_to_csv({}, stats_path)

        self.assertTrue(os.path.exists(stats_path))

        result_df = pd.read_csv(stats_path)
        self.assertEqual(len(result_df), 0)
        self.assertEqual(
            list(result_df.columns), ["Sequence", "GC %", "Sequence Length"]
        )

    def test_save_seq_dict_to_json_empty_dict(self):
        json_path = os.path.join(self.temp_dir, "empty.json")

        save_seq_dict_to_json({}, json_path)

        self.assertTrue(os.path.exists(json_path))

        with open(json_path, "r") as f:
            result = json.load(f)

        self.assertEqual(result, {})


if __name__ == "__main__":
    unittest.main()
