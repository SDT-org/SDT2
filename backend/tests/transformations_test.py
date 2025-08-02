import unittest
import numpy as np
import tempfile
import os
from pandas import DataFrame

from transformations import (
    to_triangle,
    similarity_triangle_to_matrix,
    read_csv_file,
    read_csv_matrix,
    read_stats_csv,
    read_columns_csv,
    lzani_tsv_to_distance_matrix,
)


class TestTransformations(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.sample_matrix = np.array(
            [[0.0, 10.0, 25.0], [10.0, 0.0, 15.0], [25.0, 15.0, 0.0]], dtype=float
        )

        self.sample_df = DataFrame(
            self.sample_matrix, index=["A", "B", "C"], columns=["A", "B", "C"]
        ).astype(float)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_to_triangle_dataframe_default(self):
        result = to_triangle(self.sample_df)

        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        # Check conversion to similarity and triangular form
        expected = np.array(
            [[100.0, np.nan, np.nan], [90.0, 100.0, np.nan], [75.0, 85.0, 100.0]]
        )
        np.testing.assert_array_equal(result.values, expected)
        self.assertEqual(list(result.index), ["A", "B", "C"])
        self.assertEqual(list(result.columns), ["A", "B", "C"])

    def test_to_triangle_numpy_array(self):
        result = to_triangle(self.sample_matrix)

        self.assertIsInstance(result, np.ndarray)

        expected = np.array(
            [[100.0, np.nan, np.nan], [90.0, 100.0, np.nan], [75.0, 85.0, 100.0]]
        )

        np.testing.assert_array_equal(result, expected)

    def test_to_triangle_no_similarity_conversion(self):
        result = to_triangle(self.sample_df, convert_to_similarity=False)

        expected = np.array(
            [[0.0, np.nan, np.nan], [10.0, 0.0, np.nan], [25.0, 15.0, 0.0]]
        )

        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        np.testing.assert_array_equal(result.values, expected)

    def test_to_triangle_custom_fill_value(self):
        result = to_triangle(self.sample_df, fill_value=-1)

        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        expected = np.array([[100.0, -1, -1], [90.0, 100.0, -1], [75.0, 85.0, 100.0]])

        np.testing.assert_array_equal(result.values, expected)

    def test_to_triangle_none_fill_value(self):
        result = to_triangle(self.sample_df, fill_value=None)

        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        # Upper triangle should be None
        self.assertIsNone(result.iloc[0, 1])
        self.assertIsNone(result.iloc[0, 2])
        self.assertIsNone(result.iloc[1, 2])

        # Lower triangle should have values
        self.assertEqual(result.iloc[1, 0], 90.0)
        self.assertEqual(result.iloc[2, 0], 75.0)
        self.assertEqual(result.iloc[2, 1], 85.0)

    def test_similarity_triangle_to_matrix(self):
        triangle_matrix = DataFrame(
            [[0, np.nan, np.nan], [10, 0, np.nan], [25, 15, 0]],
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )

        result = similarity_triangle_to_matrix(triangle_matrix)

        expected = np.array(
            [[100.0, 90.0, 75.0], [90.0, 100.0, 85.0], [75.0, 85.0, 100.0]]
        )

        np.testing.assert_array_equal(result.values, expected)
        self.assertEqual(list(result.index), ["A", "B", "C"])

    def test_read_csv_file_default(self):
        csv_path = os.path.join(self.temp_dir, "test.csv")
        test_data = "1,2,3\n4,5,6"  # No header for default case

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_csv_file(csv_path)
        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        # Default parameters: no header, no index_col - uses variable column handling
        self.assertEqual(result.shape, (2, 3))
        self.assertEqual(result.iloc[0, 0], 1)
        self.assertEqual(result.iloc[1, 2], 6)

    def test_read_csv_file_as_list(self):
        csv_path = os.path.join(self.temp_dir, "test.csv")
        test_data = "1,2,3\n4,5,6"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_csv_file(csv_path, header=None, as_list=True)

        expected = [[1, 2, 3], [4, 5, 6]]
        self.assertEqual(result, expected)

    def test_read_csv_file_extract_column(self):
        csv_path = os.path.join(self.temp_dir, "test.csv")
        test_data = "A,B,C\n1,2,3\n4,5,6"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_csv_file(csv_path, header=0, extract_column="B")

        expected = [2, 5]
        self.assertEqual(result, expected)

    def test_read_csv_file_variable_columns(self):
        csv_path = os.path.join(self.temp_dir, "matrix.csv")
        test_data = "A,0,10\nB,10,0,15\nC,25,15,0"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_csv_file(csv_path, header=None, index_col=0)
        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        self.assertEqual(result.shape[1], 3)  # Should handle variable columns
        self.assertEqual(list(result.index), ["A", "B", "C"])

    def test_read_csv_matrix(self):
        csv_path = os.path.join(self.temp_dir, "matrix.csv")
        test_data = "A,0,10,25\nB,10,0,15\nC,25,15,0"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_csv_matrix(csv_path)
        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        self.assertEqual(list(result.index), ["A", "B", "C"])
        self.assertEqual(result.shape, (3, 3))

    def test_read_stats_csv(self):
        csv_path = os.path.join(self.temp_dir, "stats.csv")
        test_data = "Sequence,GC %,Length\nSeq1,45.2,1000\nSeq2,52.1,1200"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_stats_csv(csv_path)
        if not isinstance(result, DataFrame):
            raise TypeError("Expected result to be a DataFrame")

        self.assertEqual(list(result.columns), ["Sequence", "GC %", "Length"])
        self.assertEqual(len(result), 2)

    def test_read_columns_csv(self):
        csv_path = os.path.join(self.temp_dir, "columns.csv")
        test_data = "First Sequence,Second Sequence,Identity Score\nSeq1,Seq2,85.5\nSeq1,Seq3,92.1"

        with open(csv_path, "w") as f:
            f.write(test_data)

        result = read_columns_csv(csv_path)

        expected = [["Seq1", "Seq2", 85.5], ["Seq1", "Seq3", 92.1]]
        self.assertEqual(result, expected)

    def test_read_tsv_file(self):
        tsv_path = os.path.join(self.temp_dir, "test.tsv")
        test_data = "id\tvalue\nA\t100\nB\t200"

        with open(tsv_path, "w") as f:
            f.write(test_data)

        # read_tsv_file calls read_csv_file with sep="\t" and extract_column
        # We need to pass the right parameters through read_csv_file
        from transformations import read_csv_file

        result = read_csv_file(tsv_path, sep="\t", header=0, extract_column="id")

        expected = ["A", "B"]
        self.assertEqual(result, expected)

    def test_lzani_tsv_to_distance_matrix_basic(self):
        # Create test files
        results_tsv = os.path.join(self.temp_dir, "results.tsv")
        ids_tsv = os.path.join(self.temp_dir, "ids.tsv")

        results_data = "query\treference\tani\nSeq1\tSeq2\t0.85\nSeq2\tSeq1\t0.83\nSeq1\tSeq3\t0.92"
        ids_data = "id\nSeq1\nSeq2\nSeq3"

        with open(results_tsv, "w") as f:
            f.write(results_data)
        with open(ids_tsv, "w") as f:
            f.write(ids_data)

        distance_matrix, ids = lzani_tsv_to_distance_matrix(results_tsv, ids_tsv)

        self.assertEqual(ids, ["Seq1", "Seq2", "Seq3"])
        self.assertEqual(distance_matrix.shape, (3, 3))

        # Check diagonal is 0 (distance from self)
        np.testing.assert_array_equal(np.diag(distance_matrix), [0, 0, 0])

        # Check symmetric averaging
        self.assertAlmostEqual(distance_matrix[0, 1], 100 - ((85 + 83) / 2), places=1)

        # Check one-way alignment
        self.assertAlmostEqual(distance_matrix[0, 2], 100 - 92, places=1)

    def test_lzani_tsv_to_distance_matrix_empty_results(self):
        results_tsv = os.path.join(self.temp_dir, "empty_results.tsv")
        ids_tsv = os.path.join(self.temp_dir, "ids.tsv")

        # Add one minimal result to avoid empty DataFrame issue
        results_data = "query\treference\tani\nSeq1\tSeq2\t0.0"  # Very low score, will be treated as 0
        ids_data = "id\nSeq1\nSeq2"

        with open(results_tsv, "w") as f:
            f.write(results_data)
        with open(ids_tsv, "w") as f:
            f.write(ids_data)

        distance_matrix, ids = lzani_tsv_to_distance_matrix(results_tsv, ids_tsv)

        self.assertEqual(ids, ["Seq1", "Seq2"])
        self.assertEqual(distance_matrix.shape, (2, 2))

        # All values should be 100 (maximum distance) except diagonal which is 0
        # The 0.0 score gets converted to distance 100
        expected = np.array([[0.0, 100.0], [100.0, 0.0]])
        np.testing.assert_array_equal(distance_matrix, expected)

    def test_lzani_tsv_to_distance_matrix_threshold_filtering(self):
        results_tsv = os.path.join(self.temp_dir, "results.tsv")
        ids_tsv = os.path.join(self.temp_dir, "ids.tsv")

        # Low scores should be set to 0, then converted to distance 100
        results_data = "query\treference\tani\nSeq1\tSeq2\t0.005\nSeq2\tSeq1\t0.008"
        ids_data = "id\nSeq1\nSeq2"

        with open(results_tsv, "w") as f:
            f.write(results_data)
        with open(ids_tsv, "w") as f:
            f.write(ids_data)

        distance_matrix, ids = lzani_tsv_to_distance_matrix(results_tsv, ids_tsv)

        # Both scores are below 1%, so distance should be 100
        self.assertEqual(distance_matrix[0, 1], 100)
        self.assertEqual(distance_matrix[1, 0], 100)


if __name__ == "__main__":
    unittest.main()
