import unittest
import numpy as np
import tempfile
import json
import shutil
from unittest.mock import patch, MagicMock, mock_open
from pandas import DataFrame

from workflow.cluster import (
    run,
    calculate_linkage,
    get_mds_coords,
    get_linkage,
    get_linkage_method_order,
    get_clusters_dataframe,
    export,
    get_cluster_data_dict,
    cluster_data_to_dataframe,
    order_clusters_sequentially,
)
from workflow.models import RunSettings, WorkflowResult


class TestCluster(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        # Sample distance matrix for testing
        self.sample_distance_matrix = np.array(
            [
                [0.0, 10.0, 25.0, 30.0],
                [10.0, 0.0, 15.0, 35.0],
                [25.0, 15.0, 0.0, 20.0],
                [30.0, 35.0, 20.0, 0.0],
            ]
        )

        self.sample_ids = ["SeqA", "SeqB", "SeqC", "SeqD"]

        # Create mock settings and result objects
        self.mock_settings = MagicMock(spec=RunSettings)
        self.mock_settings.cluster_method = "single"

        self.mock_result = MagicMock(spec=WorkflowResult)
        self.mock_result.distance_matrix = self.sample_distance_matrix
        self.mock_result.ordered_ids = self.sample_ids

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_run_with_valid_matrix(self):
        from workflow.models import WorkflowResult

        real_result = WorkflowResult(
            seq_dict={},
            ordered_ids=self.sample_ids,
            reordered_ids=[],
            min_score=0.0,
            max_sequence_length=1000,
            distance_matrix=self.sample_distance_matrix,
            similarity_matrix=DataFrame(),
            warnings=[],
            error=None,
            is_aa=False,
        )

        result = run(real_result, self.mock_settings)

        # Should return a result with reordered matrix and ids
        self.assertIsNotNone(result.distance_matrix)
        self.assertIsNotNone(result.reordered_ids)

        # Matrix should maintain same shape
        self.assertEqual(result.distance_matrix.shape, (4, 4))

    def test_run_with_none_matrix(self):
        self.mock_result.distance_matrix = None

        result = run(self.mock_result, self.mock_settings)

        # Should return original result unchanged
        self.assertEqual(result, self.mock_result)

    def test_run_with_pandas_dataframe(self):
        df_matrix = DataFrame(
            self.sample_distance_matrix, index=self.sample_ids, columns=self.sample_ids
        )
        self.mock_result.distance_matrix = df_matrix

        result = run(self.mock_result, self.mock_settings)

        self.assertIsNotNone(result.distance_matrix)
        self.assertIsNotNone(result.reordered_ids)

    def test_run_with_no_ordered_ids(self):
        from workflow.models import WorkflowResult

        real_result = WorkflowResult(
            seq_dict={},
            ordered_ids=None,
            reordered_ids=[],
            min_score=0.0,
            max_sequence_length=1000,
            distance_matrix=self.sample_distance_matrix,
            similarity_matrix=DataFrame(),
            warnings=[],
            error=None,
            is_aa=False,
        )

        result = run(real_result, self.mock_settings)

        self.assertIsNotNone(result.reordered_ids)
        # Should use indices instead of sequence IDs
        self.assertIsInstance(result.reordered_ids[0], (int, np.integer))

    def test_calculate_linkage_single_method(self):
        linkage_matrix = calculate_linkage(self.sample_distance_matrix, "single")

        # Linkage matrix should have n-1 rows and 4 columns
        self.assertEqual(linkage_matrix.shape, (3, 4))

        # All values should be finite
        self.assertTrue(np.all(np.isfinite(linkage_matrix)))

    def test_calculate_linkage_ward_method(self):
        # Ward method requires MDS transformation
        linkage_matrix = calculate_linkage(self.sample_distance_matrix, "ward")

        self.assertEqual(linkage_matrix.shape, (3, 4))
        self.assertTrue(np.all(np.isfinite(linkage_matrix)))

    def test_calculate_linkage_ward_method_mds_failure(self):
        # Test that when get_mds_coords returns None, calculate_linkage raises ValueError
        with patch("workflow.cluster.get_mds_coords", return_value=None):
            with self.assertRaises(ValueError) as cm:
                calculate_linkage(self.sample_distance_matrix, "ward")

            self.assertIn("MDS failed to compute coordinates", str(cm.exception))

    def test_get_mds_coords(self):
        coords = get_mds_coords(self.sample_distance_matrix)

        # Should return 2D coordinates for each point
        self.assertEqual(coords.shape, (4, 2))

        # All coordinates should be finite
        self.assertTrue(np.all(np.isfinite(coords)))

    def test_get_linkage(self):
        # This is a wrapper function for calculate_linkage
        linkage_matrix = get_linkage(self.sample_distance_matrix, "complete")

        self.assertEqual(linkage_matrix.shape, (3, 4))
        self.assertTrue(np.all(np.isfinite(linkage_matrix)))

    def test_get_linkage_method_order_no_threshold(self):
        order = get_linkage_method_order(
            self.sample_distance_matrix, "single", self.sample_ids
        )

        # Should return all sequence IDs in some order
        self.assertEqual(len(order), 4)
        self.assertEqual(set(order), set(self.sample_ids))

    def test_get_linkage_method_order_with_threshold(self):
        order = get_linkage_method_order(
            self.sample_distance_matrix, "single", self.sample_ids, threshold=80
        )

        # Should return all sequence IDs
        self.assertEqual(len(order), 4)
        self.assertEqual(set(order), set(self.sample_ids))

    def test_get_clusters_dataframe(self):
        df = get_clusters_dataframe(
            self.sample_distance_matrix, "single", 80, self.sample_ids
        )

        # Should return DataFrame with ID and Cluster columns
        self.assertIsInstance(df, DataFrame)
        self.assertEqual(len(df), 4)
        self.assertIn("ID", df.columns)
        self.assertTrue(any("Cluster" in col for col in df.columns))

    @patch("workflow.cluster.read_csv_matrix")
    @patch("builtins.open", new_callable=mock_open, read_data="SeqA,0,10\nSeqB,10,0")
    @patch("os.path.exists")
    @patch("os.makedirs")
    def test_export_basic(self, mock_makedirs, mock_exists, mock_file, mock_read_csv):
        # Mock the CSV matrix reading
        mock_df = DataFrame(
            [[0, 10], [10, 0]], index=["SeqA", "SeqB"], columns=["SeqA", "SeqB"]
        )
        mock_read_csv.return_value = mock_df

        # Mock sequence dictionary file
        mock_exists.return_value = True
        seq_dict = {"SeqA": "ATCG", "SeqB": "GCTA"}

        with patch("builtins.open", mock_open(read_data=json.dumps(seq_dict))):
            with patch("json.load", return_value=seq_dict):
                with patch("workflow.cluster.SeqIO.write") as mock_seqio:
                    result_df = export(
                        "/fake/matrix.csv",
                        self.temp_dir,
                        "/fake/seq_dict.json",
                        80,
                        "single",
                    )

        # Should return a DataFrame
        self.assertIsInstance(result_df, DataFrame)
        # makedirs should be called at least once (joblib may call it multiple times)
        self.assertTrue(mock_makedirs.called)

    @patch("workflow.cluster.read_csv_matrix")
    @patch("builtins.open", new_callable=mock_open, read_data="SeqA,0,10\nSeqB,10,0")
    @patch("os.path.exists")
    @patch("os.makedirs")
    def test_export_no_seq_dict(
        self, mock_makedirs, mock_exists, mock_file, mock_read_csv
    ):
        # Mock the CSV matrix reading
        mock_df = DataFrame(
            [[0, 10], [10, 0]], index=["SeqA", "SeqB"], columns=["SeqA", "SeqB"]
        )
        mock_read_csv.return_value = mock_df

        # Mock sequence dictionary file doesn't exist
        mock_exists.return_value = False

        result_df = export(
            "/fake/matrix.csv", self.temp_dir, "/fake/seq_dict.json", 80, "single"
        )

        # Should still return a DataFrame
        self.assertIsInstance(result_df, DataFrame)
        # makedirs should be called at least once (joblib may call it multiple times)
        self.assertTrue(mock_makedirs.called)

    def test_get_cluster_data_dict(self):
        # Create a simple linkage matrix for testing
        linkage_matrix = calculate_linkage(self.sample_distance_matrix, "single")

        cluster_dict = get_cluster_data_dict(linkage_matrix, 80, self.sample_ids)

        # Should return a dictionary
        self.assertIsInstance(cluster_dict, dict)

        # All sequence IDs should be accounted for
        all_sequences = []
        for cluster_seqs in cluster_dict.values():
            all_sequences.extend(cluster_seqs)
        self.assertEqual(set(all_sequences), set(self.sample_ids))

    def test_cluster_data_to_dataframe(self):
        cluster_data = {0: ["SeqA", "SeqB"], 1: ["SeqC", "SeqD"]}

        df = cluster_data_to_dataframe(cluster_data, 80)

        self.assertIsInstance(df, DataFrame)
        self.assertEqual(len(df), 4)  # 4 sequences total
        self.assertEqual(list(df.columns), ["ID", "Cluster - Threshold: 80"])

        # Check that all sequences are present
        self.assertEqual(set(df["ID"]), {"SeqA", "SeqB", "SeqC", "SeqD"})

        # Check cluster assignments
        seqa_cluster = df[df["ID"] == "SeqA"]["Cluster - Threshold: 80"].iloc[0]
        seqb_cluster = df[df["ID"] == "SeqB"]["Cluster - Threshold: 80"].iloc[0]
        self.assertEqual(seqa_cluster, seqb_cluster)  # Should be in same cluster

    def test_order_clusters_sequentially(self):
        # Test cluster label reordering
        clusters = [3, 1, 3, 1, 7, 7]

        reordered = order_clusters_sequentially(clusters)

        # Should maintain same length
        self.assertEqual(len(reordered), 6)

        # Should start from 1 and be sequential
        unique_clusters = sorted(set(reordered))
        expected_unique = list(range(1, len(set(clusters)) + 1))
        self.assertEqual(unique_clusters, expected_unique)

        # Same original cluster should have same new label
        self.assertEqual(reordered[0], reordered[2])  # Both were 3
        self.assertEqual(reordered[1], reordered[3])  # Both were 1
        self.assertEqual(reordered[4], reordered[5])  # Both were 7

    def test_order_clusters_sequentially_single_cluster(self):
        clusters = [5, 5, 5]

        reordered = order_clusters_sequentially(clusters)

        self.assertEqual(reordered, [1, 1, 1])

    def test_order_clusters_sequentially_already_sequential(self):
        clusters = [1, 1, 2, 2, 3]

        reordered = order_clusters_sequentially(clusters)

        self.assertEqual(reordered, [1, 1, 2, 2, 3])

    def test_calculate_linkage_all_methods(self):
        """Test calculate_linkage with all available methods"""
        similarity_array = [[100, 90, 75], [90, 100, 90], [75, 90, 100]]
        test_array = [[100 - val for val in row] for row in similarity_array]

        expected = dict(
            single=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 10.0, 3.0]],
            complete=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 25.0, 3.0]],
            average=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 17.5, 3.0]],
            weighted=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 17.5, 3.0]],
            centroid=[
                [1.0, 2.0, 11.672874240667145, 2.0],
                [0.0, 3.0, 17.49576226800019, 3.0],
            ],
            median=[
                [1.0, 2.0, 11.672874240667145, 2.0],
                [0.0, 3.0, 17.49576226800019, 3.0],
            ],
            ward=[
                [1.0, 2.0, 11.672874240667145, 2.0],
                [0.0, 3.0, 20.202366110215213, 3.0],
            ],
        )

        for method in expected:
            result = calculate_linkage(np.array(test_array), method)
            np.testing.assert_array_equal(result, expected[method])


if __name__ == "__main__":
    unittest.main()
