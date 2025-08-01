import unittest
import numpy
from unittest.mock import patch

# Path setup handled by PYTHONPATH in package.json

import workflow.cluster as cluster

# Unwrap the original function to bypass joblib caching
original_func = cluster.calculate_linkage.__wrapped__


class TestCluster(unittest.TestCase):
    def setUp(self):
        # Create a patcher that replaces the cached function with the unwrapped version
        self.patcher = patch("workflow.cluster.calculate_linkage", original_func)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_calculate_linkage(self):
        # Convert similarity matrix to distance matrix (diagonal should be 0)
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
            result = cluster.calculate_linkage(numpy.array(test_array), method)
            numpy.testing.assert_array_equal(result, expected[method])
