import os
import sys
import unittest
import numpy
from unittest.mock import patch

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../src/"))

import cluster

# Unwrap the original function to bypass joblib caching
original_func = cluster.calculate_linkage.__wrapped__


class TestCluster(unittest.TestCase):
    def setUp(self):
        # Create a patcher that replaces the cached function with the unwrapped version
        self.patcher = patch("cluster.calculate_linkage", original_func)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()

    def test_calculate_linkage(self):
        test_array = [[100, 90, 75], [90, 100, 90], [75, 90, 100]]
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


if __name__ == "__main__":
    unittest.main()
