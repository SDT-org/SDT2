import os
import sys
import unittest
import numpy

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../src/"))

from cluster import calculate_linkage

EXPECTED = dict(
    single=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 10.0, 3.0]],
    complete=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 25.0, 3.0]],
    average=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 17.5, 3.0]],
    weighted=[[0.0, 1.0, 10.0, 2.0], [2.0, 3.0, 17.5, 3.0]],
    centroid=[[1.0, 2.0, 11.672874240667145, 2.0], [0.0, 3.0, 17.49576226800019, 3.0]],
    median=[[1.0, 2.0, 11.672874240667145, 2.0], [0.0, 3.0, 17.49576226800019, 3.0]],
    ward=[[1.0, 2.0, 11.672874240667145, 2.0], [0.0, 3.0, 20.202366110215213, 3.0]],
)


class TestCluster(unittest.TestCase):
    def test_calculate_linkage(self):
        for method in EXPECTED:
            result = calculate_linkage(
                numpy.array([[100, 0, 0], [90, 100, 0], [75, 90, 100]]), method
            )
            numpy.testing.assert_array_equal(result, EXPECTED[method])

    def test_get_linkage(self):
        pass

    def test_get_linkage_method_order(self):
        pass


if __name__ == "__main__":
    unittest.main()
