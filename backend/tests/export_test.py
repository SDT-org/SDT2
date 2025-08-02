import unittest
import tempfile
import os
import shutil
import base64
from unittest.mock import patch, MagicMock

from export import (
    save_image_from_api,
    build_source_target_pairs,
    do_export,
    EXPORTABLE_DATA_KEYS,
    MATRIX_ONLY_EXPORTABLE_DATA_KEYS,
)
from document_state import DocState


class TestExport(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.doc = MagicMock(spec=DocState)
        self.doc.tempdir_path = self.temp_dir

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    @patch("export.build_document_paths")
    def test_save_image_from_api_svg(self, mock_build_paths):
        svg_data = '<svg><rect width="100" height="100"/></svg>'
        svg_path = os.path.join(self.temp_dir, "test.svg")

        mock_paths = MagicMock()
        mock_paths.images._asdict.return_value = {
            "heatmap": MagicMock(_asdict=lambda: {"svg": svg_path})
        }
        mock_build_paths.return_value = mock_paths

        save_image_from_api(self.doc, svg_data, "heatmap", "svg")

        self.assertTrue(os.path.exists(svg_path))
        with open(svg_path, "r") as f:
            self.assertEqual(f.read(), svg_data)

    @patch("export.build_document_paths")
    def test_save_image_from_api_base64(self, mock_build_paths):
        test_data = b"test image data"
        base64_data = "data:image/png;base64," + base64.b64encode(test_data).decode()
        png_path = os.path.join(self.temp_dir, "test.png")

        mock_paths = MagicMock()
        mock_paths.images._asdict.return_value = {
            "heatmap": MagicMock(_asdict=lambda: {"png": png_path})
        }
        mock_build_paths.return_value = mock_paths

        save_image_from_api(self.doc, base64_data, "heatmap", "png")

        self.assertTrue(os.path.exists(png_path))
        with open(png_path, "rb") as f:
            self.assertEqual(f.read(), test_data)

    @patch("export.build_document_paths")
    def test_build_source_target_pairs_default(self, mock_build_paths):
        export_path = "/export/path"
        prefix = "test_prefix"

        mock_paths = MagicMock()
        paths_dict = {
            "stats": "/src/stats.csv",
            "matrix": "/src/matrix.csv",
            "triangle": "/src/triangle.csv",
            "columns": "/src/columns.csv",
            "summary": "/src/summary.txt",
            "cluster_dir": "/src/cluster/",
            "images": MagicMock(),
        }

        mock_paths._asdict.return_value = paths_dict

        paths_dict["images"].heatmap = MagicMock(
            _asdict=lambda: {"svg": "/src/heatmap.svg"}
        )
        paths_dict["images"].violin = MagicMock(
            _asdict=lambda: {"svg": "/src/violin.svg"}
        )
        paths_dict["images"].histogram = MagicMock(
            _asdict=lambda: {"svg": "/src/histogram.svg"}
        )
        mock_build_paths.return_value = mock_paths

        with patch("export.IMAGE_KEYS", ["heatmap", "violin", "histogram"]):
            result = build_source_target_pairs("/doc/path", export_path, prefix, "svg")

        expected_data_pairs = [
            ["/src/stats.csv", "/export/path/test_prefix_stats.csv"],
            ["/src/matrix.csv", "/export/path/test_prefix_matrix.csv"],
            ["/src/triangle.csv", "/export/path/test_prefix_triangle.csv"],
            ["/src/columns.csv", "/export/path/test_prefix_columns.csv"],
            ["/src/summary.txt", "/export/path/test_prefix_summary.txt"],
            ["/src/cluster/", "/export/path/test_prefix_cluster_data"],
        ]
        expected_image_pairs = [
            ["/src/heatmap.svg", "/export/path/test_prefix_heatmap.svg"],
            ["/src/violin.svg", "/export/path/test_prefix_violin.svg"],
            ["/src/histogram.svg", "/export/path/test_prefix_histogram.svg"],
        ]

        for pair in expected_data_pairs + expected_image_pairs:
            self.assertIn(pair, result)

    @patch("export.build_document_paths")
    def test_build_source_target_pairs_matrix_only(self, mock_build_paths):
        export_path = "/export/path"
        prefix = "test_prefix"

        mock_paths = MagicMock()
        paths_dict = {
            "matrix": "/src/matrix.csv",
            "columns": "/src/columns.csv",
            "images": MagicMock(),
        }

        mock_paths._asdict.return_value = paths_dict

        paths_dict["images"].heatmap = MagicMock(
            _asdict=lambda: {"svg": "/src/heatmap.svg"}
        )
        mock_build_paths.return_value = mock_paths

        with patch("export.IMAGE_KEYS", ["heatmap"]):
            result = build_source_target_pairs(
                "/doc/path", export_path, prefix, "svg", matrix_only=True
            )

        expected_pairs = [
            ["/src/matrix.csv", "/export/path/test_prefix_matrix.csv"],
            ["/src/columns.csv", "/export/path/test_prefix_columns.csv"],
            ["/src/heatmap.svg", "/export/path/test_prefix_heatmap.svg"],
        ]

        for pair in expected_pairs:
            self.assertIn(pair, result)

        self.assertEqual(len(result), 3)

    def test_do_export_file_copy(self):
        source_file = os.path.join(self.temp_dir, "source.txt")
        target_file = os.path.join(self.temp_dir, "target.txt")

        with open(source_file, "w") as f:
            f.write("test content")

        source_target_pairs = [[source_file, target_file]]
        do_export(source_target_pairs)

        self.assertTrue(os.path.exists(target_file))
        with open(target_file, "r") as f:
            self.assertEqual(f.read(), "test content")

    def test_do_export_directory_copy(self):
        source_dir = os.path.join(self.temp_dir, "source_dir")
        target_dir = os.path.join(self.temp_dir, "target_dir")
        os.makedirs(source_dir)

        source_file = os.path.join(source_dir, "nested.txt")
        with open(source_file, "w") as f:
            f.write("nested content")

        source_target_pairs = [[source_dir, target_dir]]
        do_export(source_target_pairs)

        self.assertTrue(os.path.isdir(target_dir))
        target_nested = os.path.join(target_dir, "nested.txt")
        self.assertTrue(os.path.exists(target_nested))
        with open(target_nested, "r") as f:
            self.assertEqual(f.read(), "nested content")

    def test_do_export_replace_existing_file(self):
        source_file = os.path.join(self.temp_dir, "source.txt")
        target_file = os.path.join(self.temp_dir, "target.txt")

        with open(source_file, "w") as f:
            f.write("new content")

        with open(target_file, "w") as f:
            f.write("old content")

        source_target_pairs = [[source_file, target_file]]
        do_export(source_target_pairs)

        with open(target_file, "r") as f:
            self.assertEqual(f.read(), "new content")

    def test_do_export_replace_file_with_directory(self):
        source_dir = os.path.join(self.temp_dir, "source_dir")
        target_file = os.path.join(self.temp_dir, "target.txt")
        os.makedirs(source_dir)

        with open(os.path.join(source_dir, "file.txt"), "w") as f:
            f.write("content")

        with open(target_file, "w") as f:
            f.write("old file content")

        source_target_pairs = [[source_dir, target_file]]
        do_export(source_target_pairs)

        self.assertTrue(os.path.isdir(target_file))
        self.assertTrue(os.path.exists(os.path.join(target_file, "file.txt")))

    def test_exportable_data_keys_constants(self):
        expected_keys = [
            "stats",
            "matrix",
            "triangle",
            "columns",
            "summary",
            "cluster_dir",
        ]
        self.assertEqual(EXPORTABLE_DATA_KEYS, expected_keys)

        expected_matrix_keys = ["matrix", "columns"]
        self.assertEqual(MATRIX_ONLY_EXPORTABLE_DATA_KEYS, expected_matrix_keys)


if __name__ == "__main__":
    unittest.main()
