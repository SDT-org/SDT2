import unittest
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from document_paths import (
    build_document_paths,
    DocumentPaths,
    ImagePaths,
    ImageFormatPaths,
    DATA_FILES,
    IMAGE_KEYS,
    IMAGE_FORMATS,
)


class TestDocumentPathBuilder(unittest.TestCase):
    def setUp(self):
        self.test_doc_path = "/tmp/test_doc"
        self.result = build_document_paths(self.test_doc_path)

    def test_return_type(self):
        self.assertIsInstance(self.result, DocumentPaths)

    def test_data_file_paths(self):
        for key, filename in DATA_FILES.items():
            expected_path = os.path.join(self.test_doc_path, filename)
            self.assertEqual(getattr(self.result, key), expected_path)

    def test_images_attribute(self):
        self.assertIsInstance(self.result.images, ImagePaths)

    def test_image_paths(self):
        for key in IMAGE_KEYS:
            img_paths = getattr(self.result.images, key)
            self.assertIsInstance(img_paths, ImageFormatPaths)

            for format in IMAGE_FORMATS:
                expected_path = os.path.join(self.test_doc_path, f"{key}.{format}")
                self.assertEqual(getattr(img_paths, format), expected_path)

    def test_with_relative_path(self):
        path = "relative/path"
        result = build_document_paths(path)

        for key, filename in DATA_FILES.items():
            expected_path = os.path.join(path, filename)
            self.assertEqual(getattr(result, key), expected_path)

    def test_with_empty_path(self):
        empty_result = build_document_paths("")

        for key, filename in DATA_FILES.items():
            # with an empty path, the result should just be the filename
            self.assertEqual(getattr(empty_result, key), filename)

    def test_path_with_trailing_slash(self):
        path_with_slash = "/tmp/test_doc/"
        slash_result = build_document_paths(path_with_slash)

        for key, filename in DATA_FILES.items():
            expected_path = os.path.join(path_with_slash, filename)
            self.assertEqual(getattr(slash_result, key), expected_path)


if __name__ == "__main__":
    unittest.main()
