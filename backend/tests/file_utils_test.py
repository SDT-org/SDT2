import unittest
import json
import tempfile
import os
from unittest.mock import patch, mock_open

from file_utils import read_json_file


class TestFileUtils(unittest.TestCase):

    def test_read_json_file_valid_json(self):
        test_data = {"key": "value", "number": 42}

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            json.dump(test_data, f)
            temp_path = f.name

        try:
            result = read_json_file(temp_path)
            self.assertEqual(result, test_data)
        finally:
            os.unlink(temp_path)

    def test_read_json_file_empty_file(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            temp_path = f.name

        try:
            result = read_json_file(temp_path)
            self.assertIsNone(result)
        finally:
            os.unlink(temp_path)

    def test_read_json_file_whitespace_only(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            f.write("   \n\t  ")
            temp_path = f.name

        try:
            result = read_json_file(temp_path)
            self.assertIsNone(result)
        finally:
            os.unlink(temp_path)

    def test_read_json_file_invalid_json(self):
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            f.write("{ invalid json")
            temp_path = f.name

        try:
            with self.assertRaises(json.JSONDecodeError):
                read_json_file(temp_path)
        finally:
            os.unlink(temp_path)

    def test_read_json_file_nonexistent(self):
        with self.assertRaises(FileNotFoundError):
            read_json_file("/nonexistent/path.json")

    @patch("builtins.open", mock_open(read_data='{"test": "data"}'))
    def test_read_json_file_with_mock(self):
        result = read_json_file("mock_file.json")
        self.assertEqual(result, {"test": "data"})


if __name__ == "__main__":
    unittest.main()
