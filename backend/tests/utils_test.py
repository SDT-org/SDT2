import unittest
import uuid
from unittest.mock import patch, MagicMock

from utils import get_child_process_info, make_doc_id, open_folder, friendly_total_time


class TestUtils(unittest.TestCase):

    @patch("utils.psutil.Process")
    def test_get_child_process_info_success(self, mock_process_class):
        mock_process = MagicMock()
        mock_child1 = MagicMock()
        mock_child1.pid = 1234
        mock_child1.cpu_percent.return_value = 25.5
        mock_child1.memory_info.return_value.rss = 1024000
        mock_child1.nice.return_value = 0

        mock_child2 = MagicMock()
        mock_child2.pid = 5678
        mock_child2.cpu_percent.return_value = 10.2
        mock_child2.memory_info.return_value.rss = 2048000
        mock_child2.nice.return_value = 5

        mock_process.children.return_value = [mock_child1, mock_child2]
        mock_process_class.return_value = mock_process

        result = get_child_process_info()

        expected = [(1234, 25.5, 1024000, 0), (5678, 10.2, 2048000, 5)]
        self.assertEqual(result, expected)
        mock_child1.cpu_percent.assert_called_once_with(interval=0.1)
        mock_child2.cpu_percent.assert_called_once_with(interval=0.1)

    @patch("utils.psutil.Process")
    def test_get_child_process_info_with_exceptions(self, mock_process_class):
        mock_process = MagicMock()
        mock_child1 = MagicMock()
        mock_child1.pid = 1234
        mock_child1.cpu_percent.return_value = 25.5
        mock_child1.memory_info.return_value.rss = 1024000
        mock_child1.nice.return_value = 0

        mock_child2 = MagicMock()
        mock_child2.cpu_percent.side_effect = Exception("Process died")

        mock_process.children.return_value = [mock_child1, mock_child2]
        mock_process_class.return_value = mock_process

        result = get_child_process_info()

        expected = [(1234, 25.5, 1024000, 0)]
        self.assertEqual(result, expected)

    @patch("utils.psutil.Process")
    def test_get_child_process_info_no_children(self, mock_process_class):
        mock_process = MagicMock()
        mock_process.children.return_value = []
        mock_process_class.return_value = mock_process

        result = get_child_process_info()

        self.assertEqual(result, [])

    @patch("utils.uuid4")
    def test_make_doc_id(self, mock_uuid4):
        mock_uuid = MagicMock()
        mock_uuid.__str__ = lambda self: "test-uuid-string"  # type: ignore
        mock_uuid4.return_value = mock_uuid

        result = make_doc_id()

        self.assertEqual(result, "test-uuid-string")
        mock_uuid4.assert_called_once()

    def test_make_doc_id_is_string(self):
        result = make_doc_id()
        self.assertIsInstance(result, str)
        uuid.UUID(result)  # will raise ValueError if not valid UUID format

    @patch("utils.subprocess.Popen")
    @patch("utils.platform.system")
    def test_open_folder_darwin(self, mock_system, mock_popen):
        mock_system.return_value = "Darwin"
        path = "/test/path"

        open_folder(path)

        mock_popen.assert_called_once_with(["open", path], close_fds=True, shell=False)

    @patch("utils.subprocess.Popen")
    @patch("utils.platform.system")
    def test_open_folder_windows(self, mock_system, mock_popen):
        mock_system.return_value = "Windows"
        path = "C:\\test\\path"

        open_folder(path)

        mock_popen.assert_called_once_with(
            ["explorer", path], close_fds=True, shell=True
        )

    @patch("utils.subprocess.Popen")
    @patch("utils.shutil.which")
    @patch("utils.platform.system")
    def test_open_folder_linux_xdg_open(self, mock_system, mock_which, mock_popen):
        mock_system.return_value = "Linux"
        mock_which.side_effect = lambda cmd: cmd if cmd == "xdg-open" else None
        path = "/test/path"

        open_folder(path)

        mock_popen.assert_called_once_with(
            ["xdg-open", path], close_fds=True, shell=False
        )

    @patch("utils.subprocess.Popen")
    @patch("utils.shutil.which")
    @patch("utils.platform.system")
    def test_open_folder_linux_gnome_open(self, mock_system, mock_which, mock_popen):
        mock_system.return_value = "Linux"
        mock_which.side_effect = lambda cmd: cmd if cmd == "gnome-open" else None
        path = "/test/path"

        open_folder(path)

        mock_popen.assert_called_once_with(
            ["gnome-open", path], close_fds=True, shell=False
        )

    @patch("utils.shutil.which")
    @patch("utils.platform.system")
    def test_open_folder_linux_no_file_manager(self, mock_system, mock_which):
        mock_system.return_value = "Linux"
        mock_which.return_value = None

        with self.assertRaises(OSError) as cm:
            open_folder("/test/path")

        self.assertIn("Could not find a file manager", str(cm.exception))

    @patch("utils.platform.system")
    def test_open_folder_unsupported_platform(self, mock_system):
        mock_system.return_value = "Unknown"

        with self.assertRaises(OSError) as cm:
            open_folder("/test/path")

        self.assertIn("Unsupported platform: unknown", str(cm.exception))

    def test_friendly_total_time_seconds_only(self):
        result = friendly_total_time(45.67)
        self.assertEqual(result, "45.67s")

    def test_friendly_total_time_with_minutes(self):
        result = friendly_total_time(125.34)  # 2m 5.34s
        self.assertEqual(result, "2m 5.34s")

    def test_friendly_total_time_exact_minutes(self):
        result = friendly_total_time(120.00)  # 2m 0.00s
        self.assertEqual(result, "2m 0.00s")

    def test_friendly_total_time_zero(self):
        result = friendly_total_time(0)
        self.assertEqual(result, "0.00s")

    def test_friendly_total_time_long_duration(self):
        result = friendly_total_time(3665.78)  # 61m 5.78s
        self.assertEqual(result, "61m 5.78s")


if __name__ == "__main__":
    unittest.main()
