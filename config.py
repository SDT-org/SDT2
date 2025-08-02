from multiprocessing import cpu_count as get_cpu_count
import platform
from tempfile import TemporaryDirectory


app_version = "2.0.0-final"
dev_frontend_host = "http://localhost:5173"
is_compiled = "__compiled__" in globals()
is_macos = platform.system() == "Darwin"
is_windows = platform.system() == "Windows"
window_title_suffix = "" if is_macos else " - SDT2"
temp_dir = TemporaryDirectory()

try:
    cpu_count = get_cpu_count()
except:
    cpu_count = 1

matrix_filetypes = ("text/csv", "application/vnd.ms-excel", "text/plain")
default_window_title = "Sequence Demarcation Tool 2"

reorder_methods = [
    "single",
    "complete",
    "average",
    "weighted",
    "centroid",
    "median",
    "ward",
]
