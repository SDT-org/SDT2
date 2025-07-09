import os
import sys
import webview

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from config.config import app_version, dev_frontend_host

is_compiled = "__compiled__" in globals()


def file_exists(path):
    return os.path.exists(os.path.join(os.path.dirname(__file__), "..", path))


def get_html_path(filename="index.html"):
    if is_compiled:
        if file_exists(f"./gui/{filename}"):
            return f"./gui/{filename}"
        raise Exception(f"{filename} not found")
    else:
        return f"{dev_frontend_host}/{filename}"


def about_window():
    webview.create_window("About", get_html_path("about.html"))


def manual_window():
    webview.create_window("SDT2 Manual", get_html_path("manual.html"))
