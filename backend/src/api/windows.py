import webview

from file_utils import get_html_path
import app_globals


class Windows:
    def about_window(self):
        webview.create_window(
            "About", get_html_path("about.html"), js_api=app_globals.js_api
        )

    def manual_window(self):
        webview.create_window("SDT2 Manual", get_html_path("manual.html"))
