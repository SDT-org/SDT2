import sys
import os
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from app_state import create_app_state


class TestAppState(unittest.TestCase):
    def test_state_updates(self):
        get_state, set_state, _, _, _, _, _, _ = create_app_state()

        self.assertEqual(get_state().debug, False)
        set_state(debug=True)
        self.assertEqual(get_state().debug, True)

    def test_document_updates(self):
        get_state, _, _, new_document, _, update_document, _, _ = create_app_state()

        new_document("test")

        state = get_state().documents[0]
        self.assertEqual(state.filename, "")
        self.assertEqual(state.progress, 0)

        update_document(id="test", filename="test.txt", progress=25)
        state = get_state().documents[0]
        self.assertEqual(state.filename, "test.txt")
        self.assertEqual(state.progress, 25)


if __name__ == "__main__":
    unittest.main()
