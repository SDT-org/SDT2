import sys
import os
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from app_state import initialize_app_state, get_state, set_state, new_document, update_document


class TestAppState(unittest.TestCase):
    def setUp(self):
        # Initialize app state for each test
        initialize_app_state()

    def test_state_updates(self):
        self.assertEqual(get_state().debug, False)
        set_state(debug=True)
        self.assertEqual(get_state().debug, True)

    def test_document_updates(self):
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
