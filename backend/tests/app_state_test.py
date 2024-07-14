import sys
import os
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from app_state import create_app_state

class TestAppState(unittest.TestCase):
    def test_state_updates(self):
        get_state, set_state = create_app_state()

        state = get_state()
        self.assertEqual(state.filename, "")
        self.assertEqual(state.progress, 0)
        self.assertEqual(state.debug, False)

        set_state(filename="test.txt", progress=25, debug=True)
        state = get_state()
        self.assertEqual(state.filename, "test.txt")
        self.assertEqual(state.progress, 25)
        self.assertEqual(state.debug, True)

if __name__ == '__main__':
    unittest.main()
