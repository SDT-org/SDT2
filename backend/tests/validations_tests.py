import sys
import os
import io
import unittest
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from app import validate_fasta

def fasta_path(filename):
    return os.path.join(os.path.dirname(__file__), 'fastas', filename)

class TestValidateFasta(unittest.TestCase):
    def test_valid_fasta(self):
        result, message = validate_fasta(fasta_path('valid.fasta'))
        self.assertTrue(result)
        self.assertEqual(message, "VALID")

    def test_sequence_too_long(self):
        result, message = validate_fasta(fasta_path('sequences_too_long.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "SEQUENCE_TOO_LONG")

    def test_empty_fasta(self):
        result, message = validate_fasta(fasta_path('empty.fasta'))
        self.assertTrue(result)
        self.assertEqual(message, "VALID")

    def test_spaces(self):
        result, message = validate_fasta(fasta_path('spaces.fasta'))
        self.assertTrue(result)
        self.assertEqual(message, "VALID")

if __name__ == '__main__':
    unittest.main()
