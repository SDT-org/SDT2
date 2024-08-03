import sys
import os
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

    def test_sequence_too_long_fasta(self):
        result, message = validate_fasta(fasta_path('too_long.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "SEQUENCE_TOO_LONG")

    def test_empty_fasta(self):
        result, message = validate_fasta(fasta_path('empty.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "NOT_ENOUGH_SEQUENCES")

    def test_zero_length_fasta(self):
        result, message = validate_fasta(fasta_path('zero_length.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "ZERO_LENGTH_SEQUENCE")

    def test_one_sequence_fasta(self):
        result, message = validate_fasta(fasta_path('one.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "NOT_ENOUGH_SEQUENCES")

    def test_duplicate_name_fasta(self):
        result, message = validate_fasta(fasta_path('duplicate_name.fasta'))
        self.assertFalse(result)
        self.assertEqual(message, "DUPLICATE_SEQUENCE_NAME")

    def test_spaces_fasta(self):
        result, message = validate_fasta(fasta_path('spaces.fasta'))
        self.assertTrue(result)
        self.assertEqual(message, "VALID")

if __name__ == '__main__':
    unittest.main()
