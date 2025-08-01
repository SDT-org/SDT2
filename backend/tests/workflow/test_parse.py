import sys
import os
import unittest
import tempfile
from unittest.mock import patch

from numpy import ndarray
from pandas.core.frame import DataFrame

# Path setup handled by PYTHONPATH in package.json
from workflow.runner import run_parse
from workflow.models import WorkflowResult
from mime_setup import register_mimetypes

register_mimetypes()


class TestRunParse(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.test_fastas_dir = os.path.join(os.path.dirname(__file__), "..", "fastas")

    def test_run_parse_with_valid_fasta(self):
        """Test run_parse with a valid FASTA file."""
        valid_fasta_path = os.path.join(self.test_fastas_dir, "valid.fasta")

        result = run_parse(valid_fasta_path)

        self.assertIsNone(result.error)
        self.assertIsInstance(result.seq_dict, dict)
        self.assertGreater(len(result.seq_dict), 0)
        self.assertGreater(result.max_sequence_length, 0)
        self.assertIsInstance(result.is_aa, bool)

    def test_run_parse_with_empty_fasta(self):
        """Test run_parse with an empty FASTA file."""
        empty_fasta_path = os.path.join(self.test_fastas_dir, "empty.fasta")

        result = run_parse(empty_fasta_path)

        self.assertEqual(result.error, "NOT_ENOUGH_SEQUENCES")

    def test_run_parse_with_single_sequence(self):
        """Test run_parse with a FASTA file containing only one sequence."""
        one_fasta_path = os.path.join(self.test_fastas_dir, "one.fasta")

        result = run_parse(one_fasta_path)

        self.assertEqual(result.error, "NOT_ENOUGH_SEQUENCES")

    def test_run_parse_with_duplicate_names(self):
        """Test run_parse with a FASTA file containing duplicate sequence names."""
        duplicate_fasta_path = os.path.join(
            self.test_fastas_dir, "duplicate_name.fasta"
        )

        result = run_parse(duplicate_fasta_path)

        self.assertEqual(result.error, "DUPLICATE_SEQUENCE_NAME")

    def test_run_parse_with_zero_length_sequence(self):
        """Test run_parse with a FASTA file containing zero-length sequences."""
        zero_length_fasta_path = os.path.join(self.test_fastas_dir, "zero_length.fasta")

        result = run_parse(zero_length_fasta_path)

        self.assertEqual(result.error, "ZERO_LENGTH_SEQUENCE")

    def test_run_parse_with_gap_only_sequence(self):
        """Test run_parse with a FASTA file containing gap-only sequences."""
        gap_only_fasta_path = os.path.join(self.test_fastas_dir, "gap_only.fasta")

        result = run_parse(gap_only_fasta_path)

        self.assertEqual(result.error, "ZERO_LENGTH_SEQUENCE")

    def test_run_parse_with_spaces_in_sequence_names(self):
        """Test run_parse with a FASTA file containing spaces in sequence names."""
        spaces_fasta_path = os.path.join(self.test_fastas_dir, "spaces.fasta")

        result = run_parse(spaces_fasta_path)

        # Should succeed but format the sequence names
        self.assertIsNone(result.error)
        # Check that spaces are replaced with underscores
        for seq_id in result.seq_dict.keys():
            self.assertNotIn(" ", seq_id)

    def test_run_parse_with_too_long_sequence_names(self):
        """Test run_parse with a FASTA file containing very long sequence names."""
        too_long_fasta_path = os.path.join(self.test_fastas_dir, "too_long.fasta")

        result = run_parse(too_long_fasta_path)

        # Should succeed but truncate long names
        self.assertIsNone(result.error)
        for seq_id in result.seq_dict.keys():
            self.assertLessEqual(len(seq_id), 50)

    def test_run_parse_with_invalid_file_type(self):
        """Test run_parse with a non-FASTA file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as tmp_file:
            tmp_file.write("This is not a FASTA file")
            tmp_file_path = tmp_file.name

        try:
            result = run_parse(tmp_file_path)
            self.assertEqual(result.error, "INVALID_FILE_TYPE")
        finally:
            os.unlink(tmp_file_path)

    def test_run_parse_with_nonexistent_file(self):
        """Test run_parse with a file that doesn't exist."""
        nonexistent_path = "/path/that/does/not/exist.fasta"

        result = run_parse(nonexistent_path)

        assert result.error is not None
        self.assertTrue(result.error.startswith("INVALID_FASTA_FORMAT"))

    def test_run_parse_with_malformed_fasta(self):
        """Test run_parse with a malformed FASTA file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write("This is malformed FASTA content without proper headers")
            tmp_file_path = tmp_file.name

        try:
            result = run_parse(tmp_file_path)
            # Should either error or handle gracefully depending on Bio.SeqIO behavior
            if result.error:
                self.assertTrue(
                    result.error.startswith("INVALID_FASTA_FORMAT")
                    or result.error == "NOT_ENOUGH_SEQUENCES"
                )
        finally:
            os.unlink(tmp_file_path)

    def test_run_parse_result_structure(self):
        """Test that run_parse returns a properly structured WorkflowResult."""
        valid_fasta_path = os.path.join(self.test_fastas_dir, "valid.fasta")

        result = run_parse(valid_fasta_path)

        self.assertIsInstance(result, WorkflowResult)

        self.assertTrue(hasattr(result, "seq_dict"))
        self.assertTrue(hasattr(result, "ordered_ids"))
        self.assertTrue(hasattr(result, "reordered_ids"))
        self.assertTrue(hasattr(result, "max_sequence_length"))
        self.assertTrue(hasattr(result, "warnings"))
        self.assertTrue(hasattr(result, "error"))
        self.assertTrue(hasattr(result, "distance_matrix"))
        self.assertTrue(hasattr(result, "similarity_matrix"))
        self.assertTrue(hasattr(result, "is_aa"))
        self.assertTrue(hasattr(result, "min_score"))


class TestRunParseInitialState(unittest.TestCase):
    """Test the initial state setup in run_parse function."""

    def setUp(self):
        self.test_fastas_dir = os.path.join(os.path.dirname(__file__), "..", "fastas")

    def test_initial_workflow_result_creation(self):
        """Test that run_parse creates a proper initial WorkflowResult."""
        valid_fasta_path = os.path.join(self.test_fastas_dir, "valid.fasta")

        with patch("workflow.parse.run") as mock_parse_run:
            # Mock parse.run to return the input result unchanged
            mock_parse_run.side_effect = lambda result, _: result

            result = run_parse(valid_fasta_path)

            self.assertEqual(result.seq_dict, {})
            self.assertEqual(result.ordered_ids, [])
            self.assertEqual(result.reordered_ids, [])
            self.assertEqual(result.max_sequence_length, 0)
            self.assertEqual(result.warnings, [])
            self.assertIsNone(result.error)
            self.assertIsNone(result.is_aa)
            self.assertEqual(result.min_score, 0)

    def test_parse_module_called_correctly(self):
        """Test that run_parse calls the parse module correctly."""
        valid_fasta_path = os.path.join(self.test_fastas_dir, "valid.fasta")

        with patch("workflow.parse.run") as mock_parse_run:
            mock_result = WorkflowResult(
                seq_dict={"seq1": "ATCG", "seq2": "GCTA"},
                ordered_ids=[],
                reordered_ids=[],
                max_sequence_length=4,
                warnings=[],
                error=None,
                distance_matrix=ndarray([]),
                similarity_matrix=DataFrame(),
                is_aa=False,
                min_score=0,
            )
            mock_parse_run.return_value = mock_result

            run_parse(valid_fasta_path)

            mock_parse_run.assert_called_once()
            args, _ = mock_parse_run.call_args
            self.assertEqual(len(args), 2)
            self.assertIsInstance(args[0], WorkflowResult)
            self.assertEqual(args[1], valid_fasta_path)

    def test_error_handling_in_parse(self):
        """Test that run_parse properly handles errors from parse module."""
        valid_fasta_path = os.path.join(self.test_fastas_dir, "valid.fasta")

        with patch("workflow.parse.run") as mock_parse_run:
            error_result = WorkflowResult(
                seq_dict={},
                ordered_ids=[],
                reordered_ids=[],
                max_sequence_length=0,
                warnings=[],
                error="TEST_ERROR",
                distance_matrix=ndarray([]),
                similarity_matrix=DataFrame(),
                is_aa=None,
                min_score=0,
            )
            mock_parse_run.return_value = error_result

            result = run_parse(valid_fasta_path)

            self.assertEqual(result.error, "TEST_ERROR")
            self.assertEqual(result, error_result)
    
    def test_missing_sequence_id_error(self):
        """Test MISSING_SEQUENCE_ID error case."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_file:
            tmp_file.write(">seq1\nATCG\n>\nGCTA\n")
            tmp_file_path = tmp_file.name
        
        try:
            result = run_parse(tmp_file_path)
            self.assertEqual(result.error, "MISSING_SEQUENCE_ID")
        finally:
            os.unlink(tmp_file_path)
    
    def test_amino_acid_detection_on_third_sequence(self):
        """Test amino acid detection logic that runs every 3rd sequence."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_file:
            tmp_file.write(">seq1\nATCG\n>seq2\nGCTA\n>seq3\nEFIL\n")
            tmp_file_path = tmp_file.name
        
        try:
            result = run_parse(tmp_file_path)
            self.assertIsNone(result.error)
            self.assertTrue(result.is_aa)
        finally:
            os.unlink(tmp_file_path)
    
    def test_helper_functions(self):
        """Test unused helper functions for complete coverage."""
        from workflow.parse import format_sequence_record, residue_check
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        
        # Test format_sequence_record
        record = SeqRecord(Seq("atcg-123"), id="test seq", description="desc")
        formatted = format_sequence_record(record)
        self.assertEqual(formatted.id, "test_seq")
        self.assertEqual(str(formatted.seq), "ATCG")
        
        # Test residue_check
        self.assertTrue(residue_check("EFIL"))
        self.assertFalse(residue_check("ATCG"))


if __name__ == "__main__":
    unittest.main()
