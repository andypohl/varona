"""Test the VCF split stuff.
"""

import pathlib
import unittest

from varona.utils import cli_vcf_split


class TestVcfSplitArgs(unittest.TestCase):
    """Test the argument parser for the vcf split CLI."""

    def setUp(self):
        self.parser = cli_vcf_split.vcf_split_args()

    def test_in_vcf_argument(self):
        args = self.parser.parse_args(["input.vcf", "--chunk-size", "100"])
        self.assertEqual(args.in_vcf, pathlib.Path("input.vcf"))

    def test_default_log_level(self):
        args = self.parser.parse_args(["input.vcf", "--chunk-size", "100"])
        self.assertEqual(args.log_level, "warning")

    def test_custom_log_level(self):
        args = self.parser.parse_args(
            ["input.vcf", "--log-level", "debug", "--chunk-size", "100"]
        )
        self.assertEqual(args.log_level, "debug")

    def test_default_out_dir(self):
        args = self.parser.parse_args(
            ["input.vcf", "--chunk-size", "100", "--out-dir", str(pathlib.Path.cwd())]
        )
        self.assertEqual(args.out_dir, pathlib.Path.cwd())

    def test_custom_out_dir(self):
        custom_dir = pathlib.Path("/custom/dir")
        args = self.parser.parse_args(
            ["input.vcf", "--n-chunks", "5", "--out-dir", str(custom_dir)]
        )
        self.assertEqual(args.out_dir, custom_dir)

    def test_chunk_size_argument(self):
        args = self.parser.parse_args(["input.vcf", "--chunk-size", "100"])
        self.assertEqual(args.chunk_size, 100)
        self.assertIsNone(args.n_chunks)

    def test_n_chunks_argument(self):
        args = self.parser.parse_args(["input.vcf", "--n-chunks", "5"])
        self.assertEqual(args.n_chunks, 5)
        self.assertIsNone(args.chunk_size)

    def test_mutually_exclusive_arguments(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args(
                ["input.vcf", "--chunk-size", "100", "--n-chunks", "5"]
            )

    def test_missing_required_argument(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args(["input.vcf"])
