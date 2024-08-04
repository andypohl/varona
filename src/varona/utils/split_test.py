"""Test the VCF split stuff.
"""

import pathlib
import unittest

import pysam

from varona import fake_vcf
from varona.utils import split


class TestSplit(fake_vcf.TestWithTempDir):
    """Test splitting a VCF file."""

    records = [
        {
            "pos": 100,
            "alleles": ("A", "G", "C"),
            "info": {"DP": 100},
            "samples": [{"GT": (1, 2)}],
        },
        {
            "pos": 120,
            "alleles": ("G", "C"),
            "info": {"DP": 222},
            "samples": [{"GT": (1, 1)}],
        },
        {
            "pos": 130,
            "alleles": ("G", "1"),
            "info": {"DP": 381},
            "samples": [{"GT": (0, 1)}],
        },
        {
            "pos": 140,
            "alleles": ("A", "C"),
            "info": {"DP": 88},
            "samples": [{"GT": (1, 1)}],
        },
        {
            "pos": 150,
            "alleles": ("A", "G"),
            "info": {"DP": 97},
            "samples": [{"GT": (1, 1)}],
        },
    ]

    def _test_split_vcf(
        self, split_func, split_arg, expected_n_records, expected_files
    ):
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        out_dir = self.tmp_dir / "out"
        split_func(fake.path, out_dir, **split_arg)
        self.assertTrue(out_dir.exists())
        out_paths = [
            pathlib.Path(y) for y in sorted([str(x) for x in out_dir.iterdir()])
        ]
        self.assertEqual(len(out_paths), len(expected_files))
        for file_ix, out_path in enumerate(out_paths, start=1):
            self.assertTrue(out_path.exists())
            self.assertTrue(out_path.is_file())
            self.assertEqual(out_path.name, expected_files[file_ix - 1])
            with pysam.VariantFile(out_path, "r") as f:
                self.assertEqual(f.compression, "BGZF")
                self.assertTrue(f.is_vcf)
                n_records = sum(1 for _ in f)
                self.assertEqual(n_records, expected_n_records[file_ix - 1])

    def test_split_vcf_three_pieces(self):
        """Tests splitting a VCF file with n_chunks arg."""
        self._test_split_vcf(
            split.split_vcf,
            {"n_chunks": 3},
            [2, 2, 1],
            [f"test_0{file_ix}.vcf.gz" for file_ix in range(1, 4)],
        )

    def test_split_vcf_chunk_3(self):
        """Tests splitting a VCF file chunk arg."""
        self._test_split_vcf(
            split.split_vcf,
            {"chunk_size": 3},
            [3, 2],
            [f"test_0{file_ix}.vcf.gz" for file_ix in range(1, 3)],
        )

    def test_no_input(self):
        """Tests raising an error when no input is given."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        with self.assertRaises(ValueError):
            split.split_vcf(None, None)
        with self.assertRaisesRegex(ValueError, "Either chunk_size or n_chunks"):
            split.split_vcf(fake.path, self.tmp_dir)


class TestVcfSplitArgs(unittest.TestCase):
    """Test the argument parser for the vcf split CLI."""

    def setUp(self):
        self.parser = split.vcf_split_args()

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
            ["input.vcf", "--chunk-size", "100", "--out-dir", "."]
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
