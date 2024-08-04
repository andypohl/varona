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
        self,
        split_func,
        split_arg,
        compress_arg,
        expected_n_records,
        expected_files,
        expected_compression,
    ):
        """Helper method to test splitting a VCF file.

        :param split_func: The function to split the VCF file.
        :param split_arg: The arguments to pass to the split function.
        :param compress_arg: The arguments to pass to the compress function.
        :param expected_n_records: The expected number of records in each split.
        :param expected_files: The expected names of the split files.
        :param expected_compression: The expected compression type.
        """
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        out_dir = self.tmp_dir / "out"
        split_func(fake.path, out_dir, **split_arg, **compress_arg)
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
                self.assertEqual(f.compression, expected_compression)
                self.assertTrue(f.is_vcf)
                n_records = sum(1 for _ in f)
                self.assertEqual(n_records, expected_n_records[file_ix - 1])

    def test_split_vcf_three_pieces(self):
        """Tests splitting a VCF file with n_chunks arg."""
        self._test_split_vcf(
            split.split_vcf,
            {"n_chunks": 3},
            {"compress": True},
            [2, 2, 1],
            [f"test_0{file_ix}.vcf.gz" for file_ix in range(1, 4)],
            "BGZF",
        )

    def test_split_vcf_chunk_3(self):
        """Tests splitting a VCF file chunk arg."""
        self._test_split_vcf(
            split.split_vcf,
            {"chunk_size": 3},
            {"compress": True},
            [3, 2],
            [f"test_0{file_ix}.vcf.gz" for file_ix in range(1, 3)],
            "BGZF",
        )

    def test_split_vcf_chunk_3_no_compress(self):
        """Tests splitting a VCF file chunk/no compress args."""
        self._test_split_vcf(
            split.split_vcf,
            {"chunk_size": 3},
            {"compress": False},
            [3, 2],
            [f"test_0{file_ix}.vcf" for file_ix in range(1, 3)],
            "NONE",
        )

    def test_no_input(self):
        """Tests raising an error when no input is given."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        with self.assertRaises(ValueError):
            split.split_vcf(None, None)
        with self.assertRaisesRegex(ValueError, "Either chunk_size or n_chunks"):
            split.split_vcf(fake.path, self.tmp_dir)
