import pathlib
import shutil
import tempfile
import unittest

import pysam
from varona import fake_vcf


class TestFakeVcfFiles(fake_vcf.TestWithTempDir):

    def test_minimal_empty(self):
        fake = fake_vcf.FakeVcfFile(self.path)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            self.assertEqual(vcf.header.version, "VCFv4.2")
            records = list(vcf.fetch())
            self.assertEqual(len(records), 0)

    def test_platypus_empty_no_rewrite(self):
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.write_vcf(rewrite=False)
        self.assertTrue(fake.path.exists())
        # Check for this oddity at the beginning of the header.
        first_three_lines = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileformat=VCFv4.0""".split(
            "\n"
        )
        self.assertListEqual(fake.path.read_text().split("\n")[:3], first_three_lines)
        # Ultimately, the v4.2 supercedes the v4.0 which isn't desireable.
        with pysam.VariantFile(fake.path) as vcf:
            self.assertEqual(vcf.header.version, "VCFv4.2")
            records = list(vcf.fetch())
            self.assertEqual(len(records), 0)

    def test_platypus_empty(self):
        """Tests VCF file with no records using Platypus format with rewrite."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            self.assertEqual(vcf.header.version, "VCFv4.0")
            records = list(vcf.fetch())
            self.assertEqual(len(records), 0)

    def test_with_records(self):
        records = [
            {
                "pos": 100,
                "qual": 200,
                "alleles": ("A", "G", "C"),
                "info": {"DP": 100},
                "samples": [{"GT": (1, 2)}],
            },
        ]
        fake = fake_vcf.FakeVcfFile(self.path, records=records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf.fetch())
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].pos, 100)
            self.assertEqual(records[0].qual, 200)
            self.assertEqual(records[0].alleles, ("A", "G", "C"))
            self.assertEqual(records[0].info["DP"], 100)
