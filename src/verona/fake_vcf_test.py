import pathlib
import unittest
import tempfile
import shutil
import pysam

from verona import fake_vcf


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
        with pysam.VariantFile(fake.path) as vcf:
            self.assertEqual(vcf.header.version, "VCFv4.2")
            records = list(vcf.fetch())
            self.assertEqual(len(records), 0)

    def test_platypus_empty(self):
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
