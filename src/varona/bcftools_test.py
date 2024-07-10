"""Tests for the bcftools module.

These tests are skipped if bcftools isn't available on the system outside of
pysam. There is a test though that checks if the module raises an error if
bcftools isn't available.
"""

import unittest
from unittest import mock

from varona import bcftools, fake_vcf


class TestFillingInTags(fake_vcf.TestWithTempDir):

    @mock.patch("varona.bcftools.HAVE_BCFTOOLS", False)
    def test_bcftools_missing(self):
        """Test that an error is raised if bcftools isn't available."""
        fake = fake_vcf.FakeVcfFile(self.path)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with self.assertRaises(RuntimeError):
            bcftools.VariantFileFilledInTags(fake.path, ["MAF"])

    @unittest.skipUnless(bcftools.HAVE_BCFTOOLS, "non-pysam bcftools is not available")
    def test_fillin_maf(self):
        """Tests the filling in of the MAF tag of the VCF."""
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
        with bcftools.VariantFileFilledInTags(fake.path, ["MAF"]) as filledin_vcf:
            records = list(filledin_vcf.fetch())
            self.assertEqual(len(records), 1)
