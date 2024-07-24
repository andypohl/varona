"""Tests for the bcftools module.

These tests are skipped if bcftools isn't available on the system outside of
pysam. There is a test though that checks if the module raises an error if
bcftools isn't available.
"""

import unittest
from unittest import mock

from varona import bcftools, fake_vcf


class TestFillingInTags(fake_vcf.TestWithTempDir):

    records = [
        {
            "pos": 100,
            "qual": 200,
            "alleles": ("A", "G", "C"),
            "info": {"DP": 100},
            "samples": [{"GT": (1, 2)}],
        },
    ]

    @mock.patch("varona.bcftools.HAVE_BCFTOOLS", False)
    def test_bcftools_missing(self):
        """Test that an error is raised if bcftools isn't available."""
        fake = fake_vcf.FakeVcfFile(self.path)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with self.assertRaises(RuntimeError):
            bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"])

    @unittest.skipUnless(bcftools.HAVE_BCFTOOLS, "non-pysam bcftools is not available")
    def test_fillin_maf(self):
        """Tests the filling in of the MAF tag of the VCF."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as filledin_vcf:
            records = list(filledin_vcf.fetch())
            self.assertEqual(len(records), 1)

    @unittest.skipUnless(bcftools.HAVE_BCFTOOLS, "non-pysam bcftools is not available")
    def test_disallowed_tags(self):
        """Tests that an error is raised if a disallowed tag is used."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        bad_tag = ["BAD"]
        self.assertSetEqual(set(bad_tag) - bcftools.ALLOWED_TAGS, {"BAD"})
        with self.assertRaisesRegex(ValueError, r"Only .* tags are allowed."):
            with bcftools.VariantFileFilledInTags(
                fake.path, "r", bad_tag
            ) as filledin_vcf:
                pass
