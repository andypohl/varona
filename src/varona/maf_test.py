"""Tests for the maf module.
"""

import importlib
import unittest
from unittest import mock

import pysam

from varona import bcftools, fake_vcf, maf


class TestMafEnum(unittest.TestCase):

    def test_bcftools_missing(self):
        """Test the enum is reduced without bcftools."""
        more_values = set()
        fewer_values = set()
        with mock.patch("varona.bcftools.HAVE_BCFTOOLS", False):
            importlib.reload(maf)
            fewer_values = set(maf.MafMethod)
        with mock.patch("varona.bcftools.HAVE_BCFTOOLS", True):
            importlib.reload(maf)
            more_values = set(maf.MafMethod)
            self.assertSetEqual(more_values - fewer_values, {maf.MafMethod.BCFTOOLS})


class TestMafInInfo(fake_vcf.TestWithTempDir):

    def test_maf_in_info_already(self):
        """Tests MAF in INFO field."""
        fake = fake_vcf.FakeVcfFile(self.path)
        fake.header.add_meta(
            "INFO",
            items=[
                ("ID", "MAF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "Minor allele frequency"),
            ],
        )
        fake.add_records(
            [
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G", "C"),
                    "info": {"DP": 100, "MAF": 0.5},
                    "samples": [{"GT": (1, 2)}],
                },
            ]
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            self.assertAlmostEqual(maf.maf_from_bcftools(records[0]), 0.5)

    def test_maf_missing_from_info(self):
        """Tests MAF in INFO field first missing then filling it."""
        fake = fake_vcf.FakeVcfFile(self.path)
        fake.add_records(
            [
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G", "C"),
                    "info": {"DP": 100},
                    "samples": [{"GT": (1, 2)}],
                },
            ]
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            self.assertRaises(KeyError, maf.maf_from_bcftools, records[0])
        # put this in a new test soon. unittests shouldn't have if statements
        if bcftools.HAVE_BCFTOOLS:
            # fill in the MAF tag with bcftools
            with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
                records = list(vcf)
                self.assertEqual(len(records), 1)
                self.assertAlmostEqual(
                    maf.maf_from_method(records[0], maf.MafMethod.BCFTOOLS), 0.5
                )


class TestMafFromFr(fake_vcf.TestWithTempDir):

    def test_maf_from_fr(self):
        """Tests MAF from FR field."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.add_records_from_lines(
            """\
1\t1158631\t.\tA\tG\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-43.88,0.0:3:99:160:156
"""
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            self.assertAlmostEqual(maf.maf_from_method(records[0], maf.MafMethod.FR), 0)

    @unittest.skipUnless(bcftools.HAVE_BCFTOOLS, "bcftools not found")
    def test_maf_from_fr_bcftools_agreement(self):
        """Most of the time the FR->MAF and bcftools->MAF agree."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.add_records_from_lines(
            """\
1\t1158631\t.\tA\tG\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-43.88,0.0:3:99:160:156
"""
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        # put this in a new test soon. unittests shouldn't have if statements
        if bcftools.HAVE_BCFTOOLS:
            with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
                records = list(vcf)
                self.assertEqual(len(records), 1)
                self.assertAlmostEqual(
                    maf.maf_from_fr(records[0]), maf.maf_from_bcftools(records[0])
                )

    @unittest.skipUnless(bcftools.HAVE_BCFTOOLS, "bcftools not found")
    def test_maf_from_fr_bcftools_disagreement(self):
        """In a handful of cases the FR->MAF and bcftools->MAF are different."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        # Below are a select few records where the FR->MAF and bcftools->MAF
        # disagree. It's a bit strange, but in 11/12 of these instances (all
        # biallelic), the Genotype (GT field) in the sample is indicating
        # homozygosity for the alternate allele, but the FR field looks
        # like it's heterozygous.  That's how the first two records are below.
        # The third record is the only one where the FR field is not 1.0
        # or 0.5.
        fake.add_records_from_lines(
            """\
1\t160141286\t.\tC\tG\t2919\tPASS\tBRF=0.15;FR=0.5000;HP=2;HapScore=1;MGOF=33;MMLQ=37;MQ=59.44;NF=88;NR=47;PP=2919;QD=20;SC=GTTATTATCCCTGGGGTGAGA;SbPval=0.52;Source=Platypus;TC=135;TCF=88;TCR=47;TR=135;WE=160141294;WS=160141275	GT:GL:GOF:GQ:NR:NV\t1/1:-150.44,-39.49,0.0:33:99:135:135
5\t75923294\t.\tT\tG\t2962\tPASS\tBRF=0.19;FR=0.5000;HP=1;HapScore=1;MGOF=15;MMLQ=36;MQ=59.06;NF=73;NR=113;PP=2962;QD=20;SC=CCGTCGATGATGCCAACGTGG;SbPval=0.54;Source=Platypus;TC=186;TCF=73;TCR=113;TR=186;WE=75923315;WS=75923275	GT:GL:GOF:GQ:NR:NV\t1/1:-63.57,-52.38,0.0:15:99:186:186
7\t153109979\t.\tA\tG\t2962\tPASS\tBRF=0.19;FR=0.5121;HP=1;HapScore=1;MGOF=20;MMLQ=33;MQ=56.66;NF=88;NR=80;PP=2962;QD=20;SC=AGATGGCGGCAGTGAGGAAGC;SbPval=0.57;Source=Platypus;TC=173;TCF=89;TCR=84;TR=168;WE=153109987;WS=153109969\tGT:GL:GOF:GQ:NR:NV\t0/1:-300.0,0.0,-1.62:20:16:173:168
"""
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 3)
            for record in records:
                self.assertNotAlmostEqual(
                    maf.maf_from_method(record, maf.MafMethod.FR),
                    maf.maf_from_method(record, maf.MafMethod.BCFTOOLS),
                )


class TestMafFromSamples(fake_vcf.TestWithTempDir):

    def test_maf_from_one_sample_biallelic(self):
        """Tests MAF in INFO field."""
        fake = fake_vcf.FakeVcfFile(
            self.path,
            samples=["sample1"],
            records=[
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G"),
                    "info": {"DP": 100},
                    "samples": [{"GT": (0, 1)}],
                },
            ],
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            self.assertAlmostEqual(maf.maf_from_samples(records[0]), 0.5)
            # test high-level method too
            self.assertAlmostEqual(
                maf.maf_from_method(records[0], maf.MafMethod.SAMPLES), 0.5
            )

    def test_maf_from_multi_sample_biallelic(self):
        """Tests MAF in INFO field."""
        fake = fake_vcf.FakeVcfFile(
            self.path,
            samples=["sample1", "sample2", "sample3"],
            records=[
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G"),
                    "info": {"DP": 100},
                    "samples": [{"GT": (0, 1)}, {"GT": (1, 1)}, {"GT": (1, 1)}],
                },
            ],
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            calculated_maf = maf.maf_from_samples(records[0])
            # The expected MAF is 1/6 because of the 6 observed alleles in the
            # three samples, the allele frequencies are 5/6 and 1/6 and 1/6 is
            # the second highest.
            expected_maf = 1 / 6
            self.assertAlmostEqual(calculated_maf, expected_maf)
        # put this in a new test soon. unittests shouldn't have if statements
        if bcftools.HAVE_BCFTOOLS:
            # Check agreement with bcftools
            with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
                records = list(vcf)
                self.assertEqual(len(records), 1)
                self.assertAlmostEqual(
                    maf.maf_from_method(records[0], maf.MafMethod.SAMPLES),
                    maf.maf_from_method(records[0], maf.MafMethod.BCFTOOLS),
                    places=5,
                )

    def test_maf_from_multi_sample_multiallelic(self):
        """Tests MAF with five samples and three alleles."""
        fake = fake_vcf.FakeVcfFile(
            self.path,
            samples=["sample1", "sample2", "sample3", "sample4", "sample5"],
            records=[
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G", "C"),
                    "info": {"DP": 100},
                    "samples": [
                        {"GT": (0, 1)},
                        {"GT": (1, 1)},
                        {"GT": (1, 2)},
                        {"GT": (1, 1)},
                        {"GT": (1, 2)},
                    ],
                },
            ],
        )
        # Allele frequencies are:
        # A: 1/10
        # G: 7/10
        # C: 2/10
        # so the MAF (second-highest) is 2/10 = 1/5
        fake.write_vcf()
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            calculated_maf = maf.maf_from_samples(records[0])
            # The expected MAF is 1/6 because of the 6 observed alleles in the
            # three samples, the allele frequencies are 5/6 and 1/6 and 1/6 is
            # the second highest.
            expected_maf = 1 / 5
            self.assertAlmostEqual(calculated_maf, expected_maf)
        # put this in a new test soon. unittests shouldn't have if statements
        if bcftools.HAVE_BCFTOOLS:
            # Check agreement with bcftools
            with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
                records = list(vcf)
                self.assertEqual(len(records), 1)
                self.assertAlmostEqual(
                    maf.maf_from_method(records[0], maf.MafMethod.SAMPLES),
                    maf.maf_from_method(records[0], maf.MafMethod.BCFTOOLS),
                    places=5,
                )

    def test_maf_from_one_sample_missing_gt(self):
        """Tests trying to calculate MAF on samples missing GT."""
        fake = fake_vcf.FakeVcfFile(self.path)
        header = pysam.VariantHeader()
        header.add_line("##fileformat=VCFv4.2")
        header.contigs.add("1", 1000)
        header.add_meta(
            "INFO",
            items=[
                ("ID", "DP"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "Total Depth"),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "NR"),
                ("Number", "1"),
                ("Type", "Integer"),
                (
                    "Description",
                    "Number of reads covering variant location in this sample",
                ),
            ],
        )
        header.add_samples(["sample"])
        fake.header = header
        fake.add_records(
            [
                {
                    "pos": 100,
                    "qual": 200,
                    "alleles": ("A", "G"),
                    "info": {"DP": 100},
                    "samples": [{"NR": 300}],
                },
            ]
        )
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 1)
            self.assertRaises(KeyError, maf.maf_from_samples, records[0])
        # Note that bcftools only warns if the GT field is missing.
        # put this in a new test soon. unittests shouldn't have if statements
        if bcftools.HAVE_BCFTOOLS:
            # It'll then make a VCF with invalid MAF values, which isn't great.
            # Check agreement with bcftools
            with bcftools.VariantFileFilledInTags(fake.path, "r", ["MAF"]) as vcf:
                records = list(vcf)
                with self.assertRaisesRegex(KeyError, "Invalid INFO field: MAF"):
                    maf.maf_from_method(records[0], maf.MafMethod.BCFTOOLS)
