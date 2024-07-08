import unittest

import pysam

from verona import ensembl, fake_vcf


class TestChunkReader(fake_vcf.TestWithTempDir):

    def test_reading_plain_vcf(self):
        """Tests reading a plain VCF file with one record."""
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
        chunks = list(ensembl.get_vcf_query_data(fake.path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 1)

    def test_reading_gzipped_vcf(self):
        """Tests reading a bgzipped-VCF file with one record."""
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
        compressed_path = fake.path.with_suffix(".vcf.gz")
        pysam.tabix_compress(str(fake.path), str(compressed_path))
        self.assertTrue(compressed_path.exists())
        chunks = list(ensembl.get_vcf_query_data(compressed_path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 1)
