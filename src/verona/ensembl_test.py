import os
import unittest
from unittest import mock

import httpx
import pysam

from verona import ensembl, extract, fake_vcf

VERONA_LIVE_TESTS = os.getenv("VERONA_LIVE_TESTS", "0") == "1"
"""Set to 1 to enable live querying tests for local development.
"""


class TestQueryVepApi(unittest.TestCase):

    # A chunk of two records that'll be ignored anyway.
    chunk = [
        "1 1158631 . A G . . .",
        "1 91859795 . TATGTGA CATGTGA,CATGTGG . . .",
    ]
    # Mock responses
    response_429 = httpx.Response(
        status_code=429,
        headers={"Retry-After": "1"},
        json={"error": "Too Many Requests"},
    )
    response_429_no_retry = httpx.Response(
        status_code=429,
        json={"error": "Too Many Requests"},
    )
    response_200 = httpx.Response(status_code=200, json={"data": "success"})

    @mock.patch("httpx.post")
    @mock.patch("time.sleep", return_value=None)  # disable sleeping
    def test_query_and_retries(self, mock_sleep, mock_post):
        mock_post.side_effect = [self.response_429, self.response_200]
        # Call the function
        with self.assertLogs("verona.ensembl", level="WARNING") as cm:
            result = ensembl.query_vep_api(chunk=self.chunk)
            self.assertEqual(len(cm.output), 1)  # Verify warning message
            self.assertEqual(
                cm.output[0],
                "WARNING:verona.ensembl:API code 429 with Retry-After. Retrying after 1 seconds.",
            )
        self.assertEqual(result, {"data": "success"})  # Verify the result
        self.assertEqual(mock_post.call_count, 2)  # Verify two API calls (retry)
        mock_sleep.assert_called_once_with(1)  # Verify retry delay

    @mock.patch("httpx.post")
    @mock.patch("time.sleep", return_value=None)  # disable sleeping
    def test_too_many_retries(self, mock_sleep, mock_post):
        mock_post.side_effect = [
            self.response_429,  # mock response 1
            self.response_429,  # mock response 2
            self.response_429,  # mock response 3
            self.response_200,  # (doesn't get this far)
        ]
        # Call the function
        with self.assertLogs("verona.ensembl", level="WARNING") as cm:
            with self.assertRaises(TimeoutError):
                ensembl.query_vep_api(chunk=self.chunk, retries=2)
            self.assertEqual(len(cm.output), 3)  # Verify warning message
            for i in range(2):
                self.assertEqual(
                    cm.output[i],
                    "WARNING:verona.ensembl:API code 429 with Retry-After. Retrying after 1 seconds.",
                )
        self.assertEqual(mock_post.call_count, 3)


class TestChunkReader(fake_vcf.TestWithTempDir):

    # A few records to test with.
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

    def test_reading_plain_vcf(self):
        """Tests reading a plain VCF file with one record."""

        fake = fake_vcf.FakeVcfFile(self.path, records=[self.records[0]])
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        chunks = list(ensembl.vcf_to_vep_query_data(fake.path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 1)

    def test_reading_gzipped_vcf(self):
        """Tests reading a bgzipped-VCF file with one record."""
        fake = fake_vcf.FakeVcfFile(self.path, records=[self.records[0]])
        fake.write_vcf()
        compressed_path = fake.path.with_suffix(".vcf.gz")
        pysam.tabix_compress(str(fake.path), str(compressed_path))
        self.assertTrue(compressed_path.exists())
        chunks = list(ensembl.vcf_to_vep_query_data(compressed_path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 1)

    def test_chunking(self):
        """Tests the chunking of data collected from a VCF file."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        chunks = list(ensembl.vcf_to_vep_query_data(fake.path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 5)
        # 3 chunks if the chunk size is 2, the last chunk has 1 record
        chunks = list(ensembl.vcf_to_vep_query_data(fake.path, chunk_size=2))
        self.assertEqual(len(chunks), 3)
        self.assertEqual(len(chunks[0]), 2)
        self.assertEqual(len(chunks[1]), 2)
        self.assertEqual(len(chunks[2]), 1)

    def test_chunk_data(self):
        """Tests the list of strings returned by `func:ensembl.vcf_to_vep_query_data`

        Make sure the data is in the correct format for the Ensembl API.
        """
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())
        chunks = list(ensembl.vcf_to_vep_query_data(fake.path))
        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 5)
        expected = [
            "1 100 . A G,C . . .",
            "1 120 . G C . . .",
            "1 130 . G 1 . . .",
            "1 140 . A C . . .",
            "1 150 . A G . . .",
        ]
        self.assertListEqual(chunks[0], expected)


@unittest.skipUnless(VERONA_LIVE_TESTS, "Live querying tests are disabled by default.")
class TestLiveQuerying(unittest.TestCase):
    """Tests the live Ensembl API.

    Normally these should be disabled. They are here for development purposes and
    require opting-in by setting the environment variable `VERONA_LIVE_TESTS=1`.
    """

    def test_two_records(self):
        """Tests querying the Ensembl API with two records."""
        chunk = [
            "1 1158631 . A G . . .",  # biallelic
            "1 91859795 . TATGTGA CATGTGA,CATGTGG . . .",  # > 2 alleles
        ]
        data = ensembl.query_vep_api(
            chunk, response_extractor=extract.default_vep_response_extractor
        )
        self.assertEqual(len(data), 2)
        expected = [
            {
                "contig": "1",
                "pos": 1158631,
                "ref": "A",
                "alt": "G",
                "type": "SNV",
                "effect": "synonymous_variant",
                "gene_name": "SDF4",
                "gene_id": "ENSG00000078808",
                "transcript_id": "ENST00000360001",
            },
            {
                "contig": "1",
                "pos": 91859795,
                "ref": "TATGTGA",
                "alt": "CATGTGA,CATGTGG",
                "type": "substitution",
                "effect": "missense_variant",
                "gene_name": "HFM1",
                "gene_id": "ENSG00000162669",
                "transcript_id": "ENST00000370425",
            },
        ]
        self.assertDictEqual(data[0], expected[0])
        self.assertDictEqual(data[1], expected[1])
