"""Test module for Varona high-level functions.
"""

import unittest
from unittest import mock

import httpx
import polars as pl
import polars.testing as plt
import pysam

from varona import dataframe, ensembl, fake_vcf


class TestVcfDataFrame(fake_vcf.TestWithTempDir):
    """Test the :func:`varona.vcf_data_frame` function."""

    records = [
        {
            "pos": 100,
            "qual": 200,
            "alleles": ("A", "G", "C"),
            "info": {"DP": 100},
            "samples": [{"GT": (1, 2)}],
        },
        {
            "pos": 150,
            "qual": 200,
            "alleles": ("A", "G"),
            "info": {"DP": 90},
            "samples": [{"GT": (1, 1)}],
        },
    ]

    @staticmethod
    def vcf_extractor_a(record: pysam.VariantRecord):
        return {
            "pos": record.pos,
            "qual": record.qual,
        }

    @staticmethod
    def vcf_extractor_b(record: pysam.VariantRecord):
        return {
            "pos": record.pos,
            "alleles": "/".join(record.alleles),
        }

    schema_a = pl.Schema(
        {
            "pos": pl.UInt32,
            "qual": pl.Float32,
        }
    )
    schema_b = pl.Schema(
        {
            "pos": pl.UInt32,
            "alleles": pl.Utf8,
        }
    )

    def test_vcf_dataframe(self):
        """Tests the high-level vcf_dataframe() function."""
        fake = fake_vcf.FakeVcfFile(self.path, records=self.records)
        fake.write_vcf()
        self.assertTrue(fake.path.exists())

        # two different extractors makes two different dataframes

        df = dataframe.vcf_dataframe(
            fake.path, self.vcf_extractor_a, schema=self.schema_a
        )
        expected_df = pl.DataFrame(
            {
                "pos": [100, 150],
                "qual": [200, 200],
            },
            schema=self.schema_a,
        )
        plt.assert_frame_equal(df, expected_df)
        df = dataframe.vcf_dataframe(
            fake.path, self.vcf_extractor_b, schema=self.schema_b
        )
        expected_df = pl.DataFrame(
            {
                "pos": [100, 150],
                "alleles": ["A/G/C", "A/G"],
            },
            schema=self.schema_b,
        )
        plt.assert_frame_equal(df, expected_df)

    def test_vcf_dataframe_empty(self):
        """Try on no records."""
        fake = fake_vcf.FakeVcfFile(self.path, records=[])
        fake.write_vcf()
        self.assertTrue(fake.path.exists())

        def test_extractor_a(record: pysam.VariantRecord):
            return {
                "pos": record.pos,
                "qual": record.qual,
            }

        schema_a = pl.Schema(
            {
                "pos": pl.UInt32,
                "qual": pl.Float32,
            }
        )
        df = dataframe.vcf_dataframe(fake.path, test_extractor_a, schema=schema_a)
        expected_df = pl.DataFrame(
            {},
            schema=schema_a,
        )
        plt.assert_frame_equal(df, expected_df)


class TestVepApiDataFrame(unittest.TestCase):
    """Test the varona.vep_api_dataframe() function.

    It'll mock the HTTP client's post() method response.
    """

    def setUp(self):
        self.client = httpx.Client()

    def tearDown(self):
        self.client.close()

    response = [
        {
            "seq_region_name": "1",
            "start": 100,
            "allele_string": "A/G/C",
            "transcript_consequences": [
                {
                    "gene_symbol": "Godzilla",
                },
            ],
        },
        {
            "seq_region_name": "1",
            "start": 150,
            "allele_string": "A/G",
            "transcript_consequences": [
                {
                    "gene_symbol": "Mothra",
                },
            ],
        },
    ]

    @staticmethod
    def api_extractor_a(response_item: dict):
        return {
            "pos": response_item["start"],
            "alleles": response_item["allele_string"],
        }

    @staticmethod
    def api_extractor_b(response_item: dict):
        return {
            "pos": response_item["start"],
            "gene_name": response_item["transcript_consequences"][0]["gene_symbol"],
        }

    def check_vep_api_dataframe(self, mock_post, extractor, response, expected_df):
        """Check the varona.vep_api_dataframe() function."""
        mock_200 = mock.MagicMock()
        mock_200.status_code = 200
        mock_200.json.return_value = response
        # the test here mocks the response so the list of loci isn't consequential
        # but the input is real in any case
        loci_list = [
            "1 100 . A G,C . . .",
            "1 150 . G C . . .",
        ]
        mock_post.return_value = mock_200
        df = dataframe.vep_api_dataframe(
            self.client,
            loci_list,
            ensembl.Assembly.GRCH37,
            extractor,
            schema=expected_df.schema,
        )
        plt.assert_frame_equal(df, expected_df)

    @mock.patch("httpx.Client.post")
    def test_vep_api_dataframe(self, mock_post):
        """Test the high-level varona.vep_api_dataframe() function.

        Test a couple different extractors.
        """

        self.check_vep_api_dataframe(
            mock_post,
            self.api_extractor_a,
            self.response,
            pl.DataFrame(
                {
                    "pos": [100, 150],
                    "alleles": ["A/G/C", "A/G"],
                },
                schema=pl.Schema(
                    {
                        "pos": pl.UInt32,
                        "alleles": pl.Utf8,
                    },
                ),
            ),
        )
        mock_post.reset_mock()
        self.check_vep_api_dataframe(
            mock_post,
            self.api_extractor_b,
            self.response,
            pl.DataFrame(
                {
                    "pos": [100, 150],
                    "gene_name": ["Godzilla", "Mothra"],
                },
                schema=pl.Schema(
                    {
                        "pos": pl.UInt32,
                        "gene_name": pl.Utf8,
                    },
                ),
            ),
        )

    @mock.patch("httpx.Client.post")
    def test_vep_api_dataframe_empty(self, mock_post):
        """Test the high-level varona.vep_api_dataframe() function.

        Test with empty response, no rows.
        """

        self.check_vep_api_dataframe(
            mock_post,
            self.api_extractor_a,
            [],
            pl.DataFrame(
                {
                    "pos": [],
                    "alleles": [],
                },
                schema=pl.Schema(
                    {
                        "pos": pl.UInt32,
                        "alleles": pl.Utf8,
                    },
                ),
            ),
        )
        mock_post.reset_mock()
        self.check_vep_api_dataframe(
            mock_post,
            self.api_extractor_b,
            [],
            pl.DataFrame(
                {
                    "pos": [],
                    "gene_name": [],
                },
                schema=pl.Schema(
                    {
                        "pos": pl.UInt32,
                        "gene_name": pl.Utf8,
                    },
                ),
            ),
        )
