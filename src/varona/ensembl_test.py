import os
import unittest
from unittest import mock

import httpx
import pysam

from varona import ensembl, extract, fake_vcf

VARONA_LIVE_TESTS = os.getenv("VARONA_LIVE_TESTS", "0") == "1"
"""Set to 1 to enable live querying tests for local development.
"""


class TestQueryVepApi(unittest.TestCase):
    """Tests the :func:`varona.ensembl.query_vep_api` function.

    Still incomplete. More tests are needed.  Needed tests:

    - Code 400 (Bad Request) responses
    - Code 50x responses
    - Code 429 responses without Retry-After header

    """

    def setUp(self):
        self.client = httpx.Client()

    def tearDown(self):
        self.client.close()

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
    response_200 = httpx.Response(status_code=200, json=[{"data": "success"}])

    @mock.patch("time.sleep", return_value=None)  # disable sleeping
    def test_query_and_retries(self, mock_sleep):
        with mock.patch.object(self.client, "post") as mock_post:
            mock_post.side_effect = [self.response_429, self.response_200]
            # Call the function
            with self.assertLogs("varona.ensembl", level="WARNING") as cm:
                result = list(ensembl.query_vep_api(self.client, chunk=self.chunk))
                self.assertEqual(len(cm.output), 1)  # Verify warning message
                self.assertEqual(
                    cm.output[0],
                    "WARNING:varona.ensembl:API code 429 with Retry-After. Retrying after 1 seconds.",
                )
            self.assertListEqual(result, [{"data": "success"}])  # Verify the result
            self.assertEqual(mock_post.call_count, 2)  # Verify two API calls (retry)
            mock_sleep.assert_called_once_with(1)  # Verify retry delay

    @mock.patch("time.sleep", return_value=None)  # disable sleeping
    def test_too_many_retries(self, _):
        with mock.patch.object(self.client, "post") as mock_post:
            mock_post.side_effect = [
                self.response_429,  # mock response 1
                self.response_429,  # mock response 2
                self.response_429,  # mock response 3
                self.response_200,  # (doesn't get this far)
            ]
            # Call the function
            with self.assertLogs("varona.ensembl", level="WARNING") as cm:
                with self.assertRaises(TimeoutError):
                    list(
                        ensembl.query_vep_api(self.client, chunk=self.chunk, retries=2)
                    )
                self.assertEqual(len(cm.output), 3)  # Verify warning message
                for i in range(2):
                    self.assertEqual(
                        cm.output[i],
                        "WARNING:varona.ensembl:API code 429 with Retry-After. Retrying after 1 seconds.",
                    )
            self.assertEqual(mock_post.call_count, 3)


class TestChunkReader(fake_vcf.TestWithTempDir):
    """Tests the :func:`verona.ensembl.vcf_to_vep_query_data` function works."""

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


@unittest.skipUnless(VARONA_LIVE_TESTS, "Live querying tests are disabled by default.")
class TestLiveQuerying(unittest.TestCase):
    """Tests the live Ensembl API.

    Normally these should be disabled. They are here for development purposes and
    require opting-in by setting the environment variable `VARONA_LIVE_TESTS=1`.
    """

    def setUp(self):
        self.client = httpx.Client()

    def tearDown(self):
        self.client.close()

    def test_two_records(self):
        """Tests querying the Ensembl API with two records."""
        chunk = [
            "1 1158631 . A G . . .",  # biallelic
            "1 91859795 . TATGTGA CATGTGA,CATGTGG . . .",  # > 2 alleles
        ]
        data = list(
            ensembl.query_vep_api(
                self.client,
                chunk,
                response_extractor=extract.default_vep_response_extractor,
            )
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


class TestJsonInput(fake_vcf.TestWithTempDir):

    json_lines = """\
{"end":1158631,"start":1158631,"nearest":["SDF4"],"input":"1\t1158631\t.\tA\tG\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-43.88,0.0:3:99:160:156","strand":1,"most_severe_consequence":"synonymous_variant","seq_region_name":"1","assembly_name":"GRCh37","transcript_consequences":[{"strand":-1,"cds_end":570,"protein_start":190,"cdna_end":833,"cdna_start":833,"impact":"LOW","cds_start":570,"transcript_id":"ENST00000360001","variant_allele":"G","codons":"gaT/gaC","gene_id":"ENSG00000078808","consequence_terms":["synonymous_variant"],"protein_end":190,"amino_acids":"D"}],"variant_class":"SNV","allele_string":"A/G","id":"."}
{"nearest":["PUSL1"],"input":"1\t1246004\t.\tA\tG\t2965\tPASS\tBRF=0.09;FR=1.0000;HP=6;HapScore=1;MGOF=5;MMLQ=32;MQ=59.5;NF=101;NR=47;PP=2965;QD=20;SC=ACAGGTACGTATTTTTCCAGG;SbPval=0.62;Source=Platypus;TC=152;TCF=101;TCR=51;TR=148;WE=1246012;WS=1245994\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-41.24,0.0:5:99:152:148","strand":1,"most_severe_consequence":"splice_polypyrimidine_tract_variant","start":1246004,"variant_class":"SNV","allele_string":"A/G","id":".","seq_region_name":"1","assembly_name":"GRCh37","transcript_consequences":[{"variant_allele":"G","transcript_id":"ENST00000379031","strand":1,"impact":"LOW","gene_id":"ENSG00000169972","consequence_terms":["splice_polypyrimidine_tract_variant","intron_variant"]}],"end":1246004}
{"end":1249187,"variant_class":"SNV","id":".","allele_string":"G/A","assembly_name":"GRCh37","seq_region_name":"1","transcript_consequences":[{"amino_acids":"F","protein_end":300,"gene_id":"ENSG00000127054","consequence_terms":["synonymous_variant"],"codons":"ttC/ttT","variant_allele":"A","transcript_id":"ENST00000540437","cds_start":900,"impact":"LOW","cdna_start":1356,"cdna_end":1356,"protein_start":300,"cds_end":900,"strand":-1}],"nearest":["PUSL1"],"input":"1\t1249187\t.\tG\tA\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=3;HapScore=1;MGOF=3;MMLQ=34;MQ=59.72;NF=78;NR=57;PP=2965;QD=20;SC=TCCTCTGCACGAAAGTCTTGC;SbPval=0.53;Source=Platypus;TC=137;TCF=79;TCR=58;TR=135;WE=1249195;WS=1249177\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-37.63,0.0:3:99:137:135","strand":1,"most_severe_consequence":"synonymous_variant","start":1249187}
"""

    def test_json_uncompressed(self):
        """Tests reading JSON lines from a plain file."""
        json_path = self.path.with_suffix(".json")
        with open(json_path, "w") as f:
            f.write(self.json_lines)
        data = list(ensembl.import_vep_data(json_path))
        self.assertEqual(len(data), 3)

    def test_json_bgzipped(self):
        """Tests reading JSON lines from a bgzipped file."""
        json_path = self.path.with_suffix(".json")
        with open(json_path, "w") as f:
            f.write(self.json_lines)
        compressed_path = json_path.with_suffix(".json.gz")
        pysam.tabix_compress(str(json_path), str(compressed_path))
        data = list(ensembl.import_vep_data(compressed_path))
        self.assertEqual(len(data), 3)

    def test_platypus(self):
        """Tests reading JSON lines from a plain file."""
        json_path = self.path.with_suffix(".json")
        with open(json_path, "w") as f:
            f.write(self.json_lines)
        data = list(
            ensembl.import_vep_data(
                json_path, json_extractor=extract.default_vep_cli_json_extractor
            )
        )
        self.assertEqual(len(data), 3)
