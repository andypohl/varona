"""Tests for the platypus module."""

from unittest import mock

import polars as pl
import polars.testing as plt

from varona import cli, ensembl, fake_vcf, maf, platypus


class TestPlatypusDataFrame(fake_vcf.TestWithTempDir):
    """Tests the :func:`varona.platypus.platypus_dataframe` function works."""

    api_response = [
        {
            "assembly_name": "GRCh37",
            "variant_class": "SNV",
            "most_severe_consequence": "synonymous_variant",
            "allele_string": "A/G",
            "transcript_consequences": [
                {
                    "cdna_end": 833,
                    "impact": "LOW",
                    "protein_end": 190,
                    "protein_start": 190,
                    "variant_allele": "G",
                    "strand": -1,
                    "transcript_id": "ENST00000360001",
                    "hgnc_id": 24188,
                    "amino_acids": "D",
                    "cds_start": 570,
                    "gene_id": "ENSG00000078808",
                    "cds_end": 570,
                    "biotype": "protein_coding",
                    "gene_symbol_source": "HGNC",
                    "consequence_terms": ["synonymous_variant"],
                    "gene_symbol": "SDF4",
                    "codons": "gaT/gaC",
                    "cdna_start": 833,
                }
            ],
            "seq_region_name": "1",
            "end": 1158631,
            "start": 1158631,
            "colocated_variants": [
                {
                    "id": "COSV55420653",
                    "somatic": 1,
                    "var_synonyms": {"COSMIC": ["COSM3750257"]},
                    "start": 1158631,
                    "end": 1158631,
                    "allele_string": "COSMIC_MUTATION",
                    "seq_region_name": "1",
                    "strand": 1,
                    "phenotype_or_disease": 1,
                },
                {
                    "seq_region_name": "1",
                    "allele_string": "A/G/T",
                    "end": 1158631,
                    "start": 1158631,
                    "strand": 1,
                    "frequencies": {
                        "G": {
                            "gnomade_nfe": 0.8879,
                            "afr": 0.9985,
                            "gnomade_eas": 0.9998,
                            "gnomade_afr": 0.9834,
                            "gnomade_oth": 0.9053,
                            "gnomade_asj": 0.8568,
                            "gnomade": 0.9144,
                            "af": 0.9533,
                            "gnomade_amr": 0.9508,
                            "gnomade_fin": 0.8519,
                            "amr": 0.9294,
                            "eur": 0.8608,
                            "sas": 0.956,
                            "gnomade_sas": 0.9485,
                            "eas": 1,
                        }
                    },
                    "id": "rs6603781",
                },
            ],
            "id": ".",
            "input": "1 1158631 . A G . . .",
            "strand": 1,
        },
        {
            "end": 6219293,
            "transcript_consequences": [
                {
                    "hgnc_id": 16816,
                    "transcript_id": "ENST00000262450",
                    "strand": -1,
                    "variant_allele": "CACA",
                    "consequence_terms": ["intron_variant"],
                    "gene_symbol": "CHD5",
                    "gene_symbol_source": "HGNC",
                    "biotype": "protein_coding",
                    "impact": "MODIFIER",
                    "gene_id": "ENSG00000116254",
                }
            ],
            "seq_region_name": "1",
            "allele_string": "CACACA/CACA/-",
            "start": 6219288,
            "assembly_name": "GRCh37",
            "variant_class": "sequence_alteration",
            "most_severe_consequence": "intron_variant",
            "input": "1 6219287 . TCACACA TCACA,T . . .",
            "id": ".",
            "strand": 1,
            "colocated_variants": [
                {
                    "strand": 1,
                    "start": 6219288,
                    "allele_string": "CACACACACACACACACACACACACACACACA/CACACACACACACA/CACACACACACACACA/CACACACACACACACACA/CACACACACACACACACACA/CACACACACACACACACACACA/CACACACACACACACACACACACA/CACACACACACACACACACACACACA/CACACACACACACACACACACACACACA/CACACACACACACACACACACACACACACA/CACACACACACACACACACACACACACACACACA/CACACACACACACACACACACACACACACACACACA/CACACACACACACACACACACACACACACACACACACA/CACACACACACACACACACACACACACACACACACACACA",
                    "seq_region_name": "1",
                    "end": 6219319,
                    "id": "rs56933510",
                }
            ],
        },
    ]

    @mock.patch("httpx.Client.post")
    def test_platypus_dataframe(self, mock_post):
        """Tests the :func:`varona.platypus.platypus_dataframe` function."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.add_records_from_lines(
            """\
1\t1158631\t.\tA\tG\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-43.88,0.0:3:99:160:156
1\t6219287\t.\tTCACACA\tTCACA,T\t2997\tPASS\tBRF=0.35;FR=0.5000,0.5000;HP=1;HapScore=1;MGOF=6;MMLQ=23;MQ=51.91;NF=108,108;NR=0,0;PP=2997,2652;QD=20;SC=TGAGACTCCATCACACACACA;SbPval=1.0;Source=Platypus;TC=172;TCF=170;TCR=2;TR=108,108;WE=6219302;WS=6219277\tGT:GL:GOF:GQ:NR:NV\t1/2:-1,-1,-1:6:99:172,172:108,108
"""
        )
        fake.write_vcf()
        mock_200 = mock.MagicMock()
        mock_200.status_code = 200
        mock_200.json.return_value = self.api_response
        mock_post.return_value = mock_200
        df = platypus.platypus_dataframe(fake.path)
        expected_df = pl.DataFrame(
            {
                "contig": ["1", "1"],
                "pos": [1158631, 6219287],
                "ref": ["A", "TCACACA"],
                "alt": ["G", "TCACA,T"],
                "sequence_depth": [160, 172],
                "max_variant_reads": [156, 108],
                "variant_read_pct": [97.5, 62.7906976744186],
                "maf": [None, None],
                "type": ["SNV", None],
                "effect": ["synonymous_variant", None],
                "gene_name": ["SDF4", None],
                "gene_id": ["ENSG00000078808", None],
                "transcript_id": ["ENST00000360001", None],
            },
            schema=platypus.VCF_DF_SCHEMA | platypus.API_DF_SCHEMA,
        )
        plt.assert_frame_equal(df, expected_df)
