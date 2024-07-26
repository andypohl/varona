import unittest

import pysam

from varona import extract, fake_vcf, maf

example_response = {
    "end": 1158631,
    "transcript_consequences": [
        {
            "impact": "LOW",
            "cdna_end": 833,
            "gene_symbol": "SDF4",
            "consequence_terms": ["synonymous_variant"],
            "cds_start": 570,
            "amino_acids": "D",
            "protein_end": 190,
            "hgnc_id": 24188,
            "codons": "gaT/gaC",
            "variant_allele": "G",
            "biotype": "protein_coding",
            "strand": -1,
            "cdna_start": 833,
            "transcript_id": "ENST00000360001",
            "gene_symbol_source": "HGNC",
            "gene_id": "ENSG00000078808",
            "protein_start": 190,
            "cds_end": 570,
        }
    ],
    "colocated_variants": [
        {
            "somatic": 1,
            "strand": 1,
            "var_synonyms": {"COSMIC": ["COSM3750257"]},
            "allele_string": "COSMIC_MUTATION",
            "phenotype_or_disease": 1,
            "id": "COSV55420653",
            "seq_region_name": "1",
            "end": 1158631,
            "start": 1158631,
        },
        {
            "start": 1158631,
            "id": "rs6603781",
            "seq_region_name": "1",
            "end": 1158631,
            "allele_string": "A/G/T",
            "strand": 1,
            "frequencies": {
                "G": {
                    "af": 0.9533,
                    "afr": 0.9985,
                    "gnomade_asj": 0.8568,
                    "eur": 0.8608,
                    "gnomade_eas": 0.9998,
                    "gnomade_amr": 0.9508,
                    "gnomade_fin": 0.8519,
                    "gnomade_sas": 0.9485,
                    "gnomade": 0.9144,
                    "eas": 1,
                    "amr": 0.9294,
                    "gnomade_nfe": 0.8879,
                    "sas": 0.956,
                    "gnomade_afr": 0.9834,
                    "gnomade_oth": 0.9053,
                }
            },
        },
    ],
    "seq_region_name": "1",
    "id": ".",
    "variant_class": "SNV",
    "most_severe_consequence": "synonymous_variant",
    "input": "1 1158631 . A G . . .",
    "allele_string": "A/G",
    "strand": 1,
    "assembly_name": "GRCh37",
    "start": 1158631,
}


class TestExtractApi(unittest.TestCase):
    """Testing the extraction API functions.

    Currenly only for :func:`varona.extract.default_vep_response_extractor`.
    """

    def test_extract_example_responses(self):
        """Tests the :func:`varona.extract.default_vep_response_extractor` function.

        It's kind of a long copy/paste for test code so it's just one item.
        """
        item = example_response.copy()
        extracted = extract.default_vep_response_extractor(item)
        expected = {
            "contig": "1",
            "pos": 1158631,
            "ref": "A",
            "alt": "G",
            "type": "SNV",
            "effect": "synonymous_variant",
            "gene_name": "SDF4",
            "gene_id": "ENSG00000078808",
            "transcript_id": "ENST00000360001",
        }
        self.assertDictEqual(extracted, expected)
        # delete "transcript_consequences" key
        del item["transcript_consequences"]
        expected["gene_id"] = None
        expected["gene_name"] = None
        expected["transcript_id"] = None
        extracted = extract.default_vep_response_extractor(item)
        self.assertDictEqual(extracted, expected)


class TestExtractVcf(fake_vcf.TestWithTempDir):
    """Testing VCF extractor functions.

    Currently only Platypus-style VCFs are tested.
    """

    def test_extract_platypus_vcf(self):
        """Tests the :func:`varona.extract.platypus_vcf_record_extractor` function."""
        fake = fake_vcf.FakePlatypusVcfFile(self.path)
        fake.add_records_from_lines(
            """\
1\t1158631\t.\tA\tG\t2965\tPASS\tBRF=0.16;FR=1.0000;HP=1;HapScore=1;MGOF=3;MMLQ=33;MQ=59.75;NF=89;NR=67;PP=2965;QD=20;SC=CACTTTCCTCATCCACTTTGA;SbPval=0.58;Source=Platypus;TC=160;TCF=90;TCR=70;TR=156;WE=1158639;WS=1158621\tGT:GL:GOF:GQ:NR:NV\t1/1:-300.0,-43.88,0.0:3:99:160:156
1\t6219287\t.\tTCACACA\tTCACA,T\t2997\tPASS\tBRF=0.35;FR=0.5000,0.5000;HP=1;HapScore=1;MGOF=6;MMLQ=23;MQ=51.91;NF=108,108;NR=0,0;PP=2997,2652;QD=20;SC=TGAGACTCCATCACACACACA;SbPval=1.0;Source=Platypus;TC=172;TCF=170;TCR=2;TR=108,108;WE=6219302;WS=6219277\tGT:GL:GOF:GQ:NR:NV\t1/2:-1,-1,-1:6:99:172,172:108,108
"""
        )
        fake.write_vcf()
        with pysam.VariantFile(fake.path) as vcf:
            records = list(vcf)
            self.assertEqual(len(records), 2)
            extracted = [
                extract.platypus_vcf_record_extractor(
                    record,
                    maf=maf.maf_from_samples,
                    maf_from_fr=maf.maf_from_fr,
                )
                for record in records
            ]
            self.assertAlmostEqual(extracted[0]["maf"], 0)
            self.assertAlmostEqual(extracted[0]["maf_from_fr"], 0)
            self.assertEqual(extracted[0]["sequence_depth"], 160)
            self.assertEqual(extracted[0]["max_variant_reads"], 156)
            self.assertAlmostEqual(extracted[0]["variant_read_pct"], 97.5)
            self.assertEqual(extracted[1]["alt"], "TCACA,T")
            self.assertAlmostEqual(extracted[1]["maf"], 0.5)
            self.assertEqual(extracted[1]["sequence_depth"], 172)
            self.assertEqual(extracted[1]["max_variant_reads"], 108)
            self.assertAlmostEqual(extracted[1]["variant_read_pct"], 62.7906976744186)
