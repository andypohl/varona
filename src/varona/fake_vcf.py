"""Module for fake VCF file creation.

The classes in this module provide quick VCF file creation for testing purposes.
In one case, the file mimics a Platypus VCF file, which has an older VCF format
(v4.0).  In the other case, the file is a generic VCF file with a v4.2 format.

.. note::

    As of pysam v0.22.1 (April 23, 2024), the VCF header doesn't capture the
    fileformat line correctly sometimes, so when mimicking a Platypus VCF file,
    the file is rewritten with a Platypus-like header and the fileformat
    is set to v4.0. There is a regression for pysam in case that behvaior
    changes in the future, in which case the `rewrite=True` parameter can perhaps
    be flipped to `False`.
    
"""

import importlib.resources as pkg_resources
import pathlib
import shutil
import tempfile
import unittest

import pysam


class FakeVcfFile:
    """Make a minimal fake VCF files for testing.

    Unlike the :class:`pysam.AlignedSegment` class, the :class:`pysam.VariantRecord`
    class lacks `from_dict` and `fromstring` convenience methods. This class
    provides a way to create and add to minimal VCF files for testing.
    """

    def __init__(
        self,
        path: pathlib.Path,
        samples: list[str] = ["sample"],
        records: list[dict] = [],
    ):
        self.path = path
        self.make_header()
        self.header.add_samples(samples)
        self.records = []
        self.add_records(records)

    def make_header(self):
        """Add a basic header to the VCF file."""
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
                ("ID", "GT"),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "Genotype"),
            ],
        )
        self.header = header

    def add_records(self, records: list[dict]):
        """Add records to the VCF file."""
        for record in records:
            if "pos" in record:
                record["start"] = record.pop("pos")
                record["stop"] = record["start"]
                record["start"] -= 1
            record.setdefault("contig", "1")
            record.setdefault("id", ".")
            record.setdefault("qual", 100)
            record.setdefault("filter", "PASS")
            self.records.append(self.header.new_record(**record))

    def write_vcf(self):
        """Write the VCF file with the header and records."""
        with pysam.VariantFile(self.path, "w", header=self.header) as vcf:
            for record in self.records:
                vcf.write(record)


class FakePlatypusVcfFile(FakeVcfFile):
    """Make a fake Platypus VCF file for testing.

    Platypus VCFs have a different header format than generic VCFs.
    """

    header_preamble = """\
##fileformat=VCFv4.0
##fileDate=2016-06-21
##source=Platypus_Version_0.8.1
##platypusOptions={'assemblyRegionSize': 1500, 'trimReadFlank': 0, 'assembleBadReads': 1, 'bamFiles': ['/data/TCGA-A2-A0YC/tumor.bam'], 'minVarDist': 9, 'trimSoftClipped': 1, 'minReads': 2, 'qualBinSize': 1, 'refFile': '/data/ref_genome/human_g1k_v37.fasta', 'maxHaplotypes': 50, 'filterVarsByCoverage': 1, 'maxSize': 1500, 'originalMaxHaplotypes': 50, 'skipDifficultWindows': 0, 'parseNCBI': 0, 'skipRegionsFile': None, 'noCycles': 0, 'trimAdapter': 1, 'minPosterior': 5, 'assembleAll': 1, 'trimOverlapping': 1, 'filterDuplicates': 1, 'abThreshold': 0.001, 'minFlank': 10, 'bufferSize': 100000, 'fileCaching': 0, 'useEMLikelihoods': 0, 'coverageSamplingLevel': 30, 'calculateFlankScore': 0, 'logFileName': 'log.txt', 'nCPU': 1, 'filterReadsWithUnmappedMates': 1, 'qdThreshold': 10, 'maxVariants': 8, 'scThreshold': 0.95, 'filterReadsWithDistantMates': 1, 'maxReads': 5000000, 'badReadsWindow': 11, 'genIndels': 1, 'largeWindows': 0, 'minMapQual': 20, 'maxVarDist': 15, 'maxGOF': 30, 'rlen': 150, 'minGoodQualBases': 20, 'refCallBlockSize': 1000, 'countOnlyExactIndelMatches': 0, 'longHaps': 0, 'HLATyping': 0, 'filterReadPairsWithSmallInserts': 1, 'minBaseQual': 20, 'getVariantsFromBAMs': 1, 'genSNPs': 1, 'assemble': 0, 'assemblerKmerSize': 15, 'minVarFreq': 0.05, 'alignScoreFile': '', 'verbosity': 2, 'sourceFile': None, 'compressReads': 0, 'rmsmqThreshold': 40, 'filteredReadsFrac': 0.7, 'outputRefCalls': 0, 'badReadsThreshold': 15, 'hapScoreThreshold': 4, 'regions': None, 'sbThreshold': 0.001, 'output': '/data/TCGA-A2-A0YC/platypus_out.vcf', 'assembleBrokenPairs': 0, 'mergeClusteredVariants': 1, 'maxGenotypes': 1275, 'nInd': 1}
##filter="MQ > 50 & TC > 100 & QUAL > 2900"
##INFO=<ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">
##INFO=<ID=MMLQ,Number=1,Type=Float,Description="Median minimum base quality for bases around variant">
##INFO=<ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">
##INFO=<ID=HP,Number=1,Type=Integer,Description="Homopolymer run length around variant locus">
##INFO=<ID=WE,Number=1,Type=Integer,Description="End position of calling window">
##INFO=<ID=Source,Number=.,Type=String,Description="Was this variant suggested by Playtypus, Assembler, or from a VCF?">
##INFO=<ID=FS,Number=.,Type=Float,Description="Fisher's exact test for strand bias (Phred scale)">
##INFO=<ID=WS,Number=1,Type=Integer,Description="Starting position of calling window">
##INFO=<ID=PP,Number=.,Type=Float,Description="Posterior probability (phred scaled) that this variant segregates">
##INFO=<ID=TR,Number=.,Type=Integer,Description="Total number of reads containing this variant">
##INFO=<ID=NF,Number=.,Type=Integer,Description="Total number of forward reads containing this variant">
##INFO=<ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">
##INFO=<ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">
##INFO=<ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">
##INFO=<ID=END,Number=.,Type=Integer,Description="End position of reference call block">
##INFO=<ID=MGOF,Number=.,Type=Integer,Description="Worst goodness-of-fit value reported across all samples">
##INFO=<ID=SbPval,Number=.,Type=Float,Description="Binomial P-value for strand bias test">
##INFO=<ID=START,Number=.,Type=Integer,Description="Start position of reference call block">
##INFO=<ID=ReadPosRankSum,Number=.,Type=Float,Description="Mann-Whitney Rank sum test for difference between in positions of variants in reads from ref and alt">
##INFO=<ID=MQ,Number=.,Type=Float,Description="Root mean square of mapping qualities of reads at the variant position">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">
##INFO=<ID=SC,Number=1,Type=String,Description="Genomic sequence 10 bases either side of variant position">
##INFO=<ID=BRF,Number=1,Type=Float,Description="Fraction of reads around this variant that failed filters">
##INFO=<ID=HapScore,Number=.,Type=Integer,Description="Haplotype score measuring the number of haplotypes the variant is segregating into in a window">
##INFO=<ID=Size,Number=.,Type=Integer,Description="Size of reference call block">
##FILTER=<ID=GOF,Description="Variant fails goodness-of-fit test.">
##FILTER=<ID=badReads,Description="Variant supported only by reads with low quality bases close to variant position, and not present on both strands.">
##FILTER=<ID=alleleBias,Description="Variant frequency is lower than expected for het">
##FILTER=<ID=hp10,Description="Flanking sequence contains homopolymer of length 10 or greater">
##FILTER=<ID=Q20,Description="Variant quality is below 20.">
##FILTER=<ID=HapScore,Description="Too many haplotypes are supported by the data in this region.">
##FILTER=<ID=MQ,Description="Root-mean-square mapping quality across calling region is low.">
##FILTER=<ID=strandBias,Description="Variant fails strand-bias filter">
##FILTER=<ID=SC,Description="Variants fail sequence-context filter. Surrounding sequence is low-complexity">
##FILTER=<ID=QualDepth,Description="Variant quality/Read depth ratio is low.">
##FILTER=<ID=REFCALL,Description="This line represents a homozygous reference call">
##FILTER=<ID=QD,Description="Variants fail quality/depth filter.">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
##FORMAT=<ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">
##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">\
"""

    def make_header(self):
        """Add a Platypus header to the VCF file."""
        header = pysam.VariantHeader()
        for line in self.header_preamble.split("\n"):
            header.add_line(line)
        with pkg_resources.open_text("varona.data", "human_g1k_v37.fasta.fai") as f:
            for line in f:
                name, length = line.split("\t")[:2]
                header.contigs.add(name, length)
        self.header = header

    def add_records_from_lines(self, lines):
        """Add records to the VCF file from a list of lines.

        These lines are expected to be in the format of a Platypus VCF file and
        may be copy/pasted from a real file.

        One note is that the this method writes a temporary file for the sake
        of parsing the records into memory. This is a bit silly, but it works.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            tmp_path = pathlib.Path(tempdir) / "temp.vcf"
            cur_records = self.records
            self.records = []
            self.write_vcf(rewrite=True)
            shutil.copy(self.path, tmp_path)
            self.path.unlink()
            with tmp_path.open("a", encoding="utf-8") as f:
                f.write(lines)
            with pysam.VariantFile(str(tmp_path)) as vcf:
                self.records.extend(list(vcf))

    def write_vcf(self, rewrite=True):
        """Write the VCF file with the header and records.

        The Platypus VCF file has a slightly different header format than generic
        VCF files written by pysam. The VCFv4.0 format doesn't have contig information
        mandatory.

        In the Platypus VCF, pysam seems to ignore the fileformat line, and otherwise writes
        out in a different order as the test file.

        :param rewrite: If True, rewrite the file to ensure the header is correct.
        """
        super().write_vcf()
        if rewrite:
            tmp_path = self.path.with_suffix(".tmp")
            with open(self.path, "r") as infile, open(tmp_path, "w") as outfile:
                # skip the fileformat line and the PASS filter line pysam puts at the top
                for _ in range(2):
                    next(infile)
                while True:
                    line = next(infile)
                    if line.startswith("#CHROM"):
                        break
                    elif line.startswith("##contig"):
                        continue
                    outfile.write(line)
                # don't forget the CHROM line
                outfile.write(line)
                # write the rest of the file
                for line in infile:
                    outfile.write(line)
            tmp_path.replace(self.path)


class TestWithTempDir(unittest.TestCase):
    """Boilerplate for tests that need a temporary directory."""

    def setUp(self):
        self.tmp_dir = pathlib.Path(tempfile.mkdtemp())
        self.path = self.tmp_dir / "test.vcf"

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)
