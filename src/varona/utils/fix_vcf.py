"""VCF file utilities, mainly to fix problems.
"""

import importlib.resources as pkg_resources

import pysam

from varona import ensembl


def vcf_header_inject_contigs(
    header: pysam.VariantHeader, assembly: ensembl.Assembly
) -> pysam.VariantHeader:
    """Inject contigs into a VCF header.

    Some VCF files do not have contigs in the header. This causes problems with
    pysam when trying to write out a new VCF. This function shoe-horns contigs
    into the headers.

    .. warning::

        This function only works for human assemblies GRCh37 and GRCh38, and
        for GRCh37 it uses the 1000 genomes project contigs. For that reason,
        non-1000-genome-project-GRCh37 VCF files are not quite compatible at
        this time.

    :param header: The VCF header.
    :return: The VCF header with contigs injected.
    """
    ret = header.copy()
    if assembly == ensembl.Assembly.GRCH37:
        data_path = "human_g1k_v37.fasta.fai"
    elif assembly == ensembl.Assembly.GRCH38:
        data_path = "human_grch38.fasta.fai"
    else:
        raise ValueError(f"Unknown assembly: {assembly}")
    with pkg_resources.open_text("varona.data", "human_g1k_v37.fasta.fai") as f:
        contigs = {}
        for line in f:
            name, length = line.split("\t")[:2]
            contigs[name] = int(length)
        # check if the contigs are already in the header
        for name, length in contigs.items():
            if name in ret.contigs:
                if ret.contigs[name] != length:
                    raise ValueError(
                        f"Contig {name} already in header with different length"
                    )
            else:
                ret.contigs.add(name, length)
    return ret
