"""Minor allele frequency (MAF) calculation.

This module provides a few competing strategies for obtaining MAF data from
a VFC file.

"""

import collections

import pysam

from varona import bcftools, enum


class MafMethod(enum.CiStrEnum):
    """Enum for the MAF calculation method.

    Values:

    +-----------+----------------------------------------------+
    | Value     | Description                                  |
    +===========+==============================================+
    | FR        | Use the FR field in the info section of the  |
    |           | variant to calculate the MAF.                |
    +-----------+----------------------------------------------+
    | BCFTOOLS  | Only available if non-pysam bcftools is      |
    |           | installed. Use the MAF                       |
    |           |                                              |
    |           | field in the info section of the variant to  |
    |           | calculate the MAF, which                     |
    |           |                                              |
    |           | is added by ``bcftools +fill-tags``.         |
    +-----------+----------------------------------------------+
    | SAMPLES   | Description of SAMPLES method                |
    +-----------+----------------------------------------------+

    Note: The BCFTOOLS value is conditionally available based on the
    presence of non-pysam bcftools.
    """

    FR = enum.auto()
    # Conditionally add the INFO method if non-pysam bcftools
    # is available.
    if bcftools.HAVE_BCFTOOLS:
        BCFTOOLS = enum.auto()
    SAMPLES = enum.auto()


def maf_from_fr(variant: pysam.VariantRecord) -> float:
    """Compute the MAF from the "FR" info value in a variant.

    :param variant: A :class:`pysam.VariantRecord` object.
    :return: The MAF value derived from the "FR" info value.
    :raises: :class:`KeyError` if the "FR" key is not found in the info section.
    """
    # the FR values focus on the alt alleles but we need the
    # reference allele frequency too
    af_list = list(variant.info["FR"])
    ref_af = 1 - sum(af_list)
    af_list.append(ref_af)
    af_list.sort(reverse=True)
    # it'll be the second-highest
    return af_list[1]


def maf_from_info(variant: pysam.VariantRecord) -> float:
    """Extract the "MAF" value from the info section in a variant.

    This function assumes the presence of an "MAF" key in the info section.
    `bcftools <https://samtools.github.io/bcftools/bcftools.html>`_, if
    installed, may be used to add this key to the VCF file using a command
    like ``bcftools +fill-tags``.

    :param variant: A :class:`pysam.VariantRecord` object.
    :return: The MAF value.
    :raises: KeyError if the "MAF" key is not found in the info section.
    """
    return variant.info["MAF"]


def maf_from_samples(variant: pysam.VariantRecord) -> float:
    """Compute the MAF from the sample genotypes in a variant.

    This function assumes that the genotypes are in the "GT" key of the
    samples section of the variant.

    :param variant: A :class:`pysam.VariantRecord` object.
    :return: The MAF value.
    :raises: KeyError if the "GT" key is not found in the samples section.
    """
    var_allele_count = len(variant.alleles)
    samp_allele_counts = collections.Counter()
    # for some reason, iterating over samples just gives the sample name
    for ix in range(len(variant.samples)):
        samp_allele_counts.update(variant.samples[ix]["GT"])
    num_alleles = samp_allele_counts.total()
    af_list = [
        float(samp_allele_counts.get(i, 0)) / float(num_alleles)
        for i in range(var_allele_count)
    ]
    af_list.sort(reverse=True)
    return af_list[1]


def maf_from_method(variant: pysam.VariantRecord, method: MafMethod) -> float:
    """A dispatcher for :func:`maf_from_fr`, :func:`maf_from_info`, and :func:`maf_from_samples`.

    :param variant: A :class:`pysam.VariantRecord` object.
    :param method: The MAF calculation method.
    :return: The MAF value.
    :raises: KeyError if the method is not recognized.
    """
    if method == MafMethod.FR:
        return maf_from_fr(variant)
    elif method == MafMethod.BCFTOOLS:
        return maf_from_info(variant)
    return maf_from_samples(variant)
