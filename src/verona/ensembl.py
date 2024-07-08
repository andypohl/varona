"""Module for querying Ensembl API.

This module also contains code to read the VCF file and prepare
the data for querying the Ensembl API.
"""

import contextlib
import csv
import io
import logging
import pathlib
import typing

import pysam

HTS_COMPRESSED_VALS = ("GZIP", "BGZF")
"""Possible values for the compression field in a :class:`pysam.VariantFile`.
"""

VCF_MASK_COLS = [2, 5, 6, 7]
"""Columns to mask in the VCF file.

Values in these columns are replaced with '.' before querying the
Ensembl API.
"""


@contextlib.contextmanager
def _text_from_bgzf(file_path, mode="r"):
    """Wraps a UTF-8 byte decoder around a BGZF-decompressor.

    This is just a helper to open .vcf and .vcf.gz files in the
    """
    with pysam.BGZFile(file_path, mode) as bgzf, io.TextIOWrapper(
        bgzf, encoding="utf-8"
    ) as text:
        yield text


def vcf_rows(vcf_path: pathlib.Path) -> typing.Iterator[list[str]]:
    """Yields rows from a VCF file.

    This function opens a VCF file that's either plain text or compressed.
    Normally the `:class:pysam.VariantFile` class would be used to read
    VCF files but in this case there's no need to parse the VCF fully as
    VariantFile does, it's just a matter of reading tab-separated lines
    of text.

    This function could be moved to a more general module in the future,
    but currently it's only needed by the Ensembl API querying code.

    :param vcf_path: Path to the VCF file.
    """
    file_opener = open
    with pysam.VariantFile(vcf_path, "r") as vf:
        if not vf.is_vcf:
            raise ValueError(f"{vcf_path} does not appear to be a VCF file.")
        if vf.compression in HTS_COMPRESSED_VALS:
            file_opener = _text_from_bgzf
    with file_opener(vcf_path, "r") as vcf:
        reader = csv.reader(vcf, delimiter="\t")
        for row in reader:
            if row and len(row) >= 1 and row[0].startswith("#"):
                continue
            yield row


def get_vcf_query_data(
    vcf_path: pathlib.Path, chunk_size=1000
) -> typing.Iterator[list[str]]:
    """Gets lists of data in chunks from the VCF file for querying the Ensembl API.

    Specifically, the `POST vep/:species/region <https://rest.ensembl.org/documentation/info/vep_region_post>`_
    endpoint.
    """
    chunk = []
    for row in vcf_rows(vcf_path):
        for col in VCF_MASK_COLS:
            row[col] = "."
        chunk.append(" ".join(row[:8]))
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
