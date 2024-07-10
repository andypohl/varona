"""Module with functions for extracting data from lower-level structures.

This module contains functions that can be used to extract data from
lower-level structures like VCF records and VEP API responses (dictionaries).
These functions are called from the Verona top-level functions in the
:mod:`verona.verona` module and can be replaced with custom functions 
when using those, if desired.
"""

import logging
import typing

import pysam

logger = logging.getLogger("verona.extract")


def default_vep_response_extractor(response_item: dict) -> dict:
    """An example function to extract VEP data from the response.

    Care should be taken to handle the data in the response with correct
    spelling of keys, indexing, etc. :class:`KeyError` and :class:`IndexError`
    exceptions can either be handled by this function or by the caller to
    :func:`verona.ensemble.query_vep_api`.  For this specific example, no
    exceptions are caught and so will cause the query to fail.

    Input item (example):

    .. code-block:: json

        {
        }

    Output item (example):

    .. code-block:: json

        {
        }

    :param response_item: The VEP API response is a list of dictionaries,
        and this is one of the items in the list.
    :return: A dictionary with the VEP data transformed into a preferable format.
    """
    new_item = {}
    new_item["contig"] = response_item["seq_region_name"]
    new_item["pos"] = response_item["start"]
    alleles = response_item["allele_string"].split("/")
    new_item["ref"] = alleles[0]
    new_item["alt"] = ",".join(alleles[1:])
    new_item["type"] = response_item["variant_class"]
    new_item["effect"] = response_item["most_severe_consequence"]
    if "transcript_consequences" not in response_item:
        new_item["gene_name"] = None
        new_item["gene_id"] = None
        new_item["transcript_id"] = None
    else:
        new_item["gene_name"] = response_item["transcript_consequences"][0][
            "gene_symbol"
        ]
        new_item["gene_id"] = response_item["transcript_consequences"][0]["gene_id"]
        new_item["transcript_id"] = response_item["transcript_consequences"][0][
            "transcript_id"
        ]
    return new_item


def platypus_vcf_record_extractor(
    record: pysam.VariantRecord,
    **addl_cols: typing.Callable[[pysam.VariantRecord], typing.Union[int, float, str]]
) -> dict:
    """Example function to extract data from a VCF record."""
    new_item = {}
    new_item["contig"] = record.contig
    new_item["pos"] = record.pos
    new_item["ref"] = record.ref
    new_item["alt"] = ",".join(record.alts)
    new_item["sequence_depth"] = record.info["TC"]  # in some VCFs this is "DP"
    # for multiallelic TR is a list of the number of reads supporting each allele
    # so the max is taken.
    new_item["max_variant_reads"] = max(record.info["TR"])
    new_item["variant_read_pct"] = (
        new_item["max_variant_reads"] / new_item["sequence_depth"] * 100
    )
    for col, func in addl_cols.items():
        new_item[col] = func(record)
    return new_item
