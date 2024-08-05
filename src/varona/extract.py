"""Module with functions for extracting data from lower-level structures.

This module contains functions that can be used to extract data from
lower-level structures like VCF records and VEP API responses (dictionaries).
These functions are called from the Varona top-level functions in the
:mod:`varona.varona` module and can be replaced with custom functions 
when using those, if desired.
"""

import logging
import typing

import pysam

logger = logging.getLogger("varona.extract")


def default_vep_response_extractor(response_item: dict) -> dict:
    """An example function to extract VEP data from the response.

    The response item from the API is a dictionary (not a flat one),
    and this function extracts the data we're interested in
    into a flat dictionary.

    +-------------------+--------------------------------------------+
    | extracted key     | original response_item key                 |
    +===================+============================================+
    | contig            | seq_region_name                            |
    +-------------------+--------------------------------------------+
    | pos               | start                                      |
    +-------------------+--------------------------------------------+
    | ref               | allele_string (first allele, '/' sep)      |
    +-------------------+--------------------------------------------+
    | alt (comma-sep)   | allele_string (all other alleles, '/' sep) |
    +-------------------+--------------------------------------------+
    | type              | variant_class                              |
    +-------------------+--------------------------------------------+
    | effect            | most_severe_consequence                    |
    +-------------------+--------------------------------------------+
    | gene_name         | transcript_consequences[0].gene_symbol     |
    +-------------------+--------------------------------------------+
    | gene_id           | transcript_consequences[0].gene_id         |
    +-------------------+--------------------------------------------+
    | transcript_id     | transcript_consequences[0].transcript_id   |
    +-------------------+--------------------------------------------+

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


def default_vep_cli_json_extractor(response_item: dict) -> dict:
    """An example function to extract VEP data from the CLI.

    When using the VEP CLI, the output is a file where each line is a JSON
    record.  The format of the JSON record is similar to the API response,
    but not identical.  Using the command below with v112 Ensemble VEP,
    the translations are listed in the table below the command:

    .. code-block:: bash

        vep \\
            -i some.vcf \\
            -o some.json.gz \\
            --variant_class \\
            --nearest symbol \\
            --pick \\
            --stats_text \\
            --stats_file some_vep.txt \\
            --compress_output bgzip \\
            --json --assembly GRCh37 \\
            --species homo_sapiens \\
            --cache --cache_version 112 \\
            --dir_cache /data/vep_cache/112/GRCh37

    
    +-------------------+--------------------------------------------+
    | extracted key     | original response_item key                 |
    +===================+============================================+
    | contig            | seq_region_name                            |
    +-------------------+--------------------------------------------+
    | pos               | start                                      |
    +-------------------+--------------------------------------------+
    | ref               | allele_string (first allele, '/' sep)      |
    +-------------------+--------------------------------------------+
    | alt (comma-sep)   | allele_string (all other alleles, '/' sep) |
    +-------------------+--------------------------------------------+
    | type              | variant_class                              |
    +-------------------+--------------------------------------------+
    | effect            | most_severe_consequence                    |
    +-------------------+--------------------------------------------+
    | gene_name         | nearest[0]                                 |
    +-------------------+--------------------------------------------+
    | gene_id           | transcript_consequences[0].gene_id         |
    +-------------------+--------------------------------------------+
    | transcript_id     | transcript_consequences[0].transcript_id   |
    +-------------------+--------------------------------------------+

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
    if "nearest" not in response_item:
        new_item["gene_name"] = None
    else:
        try:
            new_item["gene_name"] = response_item["nearest"][0]
        except IndexError:
            new_item["gene_name"] = None
    if "transcript_consequences" not in response_item:
        new_item["gene_id"] = None
        new_item["transcript_id"] = None
    else:
        new_item["gene_id"] = response_item["transcript_consequences"][0]["gene_id"]
        new_item["transcript_id"] = response_item["transcript_consequences"][0][
            "transcript_id"
        ]
    return new_item


def platypus_vcf_record_extractor(
    record: pysam.VariantRecord,
    **addl_cols: typing.Callable[[pysam.VariantRecord], typing.Union[int, float, str]],
) -> dict:
    """Example function to extract data from a VCF record.

    Currently there actually isn't a higher-level function to call this in an
    analogous way to the VEP API response extractor, but below is the general
    pattern:

    .. code-block:: python

        import pysam
        import pathlib
        import polars as pl

        from varona import extract

        vcf_path = pathlib.Path("/path/to/file.vcf")
        data = []
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                extracted_data = extract.platypus_vcf_record_extractor(record)
                data.append(extracted_data)
        df = pl.DataFrame(data)

    The idea is that the "extractor" will transform the VCF record into a
    dictionary that will in turn be used as a row in a DataFrame.  This function
    is therefore just and example extractor, and the one that the Varona
    command-line tool uses.  If alternative fields from the VCF record are
    needed, the goal is to make substituting this extractor with a different
    one as easy as possible.

    Below are the members of the VCF record object that are extracted by this
    function, along with the key names they are assigned to in the returned
    dictionary.

    +-------------------+--------------------------------------------+
    | extracted key     | VCF record member                          |
    +===================+============================================+
    | contig            | contig                                     |
    +-------------------+--------------------------------------------+
    | pos               | pos                                        |
    +-------------------+--------------------------------------------+
    | ref               | ref                                        |
    +-------------------+--------------------------------------------+
    | alt (comma-sep)   | alts (list)                                |
    +-------------------+--------------------------------------------+
    | sequence_depth    | info["TC"]                                 |
    +-------------------+--------------------------------------------+
    | max_variant_reads | max(info["TR"])                            |
    +-------------------+--------------------------------------------+
    | variant_read_pct  | max_variant_reads / sequence_depth * 100   |
    +-------------------+--------------------------------------------+

    :param record: A :class:`pysam.VariantRecord` object.
    :param addl_cols: Additional columns to extract from the record.
        The keys are the column names and the values are functions that
        take the record and return the value for that column.  The reason
        for having this additional layer of callback is to accommodate
        three competing MAF calculations in a single function.  See
        the :mod:`varona.maf` module for more information on those
        functions, but these are a variable number of keyword arguments,
        where the argument name will become a key in the returned
        dictionary, and the argument value is a function also taking
        a :class:`pysam.VariantRecord` object, but rather than returning
        a dictionary, returning a scalar such as in :func:`varona.maf.maf_from_fr`.

    """
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
