"""High-level module for building a DataFrame from Platypus-style VCF.

This module calls functions from the :mod:`varona.dataframe` module to build
a DataFrame from a Platypus-style VCF file.
"""

import functools
import logging
import pathlib

import httpx
import polars as pl

from varona import bcftools, dataframe, ensembl, extract, maf

logger = logging.getLogger("varona.platypus")

API_DF_SCHEMA = {
    "contig": pl.Utf8,
    "pos": pl.UInt32,
    "ref": pl.Utf8,
    "alt": pl.Utf8,
    "type": pl.Utf8,
    "effect": pl.Utf8,
    "gene_name": pl.Utf8,
    "gene_id": pl.Utf8,
    "transcript_id": pl.Utf8,
}
"""Polars schema for the API DataFrame."""

VCF_DF_SCHEMA = {
    "contig": pl.Utf8,
    "pos": pl.UInt32,
    "ref": pl.Utf8,
    "alt": pl.Utf8,
    "sequence_depth": pl.UInt64,
    "max_variant_reads": pl.UInt64,
    "variant_read_pct": pl.Float64,
    "maf": pl.Float64,
}
"""Polars schema for the VCF DataFrame."""


def platypus_dataframe(
    vcf_path: pathlib.Path,
    maf_method: maf.MafMethod = maf.MafMethod.SAMPLES,
    timeout: int = 300,
    genome_assembly: ensembl.Assembly = ensembl.Assembly.GRCH37,
    vcf_extractor=extract.platypus_vcf_record_extractor,
    api_extractor=extract.default_vep_response_extractor,
    no_vep: bool = False,
    vep_json_path: pathlib.Path | None = None,
) -> pl.DataFrame:
    """Read a Platypus VCF file into a DataFrame.

    :param vcf_path: The path to the Platypus VCF file.
    :param maf_method: The method to use for calculating the MAF.
    :param timeout: The timeout (seconds) for the VEP API query.
    :param genome_assembly: The genome assembly used in the Ensembl VEP API.
    :param vcf_extractor: The function to extract data from the VCF.
    :param api_extractor: The function to extract data from the VEP API response.
    :param no_vep: Skip querying the VEP API.
    :param vep_json_path: Path to the VEP output file from running VEP locally,
        (bypasses querying API).
    :return: A DataFrame with the VCF data.
    """
    # VCF part
    lst = []
    maf_func = functools.partial(maf.maf_from_method, method=maf_method)
    if bcftools.HAVE_BCFTOOLS and maf_method == maf.MafMethod.BCFTOOLS:
        with bcftools.VariantFileFilledInTags(vcf_path, fillin_tags=["MAF"]) as vf:
            for record in vf:
                new_item = vcf_extractor(record, maf=maf_func)
                lst.append(new_item)
    else:
        lst = list(dataframe._vcf_rows(vcf_path, vcf_extractor))
    vcf_df = pl.DataFrame(lst, schema=VCF_DF_SCHEMA)
    if no_vep and vep_json_path is None:
        # no_vep ignored if vep_json_path is provided
        return vcf_df
    vep_df = None
    if vep_json_path:
        vep_df = ensembl.import_vep_data(
            vep_json_path, extractor=extract.default_vep_cli_json_extractor
        )
    else:
        # API part
        chunks = list(ensembl.vcf_to_vep_query_data(vcf_path))
        n_chunks = len(chunks)
        vep_df = pl.DataFrame(
            {k: [] for k in API_DF_SCHEMA.keys()}, schema=API_DF_SCHEMA
        )
        with httpx.Client(
            limits=httpx.Limits(max_connections=5, max_keepalive_connections=5),
            timeout=httpx.Timeout(float(timeout)),
        ) as client:
            for ix, chunk in enumerate(chunks, start=1):
                chunk_df = dataframe.vep_api_dataframe(
                    client, chunk, genome_assembly, api_extractor, schema=API_DF_SCHEMA
                )
                vep_df = vep_df.vstack(chunk_df)
                logger.info(f"processed {ix}/{n_chunks} chunks from VEP API")
        vep_df.rechunk()
    combined_df = vcf_df.join(vep_df, on=["contig", "pos", "ref", "alt"], how="left")
    return combined_df
