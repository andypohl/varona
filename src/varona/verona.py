"""High-level routines for the Varona library.

All of the functions and classes in this module are imported into the
top-level library namespace.
"""

import logging
import pathlib

import httpx
import polars as pl
import pysam
from varona import ensembl, extract, maf

logger = logging.getLogger("varona.varona")

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


def platypus_vcf_dataframe(
    vcf_path: pathlib.Path, maf_method: maf.MafMethod = maf.MafMethod.SAMPLES
) -> pl.DataFrame:
    """Read a Platypus VCF file into a DataFrame.

    :param vcf_path: The path to the Platypus VCF file.
    :param maf_method: The method to use for calculating the MAF.
    :return: A DataFrame with the VCF data.
    """
    lst = []
    VarFile = pysam.VariantFile
    with VarFile(vcf_path, "r") as vf:
        for record in vf:
            new_item = extract.platypus_vcf_record_extractor(
                record, maf=maf.maf_from_samples
            )
            lst.append(new_item)
    vcf_df = pl.DataFrame(lst, schema=VCF_DF_SCHEMA)
    api_df = pl.DataFrame(schema=API_DF_SCHEMA)
    limits = httpx.Limits(max_connections=5, max_keepalive_connections=5)
    timeout = httpx.Timeout(300.0)
    chunks = list(ensembl.vcf_to_vep_query_data(vcf_path))[:10]
    n_chunks = len(chunks)
    with httpx.Client(limits=limits, timeout=timeout) as client:
        for ix, chunk in enumerate(chunks, start=1):
            data = ensembl.query_vep_api(
                client, chunk, response_extractor=extract.default_vep_response_extractor
            )
            api_df = api_df.vstack(pl.DataFrame(data, schema=API_DF_SCHEMA))
            logger.info(f"processed {ix}/{n_chunks} chunks from VEP API")
    api_df.rechunk()
    combined_df = vcf_df.join(api_df, on=["contig", "pos", "ref", "alt"], how="left")
    return combined_df
