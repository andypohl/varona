"""Module for argument parsing and CLI functionality.
"""

import argparse
import functools
import logging
import pathlib

import httpx
import polars as pl

import varona
from varona import bcftools, ensembl, extract, maf, varona

logger = logging.getLogger("varona.cli")

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


def varona_args_parser() -> argparse.ArgumentParser:
    """Setup for the varona CLI.

    The CLI ``--help`` message is shown below:

    .. code-block:: text

        usage: varona
            [-h]
            [--log-level {debug,info,warning,error}]
            [--assembly {GRCh37,GRCh38}]
            [--maf {FR,SAMPLES,BCFTOOLS}]
            input_vcf output_csv

        Annotate a VCF file.

        positional arguments:
          input_vcf             Path to the input VCF file
          output_csv            Path to the output CSV file

        options:
          -h, --help            show this help message and exit
          --log-level {debug,info,warning,error}
                                Set the logging level (default: WARNING)
          --assembly {GRCh37,GRCh38}
                                genome assembly used in Ensembl VEP API (default: GRCh37)
          --maf {FR,SAMPLES,BCFTOOLS}
                                MAF calculation method (default: SAMPLES)

    :return: the configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Annotate a VCF file.", prog="varona")
    # Positional arguments
    parser.add_argument(
        "input_vcf", type=pathlib.Path, help="Path to the input VCF file"
    )
    parser.add_argument(
        "output_csv", type=pathlib.Path, help="Path to the output CSV file"
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["debug", "info", "warning", "error"],
        default="warning",
        help="Set the logging level (default: WARNING)",
    )
    # Optional arguments using Enum classes
    parser.add_argument(
        "--assembly",
        type=ensembl.Assembly,
        choices=list(ensembl.Assembly),
        default=ensembl.Assembly.GRCH37,
        help="genome assembly used in Ensembl VEP API (default: GRCh37)",
    )
    parser.add_argument(
        "--maf",
        type=maf.MafMethod,
        choices=list(maf.MafMethod),
        default=maf.MafMethod.SAMPLES,
        help="MAF calculation method (default: SAMPLES)",
    )
    return parser


def platypus_dataframe(
    vcf_path: pathlib.Path,
    maf_method: maf.MafMethod = maf.MafMethod.SAMPLES,
    timeout: int = 300,
    genome_assembly: ensembl.Assembly = ensembl.Assembly.GRCH37,
    vcf_extractor=extract.platypus_vcf_record_extractor,
    api_extractor=extract.default_vep_response_extractor,
) -> pl.DataFrame:
    """Read a Platypus VCF file into a DataFrame.

    :param vcf_path: The path to the Platypus VCF file.
    :param maf_method: The method to use for calculating the MAF.
    :param timeout: The timeout (seconds) for the VEP API query.
    :param genome_assembly: The genome assembly used in the Ensembl VEP API.
    :param vcf_extractor: The function to extract data from the VCF.
    :param api_extractor: The function to extract data from the VEP API response.
    :return: A DataFrame with the VCF data.
    """
    # VCF part
    lst = []
    maf_func = functools.partial(maf.maf_from_method, method=maf_method)
    if maf_method == maf.MafMethod.BCFTOOLS:
        with bcftools.VariantFileFilledInTags(vcf_path, fillin_tags=["MAF"]) as vf:
            for record in vf:
                new_item = vcf_extractor(record, maf=maf_func)
                lst.append(new_item)
    else:
        lst = list(varona._vcf_rows(vcf_path, vcf_extractor))
    vcf_df = pl.DataFrame(lst, schema=VCF_DF_SCHEMA)
    # API part
    chunks = list(ensembl.vcf_to_vep_query_data(vcf_path))
    n_chunks = len(chunks)
    api_df = pl.DataFrame({k: [] for k in API_DF_SCHEMA.keys()}, schema=API_DF_SCHEMA)
    with httpx.Client(
        limits=httpx.Limits(max_connections=5, max_keepalive_connections=5),
        timeout=httpx.Timeout(float(timeout)),
    ) as client:
        for ix, chunk in enumerate(chunks, start=1):
            chunk_df = varona.vep_api_dataframe(
                client, chunk, genome_assembly, api_extractor, schema=API_DF_SCHEMA
            )
            api_df = api_df.vstack(chunk_df)
            logger.info(f"processed {ix}/{n_chunks} chunks from VEP API")
    api_df.rechunk()
    combined_df = vcf_df.join(api_df, on=["contig", "pos", "ref", "alt"], how="left")
    return combined_df


def main():
    args = varona_args_parser().parse_args()
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", varona.__version__)
    df = platypus_dataframe(args.input_vcf, maf_method=args.maf)
    logger.info("writing CSV file: %s", str(args.output_csv))
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
