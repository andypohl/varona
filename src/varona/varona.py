"""High-level routines for the Varona library.

All of the functions and classes in this module are imported into the
top-level library namespace.
"""

import functools
import logging
import pathlib
import typing

import httpx
import polars as pl
import pysam

from varona import bcftools, ensembl, extract, maf

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


def _vcf_rows(
    vcf_path: pathlib.Path, vcf_extractor: typing.Callable[[pysam.VariantRecord], dict]
):
    """Helper function to extract rows from a VCF file.

    :param vcf_path: The path to the VCF file.
    :param vcf_extractor: The function to extract data from the VCF.
    :yields: A dictionary of extracted data from the VCF.
    """
    with pysam.VariantFile(vcf_path, "r") as vf:
        for record in vf:
            new_item = vcf_extractor(record)
            yield new_item


def vcf_dataframe(
    vcf_path: pathlib.Path,
    vcf_extractor: typing.Callable[[pysam.VariantRecord], dict],
    schema: dict[str, typing.Any] | None = None,
) -> pl.DataFrame:
    """From the records in a VCF file, make a dataframe given an extractor.

    .. code-block:: python

        import pathlib
        import polars as pl
        import pysam
        from varona import vcf_dataframe

        def example_extractor(record: pysam.VariantRecord) -> dict:
            return {
                "contig": record.contig,
                "pos": record.pos,
                "ref": record.ref,
                "alt": record.alts[0]
            }

        # Make a DataFrame from the VCF file.  The columns laid out by
        # the extractor function.

        vcf_path = pathlib.Path("/path/to/file.vcf")
        df = vcf_dataframe(vcf_path, example_extractor)
        print(df)
        ##shape: (5, 4)
        ##┌────────┬─────────┬─────┬─────┐
        ##│ contig ┆ pos     ┆ ref ┆ alt │
        ##│ ---    ┆ ---     ┆ --- ┆ --- │
        ##│ str    ┆ i64     ┆ str ┆ str │
        ##╞════════╪═════════╪═════╪═════╡
        ##│ 1      ┆ 1158631 ┆ A   ┆ G   │
        ##│ 1      ┆ 1246004 ┆ A   ┆ G   │
        ##│ 1      ┆ 1249187 ┆ G   ┆ A   │
        ##│ 1      ┆ 1261824 ┆ G   ┆ C   │
        ##│ 1      ┆ 1387667 ┆ C   ┆ G   │
        ##└────────┴─────────┴─────┴─────┘

    :param vcf_path: The path to the VCF file.
    :param vcf_extractor: The function to extract data from the VCF.
    :param schema: Optional schema for the DataFrame to help enforce column types.
    :return: DataFrame with the extracted data.
    """
    return pl.LazyFrame(_vcf_rows(vcf_path, vcf_extractor), schema=schema).collect()


def vep_api_dataframe(
    client: httpx.Client,
    loci_list: list[str],
    genome_assembly: ensembl.Assembly,
    api_extractor: typing.Callable[[dict], dict],
    schema: dict[str, typing.Any] | None = None,
) -> pl.DataFrame:
    """Query the Ensembl VEP API and make a DataFrame using a provided extractor.

    Like :func:`vcf_dataframe`, this is a vehicle for a custom extractor function
    to be used on the response dictionaries from the Ensembl VEP API.  Below is an
    example of how to use this function.  A :class:`httpx.Client` still needs to
    be supplied.

    .. code-block:: python

            import pathlib
            import polars as pl
            import httpx
            from varona import vep_api_dataframe, ensembl

            def example_extractor(response: dict) -> dict:
                return {
                    "contig": response["seq_region_name"],
                    "pos": response["start"],
                    "type": response["variant_class"]
                }

            loci_list = [
                "1 1158631 . A G . . .",
                "1 91859795 . TATGTGA CATGTGA,CATGTGG . . .",
            ]
            with httpx.Client(
                limits=httpx.Limits(
                    max_connections=5,
                    max_keepalive_connections=5
                ),
                timeout=httpx.Timeout(float(300)),
            ) as client:
                api_df = vep_api_dataframe(
                    client,
                    loci_list,
                    ensembl.Assembly.GRCh37,
                    example_extractor
                )
                print(api_df)
            ##shape: (2, 3)
            ##┌────────┬──────────┬──────────────┐
            ##│ contig   ┆ pos        ┆ type            │
            ##│ ---      ┆ ---        ┆ ---             │
            ##│ str      ┆ i64        ┆ str             │
            ##╞════════╪══════════╪══════════════╡
            ##│ 1        ┆ 1158631    ┆ SNV             │
            ##│ 1        ┆ 91859795   ┆ substitution    │
            ##└────────┴──────────┴──────────────┘

    :param client: The HTTPX client to use for the API query.
    :param loci_list: The list of loci to query the API.
    :param genome_assembly: The genome assembly used in the Ensembl VEP API.
    :param api_extractor: The function to extract data from the VEP API response.
    :param schema: Optional schema for the DataFrame to help enforce column types.
    :return: A DataFrame with the data from the VEP API.
    """
    data = ensembl.query_vep_api(
        client, loci_list, genome_assembly, response_extractor=api_extractor
    )
    return pl.DataFrame(data, schema=schema)


def varona_dataframe(
    vcf_path: pathlib.Path,
    maf_method: maf.MafMethod = maf.MafMethod.SAMPLES,
    timeout: int = 300,
    genome_assembly: ensembl.Assembly = ensembl.Assembly.GRCh37,
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
        lst = list(_vcf_rows(vcf_path, vcf_extractor))
    vcf_df = pl.DataFrame(lst, schema=VCF_DF_SCHEMA)
    # API part
    chunks = list(ensembl.vcf_to_vep_query_data(vcf_path))
    n_chunks = len(chunks)
    with httpx.Client(
        limits=httpx.Limits(max_connections=5, max_keepalive_connections=5),
        timeout=httpx.Timeout(float(timeout)),
    ) as client:
        for ix, chunk in enumerate(chunks, start=1):
            chunk_df = vep_api_dataframe(
                client, chunk, genome_assembly, api_extractor, schema=API_DF_SCHEMA
            )
            api_df = api_df.vstack(chunk_df)
            logger.info(f"processed {ix}/{n_chunks} chunks from VEP API")
    api_df.rechunk()
    combined_df = vcf_df.join(api_df, on=["contig", "pos", "ref", "alt"], how="left")
    return combined_df
