"""Module for argument parsing and CLI functionality.
"""

import argparse
import logging
import pathlib

import varona
from varona import dataframe, ensembl, maf, platypus

logger = logging.getLogger("varona.cli")


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


def main():
    args = varona_args_parser().parse_args()
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", varona.__version__)
    df = platypus.platypus_dataframe(args.input_vcf, maf_method=args.maf)
    logger.info("writing CSV file: %s", str(args.output_csv))
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
