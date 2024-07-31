"""Module for argument parsing and CLI functionality.
"""

import argparse
import logging
import os
import pathlib

import varona
from varona import ensembl, maf, platypus

logger = logging.getLogger("varona.cli")


def varona_args_parser() -> argparse.ArgumentParser:
    """Setup for the varona CLI.

    The CLI ``--help`` message is shown below:

    .. code-block:: text

        usage: varona [-h]
            [--log-level {debug,info,warning,error}]
            [--assembly {GRCH37,GRCH38}]
            [--maf {FR,SAMPLES,BCFTOOLS}]
            [--no-vep]
            [--vep-data VEP_DATA]
            input_vcf output_csv

        Annotate a VCF file.

        positional arguments:
          input_vcf             Path to the input VCF file
          output_csv            Path to the output CSV file

        options:
          -h, --help            show this help message and exit
          --log-level {debug,info,warning,error}
                                Set the logging level (default: WARNING)
          --assembly {GRCH37,GRCH38}
                                genome assembly used in Ensembl VEP API (default: GRCh37)
          --maf {FR,SAMPLES,BCFTOOLS}
                                MAF calculation method (default: SAMPLES)
          --no-vep              Skip VEP API querying (no effect if --vep-data is provided)
          --vep-data VEP_DATA   Path to VEP output file (currently unused)
          --version             Show program's version number and exit

    :return: the configured argument parser.
    """
    parser = argparse.ArgumentParser(description="Annotate a VCF file.", prog="varona")
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
    parser.add_argument(
        "--no-vep",
        action="store_true",
        help="Skip VEP API querying (no effect if --vep-data is provided)",
    )
    parser.add_argument(
        "--vep-data",
        type=pathlib.Path,
        help="Path to VEP output file (currently unused)",
        default=None,
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show program's version number and exit",
    )
    initial_args, _ = parser.parse_known_args()
    if initial_args.version:
        return parser
    # Positional arguments added if --version is not used
    parser.add_argument(
        "input_vcf",
        type=pathlib.Path,
        help="Path to the input VCF file",
    )
    parser.add_argument(
        "output_csv", type=pathlib.Path, help="Path to the output CSV file"
    )
    return parser


def main():
    parser = varona_args_parser()
    # Extra gymnasics to handle the case where --version is used
    # seems a bit delicate
    if parser.parse_known_args()[0].version:
        print(f"{varona.__version__}")
        return
    args = parser.parse_args()
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", varona.__version__)
    df = platypus.platypus_dataframe(
        args.input_vcf, maf_method=args.maf, no_vep=args.no_vep
    )
    logger.info("writing CSV file: %s", str(args.output_csv))
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
