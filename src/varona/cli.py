"""Module for argument parsing and CLI functionality.
"""

import argparse
import logging
import pathlib

import varona
from varona import ensembl, maf

logger = logging.getLogger("varona.cli")


def main():
    # Create the argument parser
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
        default=ensembl.Assembly.GRCh37,
        help="genome assembly used in Ensembl VEP API (default: GRCh37)",
    )
    parser.add_argument(
        "--maf",
        type=maf.MafMethod,
        choices=list(maf.MafMethod),
        default=maf.MafMethod.SAMPLES,
        help="MAF calculation method (default: SAMPLES)",
    )
    args = parser.parse_args()
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", varona.__version__)
    df = varona.varona_dataframe(args.input_vcf, maf_method=args.maf)
    logger.info("writing CSV file: %s", str(args.output_csv))
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
