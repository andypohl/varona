"""Module for argument parsing and CLI functionality.
"""

import argparse
import logging
import pathlib

import verona
from verona import ensembl, maf

logger = logging.getLogger("verona.cli")


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
        help="genome assembly used in VEP API (default: GRCh37)",
    )
    parser.add_argument(
        "--maf",
        type=maf.MafMethod,
        choices=list(maf.MafMethod),
        default=maf.MafMethod.SAMPLES,
        help="MAF calculation method (default: SAMPLES)",
    )
    # Parse the arguments
    args = parser.parse_args()
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", verona.__version__)
    # Print the parsed arguments
    df = verona.platypus_vcf_dataframe(args.input_vcf)
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
