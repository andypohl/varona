import argparse
import logging
import pathlib

from varona.utils import split


def vcf_split_args() -> argparse.ArgumentParser:
    """Set up the argument parser for the vcf_split CLI."""
    parser = argparse.ArgumentParser(
        description="Split a VCF file into smaller pieces.", prog="vcf_split"
    )
    parser.add_argument("in_vcf", type=pathlib.Path, help="The VCF file to split.")
    parser.add_argument(
        "--log-level",
        type=str,
        choices=["debug", "info", "warning", "error"],
        default="warning",
        help="Set the logging level (default: WARNING)",
    )
    parser.add_argument(
        "--out-dir",
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help="The directory to save the split VCF files.",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        default=True,
        help="Whether to compress the output files.",
    )
    split_group = parser.add_mutually_exclusive_group(required=True)
    split_group.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="The number of records per chunk.",
    )
    split_group.add_argument(
        "--n-chunks",
        type=int,
        default=None,
        help="The number of chunks to split the VCF into.",
    )
    return parser


def cli_main():
    """Entrypoint for the vcf_split CLI."""
    args = vcf_split_args().parse_args()
    logging.basicConfig(level=args.log_level.upper())
    split.split_vcf(
        args.in_vcf,
        args.out_dir,
        chunk_size=args.chunk_size,
        n_chunks=args.n_chunks,
        compress=args.compress,
    )


if __name__ == "__main__":
    cli_main()
