"""Module for argument parsing and CLI functionality.
"""

import argparse
import logging
import os
import pathlib

import varona
from varona import ensembl, maf, platypus

logger = logging.getLogger("varona.cli")


class VaronaArgumentParser(argparse.ArgumentParser):
    """Varona argument parser.

    Subclass of :class:`argparse.ArgumentParser` that essentially allows the
    required positional arguments to be optional in the presence of the
    ``--version`` flag.

    .. warning::

        Due to the way this behavior is implemented, the
        :meth:`parse_args` method can only be called once, and calling
        twice will raise a :class:`RuntimeError`.

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

    """

    def __init__(self):
        super(VaronaArgumentParser, self).__init__(
            description="Annotate a VCF file.", prog="varona"
        )
        self._add_optional_args()
        # hacky flag to make sure .parse_args() is called only once
        self.args_parsed = False

    def _add_optional_args(self):
        """Add optional arguments to the parser."""
        self.add_argument(
            "--log-level",
            type=str,
            choices=["debug", "info", "warning", "error"],
            default="warning",
            help="Set the logging level (default: WARNING)",
        )
        # Optional arguments using Enum classes
        self.add_argument(
            "--assembly",
            type=ensembl.Assembly,
            choices=list(ensembl.Assembly),
            default=ensembl.Assembly.GRCH37,
            help="genome assembly used in Ensembl VEP API (default: GRCh37)",
        )
        self.add_argument(
            "--maf",
            type=maf.MafMethod,
            choices=list(maf.MafMethod),
            default=maf.MafMethod.SAMPLES,
            help="MAF calculation method (default: SAMPLES)",
        )
        self.add_argument(
            "--no-vep",
            action="store_true",
            help="Skip VEP API querying (no effect if --vep-data is provided)",
        )
        self.add_argument(
            "--vep-data",
            type=pathlib.Path,
            help="Path to VEP output file (currently unused)",
            default=None,
        )
        self.add_argument(
            "--version",
            action="store_true",
            help="Show program's version number and exit",
        )
        initial_args, _ = self.parse_known_args()

    def _add_positional_args(self):
        # Positional arguments added if --version is not used
        self.add_argument(
            "input_vcf",
            type=pathlib.Path,
            help="Path to the input VCF file",
        )
        self.add_argument(
            "output_csv", type=pathlib.Path, help="Path to the output CSV file"
        )

    def parse_args(self, args=None, namespace=None):
        """Parse command-line arguments.

        If the ``--version`` flag is not used, then the positional
        arguments are added to the parser.

        :param args: List of strings to parse. Default is taken from sys.argv.
        :param namespace: An object to take the attributes. Default is a new empty namespace.
        :return: The parsed arguments.
        """
        # Add positional arguments if --version is not used
        if self.args_parsed:
            raise RuntimeError("parse_args() called more than once")
        initial_args, _ = self.parse_known_args(args=args, namespace=namespace)
        if initial_args.version:
            return initial_args
        self._add_positional_args()
        # something here isn't quite right, the
        parsed = super().parse_args(args=args, namespace=namespace)
        self.args_parsed = True
        return parsed


def main():
    parser = VaronaArgumentParser()
    args = parser.parse_args()
    if args.version:
        print(f"{varona.__version__}")
        return
    # Set the logging level
    logging.basicConfig(level=args.log_level.upper())
    logger.info("varona version: %s", varona.__version__)
    df = platypus.platypus_dataframe(
        args.input_vcf,
        maf_method=args.maf,
        no_vep=args.no_vep,
        vep_data=args.vep_data,
    )
    logger.info("writing CSV file: %s", str(args.output_csv))
    df.write_csv(args.output_csv)


if __name__ == "__main__":
    main()
