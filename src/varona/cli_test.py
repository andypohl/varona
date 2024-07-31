"""Tests for the CLI module.
"""

import io
import pathlib
import sys
import unittest
from unittest import mock

import varona
from varona import cli, ensembl, maf


class TestVaronaArgsParser(unittest.TestCase):

    @unittest.mock.patch("sys.argv", ["varona", "--help"])
    def test_help_message(self):
        with self.assertRaises(SystemExit):
            cli.VaronaArgumentParser().parse_args()

    @unittest.mock.patch("sys.argv", ["varona", "input.vcf", "output.csv"])
    def test_positional_arguments(self):
        args = cli.VaronaArgumentParser().parse_args()
        self.assertEqual(args.input_vcf, pathlib.Path("input.vcf"))
        self.assertEqual(args.output_csv, pathlib.Path("output.csv"))

    @unittest.mock.patch("sys.argv", ["varona", "input.vcf", "output.csv"])
    def test_log_level_default(self):
        args = cli.VaronaArgumentParser().parse_args()
        self.assertEqual(args.log_level, "warning")

    def test_log_level_patching(self):
        for level in ["debug", "info", "warning", "error"]:
            with mock.patch(
                "sys.argv", ["varona", "--log-level", level, "input.vcf", "output.csv"]
            ):
                parser = cli.VaronaArgumentParser()
                args = parser.parse_args()
                self.assertEqual(args.log_level, level)

    def test_log_level(self):
        for level in ["debug", "info", "warning", "error"]:
            parser = cli.VaronaArgumentParser()
            args = parser.parse_args(["--log-level", level, "input.vcf", "output.csv"])
            self.assertEqual(args.log_level, level)

    def test_supplied_args_fail_2(self):
        parser = cli.VaronaArgumentParser()
        # with self.assertRaises(SystemExit):
        parser.parse_args(["--log-level", "debug", "input.vcf", "output.csv"])

    def test_assembly_default(self):
        args = self.parser.parse_args(["input.vcf", "output.csv"])
        self.assertEqual(args.assembly, ensembl.Assembly.GRCH37)

    def test_assembly_choices(self):
        for assembly in ensembl.Assembly:
            args = self.parser.parse_args(
                ["--assembly", assembly.name, "input.vcf", "output.csv"]
            )
            self.assertEqual(args.assembly, assembly)

    def test_assembly_choices_lower(self):
        for assembly in ensembl.Assembly:
            args = self.parser.parse_args(
                ["--assembly", assembly.name.lower(), "input.vcf", "output.csv"]
            )
            self.assertEqual(args.assembly, assembly)

    def test_maf_default(self):
        args = self.parser.parse_args(["input.vcf", "output.csv"])
        self.assertEqual(args.maf, maf.MafMethod.SAMPLES)

    def test_maf_choices(self):
        for method in maf.MafMethod:
            args = self.parser.parse_args(
                ["--maf", method.name, "input.vcf", "output.csv"]
            )
            self.assertEqual(args.maf, method)

    def test_version(self):
        args = self.parser.parse_args(["--version"])
        self.assertEqual(args.version, True)


class TestMainVersion(unittest.TestCase):

    @mock.patch("sys.argv", ["varona", "--version"])
    def test_version(self, *_):
        with mock.patch("sys.stdout", new=io.StringIO()) as mock_stdout:
            cli.main()
            self.assertEqual(mock_stdout.getvalue().strip(), varona.__version__)
