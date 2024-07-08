"""Module specifically for using (non-pysam) bcftools.
"""

import pathlib
import shutil
import subprocess
import tempfile

import pysam

HAVE_BCFTOOLS = shutil.which("bcftools") is not None
"""True if the bcftools command is available on the system."""

ALLOWED_TAGS = {"AN", "AC", "AF", "MAF"}
"""Set of allowed tags that can be added to the VCF file."""


def run_bcftools_fill_tags(
    input_vcf: pathlib.Path, output_vcf: pathlib.Path, tags: list[str]
) -> None:
    """
    Run the bcftools fill-tags plugin to add AN, AC, AF, and MAF tags.

    :param input_vcf: Path to the input VCF file (compressed with bgzip).
    :param output_vcf: Path to the output VCF file (compressed with bgzip).
    :param tags: List of tags to add to the VCF file can be any of AN, AC, AF, MAF.
    :raises subprocess.CalledProcessError: If the bcftools command fails.
    """
    if set(tags) - ALLOWED_TAGS:
        raise ValueError("Only AN, AC, AF, and MAF tags are allowed.")
    tag_str = ",".join(tags)
    # Construct the bcftools command
    cmd = [
        "bcftools",
        "plugin",
        "fill-tags",
        str(input_vcf),
        "-o",
        str(output_vcf),
        "-O",
        "z",
        "--",
        "-t",
        tag_str,
    ]
    # Run the command using subprocess
    subprocess.run(cmd, check=True)


class VariantFileFilledInTags(pysam.VariantFile):

    def __init__(self, filename, mode, fillin_tags: list[str], *args, **kwargs):
        if mode != "r":
            raise NotImplementedError("Only VCF read mode is supported (currently).")
        if not HAVE_BCFTOOLS:
            raise RuntimeError(
                "Non-pysam bcftools is needed for preprocess operations."
            )
        self.tmp_dir = pathlib.Path(tempfile.mkdtemp())
        unprocessed = self.tmp_dir / "unprocessed.vcf"
        shutil.copy(filename, unprocessed)
        compressed = self.tmp_dir / "compressed.vcf.gz"
        pysam.tabix_compress(str(unprocessed), str(compressed))
        pysam.tabix_index(str(compressed), preset="vcf")
        processed = self.tmp_dir / "processed.vcf.gz"
        run_bcftools_fill_tags(compressed, processed, fillin_tags)
        super().__init__(processed, mode, *args, **kwargs)

    def close(self):
        super().close()
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
