"""Splits up a VCF file specifically for Varona+Nextflow.
"""

import logging
import math
import pathlib

import pysam

from varona.utils import misc

logger = logging.getLogger("varona.utils.split")


def split_vcf(
    vcf_path: pathlib.Path,
    out_dir: pathlib.Path,
    chunk_size: int | None = None,
    n_chunks: int | None = None,
) -> list[pathlib.Path]:
    """Split a VCF file into smaller chunks.

    :param vcf_path: The path to the VCF file.
    :param out_dir: The directory to save the chunks.
    :param chunk_size: The number of records per chunk.
    :param n_chunks: The number of chunks to split the VCF into.  If set,
        then ``chunk_size`` is ignored.
    :return: The list of paths to the split VCF files.
    """
    if not chunk_size and not n_chunks:
        raise ValueError("Either chunk_size or n_chunks must be set")
    # quick pass to get the number of records
    splits = []
    with pysam.VariantFile(vcf_path) as vcf:
        n_records = sum(1 for _ in vcf)
    if n_chunks:
        chunk_size = int(math.ceil(n_records / n_chunks))
    else:
        n_chunks = int(math.ceil(n_records / chunk_size))
    logger.info(f"Splitting {vcf_path} into {n_chunks} pieces")
    n_zeros = int(math.ceil(math.log10(n_chunks))) + 1
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = misc.multi_suffix_stem(vcf_path)
    with pysam.VariantFile(vcf_path) as in_vcf:
        records = in_vcf.fetch()
        for file_ix in range(1, n_chunks + 1):
            out_path = out_dir / f"{stem}_{str(file_ix).zfill(n_zeros)}.vcf"
            with pysam.VariantFile(
                out_path,
                "w",
                header=in_vcf.header,
            ) as out_vcf:
                for _ in range(chunk_size):
                    try:
                        record = next(records)
                        out_vcf.write(record)
                    except StopIteration:
                        break
            gz_path = out_path.with_suffix(".vcf.gz")
            pysam.tabix_compress(out_path, gz_path)
            splits.append(gz_path)
            out_path.unlink()
            logger.info(f"Wrote {gz_path}")
    logger.info(f"Finished splitting {vcf_path}")
    return splits
