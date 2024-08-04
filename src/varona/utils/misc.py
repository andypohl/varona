"""Miscellaneous utility functions."""

import pathlib


def multi_suffix_stem(file_path: pathlib.Path) -> str:
    """For paths like ``/path/to/file.vcf.gz``, return ``file``.

    The :attr:`pathlib.Path.stem` attribute is annoying that it
    only has the last suffix.  This function removes all
    suffixes.

    This function is a bit aggressive, removing everything
    after the first ".".  Be mindful of this when formulating
    the file naming strategy.

    Examples:

    >>> from pathlib import Path
    >>> multi_suffix_stem(Path("file.txt"))
    'file'
    >>> multi_suffix_stem(Path("file.vcf.gz"))
    'file'

    :param file_path: The file path.
    :return: The stem with the multi suffix removed.
    """
    return file_path.name.removesuffix("".join(file_path.suffixes))
