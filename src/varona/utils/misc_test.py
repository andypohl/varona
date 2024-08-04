"""Test the misc module.
"""

import pathlib
import unittest

from varona.utils import misc


class TestStems(unittest.TestCase):

    def test_one_suffix(self):
        file_path = pathlib.PurePath("file.txt")
        self.assertEqual(misc.multi_suffix_stem(file_path), "file")

    def test_multi_suffix(self):
        file_path = pathlib.PurePath("file.vcf.gz")
        self.assertEqual(misc.multi_suffix_stem(file_path), "file")
        file_path = pathlib.PurePath("file.vcf.gz.gz")
        self.assertEqual(misc.multi_suffix_stem(file_path), "file")
