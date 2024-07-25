"""Tests for the enum module.
"""

import enum as py_enum
import importlib
import sys
import unittest
from unittest import mock

from varona import enum as varona_enum


class TestBackportStrEnum(unittest.TestCase):

    @unittest.skipIf(sys.version_info < (3, 11), "requires Python 3.11 or higher")
    def test_backport_equiv(self):
        """Test that the backport is equivalent (for our use) as Python's StrEnum."""

        class TestEnumA(varona_enum.BackportStrEnum):
            A = "A"
            B = "B"
            C = "C"

        class TestEnumB(py_enum.StrEnum):
            A = "A"
            B = "B"
            C = "C"

        self.assertEqual(TestEnumA.A, TestEnumB.A)
        self.assertListEqual(list(TestEnumA), list(TestEnumB))
        self.assertEqual(TestEnumA.A, TestEnumB.A)
        self.assertEqual(TestEnumA.A.value, TestEnumB.A.value)
        self.assertEqual(TestEnumA.A.name, TestEnumB.A.name)
        self.assertRaises(ValueError, lambda: TestEnumA("a"))
        self.assertRaises(ValueError, lambda: TestEnumB("a"))

    @unittest.skipIf(sys.version_info < (3, 11), "requires Python 3.11 or higher")
    def test_backport(self):
        """Mock deleting the StrEnum class and test that the backport works."""
        with mock.patch.object(py_enum, "StrEnum", create=False):
            delattr(py_enum, "StrEnum")
            importlib.reload(varona_enum)

    def test_backport_multenum(self):
        """Really not likely just check multivalue enums."""

        class TestMultiEnum(varona_enum.BackportStrEnum):
            ORIGIN = (0, 0)
            POINT_A = (1, 2)
            POINT_B = (3, 4)

        self.assertEqual(len(TestMultiEnum), 3)


class TestCiStrEnum(unittest.TestCase):

    def test_ci_str_enum(self):
        """Test the case-insensitive string enum."""

        class TestEnum(varona_enum.CiStrEnum):
            A = "A"
            B = "B"
            C = "C"

        self.assertEqual(TestEnum("A"), TestEnum("a"))
        self.assertEqual(TestEnum("A").value, TestEnum("a").value)
        self.assertEqual(TestEnum("A").value, "A")
        self.assertEqual(TestEnum("a").name, "A")
        self.assertEqual(TestEnum("A").name, str(TestEnum("A")))
        self.assertEqual(TestEnum("A").name, str(TestEnum("a")))
        self.assertRaises(ValueError, lambda: TestEnum("D"))
