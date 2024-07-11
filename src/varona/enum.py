"""Helper class for case-insensitive string enums.
"""

import enum
import sys
from enum import auto


class BackportStrEnum(str, enum.Enum):
    """
    Enum where members are also instances of str.

    This is available in Python 3.11+ so it's a backport for older Python versions.
    """

    def __new__(cls, *values):
        if len(values) == 1:
            value = values[0]
        else:
            value = values
        obj = str.__new__(cls, value)
        obj._value_ = value
        return obj


try:
    StrEnum = enum.StrEnum
except AttributeError:
    StrEnum = BackportStrEnum


class CiStrEnum(StrEnum):
    """Gets the enum member by case-insensitive string value.

    From Python `docs <https://docs.python.org/3/library/enum.html#enum.Enum._missing_>`_.
    """

    @classmethod
    def _missing_(cls, value):
        value = value.lower()
        for member in cls:
            if member.value == value:
                return member
        return None

    def __str__(self):
        """This helps list(cls) return the string values (not the enum members)"""
        return self.name
