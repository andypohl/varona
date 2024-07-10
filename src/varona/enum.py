"""Helper class for case-insensitive string enums.
"""

import enum
from enum import auto


class CiStrEnum(enum.StrEnum):
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
