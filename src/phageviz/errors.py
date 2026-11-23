"""Exception types. Every message is written for the end user and is safe to show in the UI."""

from __future__ import annotations


class PhageVizError(Exception):
    """Base class for all expected errors. ``str(exc)`` is a user-facing message."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


class InvalidRecord(PhageVizError):
    """A genome record or its features table breaks the data model (SPEC section 4)."""


class UnsupportedFormat(PhageVizError):
    """The input isn't a format (or kind of file) PhageViz can load."""


class ParseError(PhageVizError):
    """The input looks like a supported format but can't be read."""


class MissingSequence(PhageVizError):
    """A record has annotation but no DNA sequence (e.g. a GenBank CON record)."""


class RecordTooLarge(PhageVizError):
    """A genome is longer than the supported maximum (SPEC 6.5)."""


class RecordTooShort(PhageVizError):
    """A genome is shorter than the supported minimum (SPEC 6.5)."""
