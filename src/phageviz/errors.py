"""Exception types. Every message is written for the end user and is safe to show in the UI."""

from __future__ import annotations


class PhageVizError(Exception):
    """Base class for all expected errors. ``str(exc)`` is a user-facing message."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


class InvalidRecord(PhageVizError):
    """A genome record or its features table breaks the data model (SPEC section 4)."""
