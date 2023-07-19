"""
This module defines the card classes
"""
from abc import ABC, abstractmethod
import logging
import re
from enum import Enum

from vesta_espresso.utils import parse_pwvals


class InputCard(ABC):
    @property
    @abstractmethod
    def name(self, value):
        pass

    @property
    @abstractmethod
    def opts(self, value):
        pass

    @property
    @abstractmethod
    def default_option(self, value):
        pass

    @property
    @abstractmethod
    def default_deprecated(self, value):
        pass

    @classmethod
    def get_option(cls, option):
        """Initializes a card's options"""
        if option is not None:
            return cls.opts.from_string(option)
        if cls.default_deprecated:
            logging.warning(
                f"No option specified for {cls.name} card. This is deprecated, but {cls.default_option} will be used by default."
            )
        return cls.default_option

    @classmethod
    def split_card_string(cls, s: str):
        """
        Splits a card into an option and a list of values of the correct type.
        :param s: String containing a card (as it would appear in a PWin file)
        :return: option: string for the card's option or None
                 values: list of lists of values for the card

        Example:
        >>> s = "ATOMIC_SPECIES\nH 1.00794 H.UPF\nO 15.9994 O.UPF"
        >>> option, values = InputCard.split_card_string_string(s)
        >>> option, values
        >>> (None, [["H", 1.00794, "H.UPF"], ["O", 15.9994, "O.UPF"]])
        """
        header = s.strip().split("\n")[0]
        body = s.strip().split("\n")[1:]
        if len(header.split()) > 1:
            option = re.sub(r"[()]", "", header.split()[1])
            option = option.lower()
            option = re.sub(r"[()]", "", option)
            option = re.sub(r"[{}]", "", option)
        else:
            option = None
        return cls.get_option(option), parse_pwvals(body)


class CardOptions(Enum):
    """Enum type of all supported options for a PWin card."""

    def __str__(self):
        return str(self.value)

    @classmethod
    def from_string(cls, s: str):
        """
        :param s: String
        :return: SupportedOptions
        """
        for m in cls:
            if m.value.lower() == s.lower():
                return m
        raise ValueError(f"Can't interpret option {s}.")


class AtomicPositionsCard(InputCard):
    """ATOMIC_POSITIONS card"""

    class AtomicPositionsOptions(CardOptions):
        alat = "alat"
        bohr = "bohr"
        angstrom = "angstrom"
        crystal = "crystal"
        crystal_sg = "crystal_sg"

    name = "atomic_positions"
    required = True
    opts = AtomicPositionsOptions
    default_option = opts.alat
    default_deprecated = True

    def __init__(self, option, symbols, positions):
        self.option = option
        self.symbols = symbols
        self.positions = positions

    @classmethod
    def from_string(cls, s: str):
        """Parse a string containing an ATOMIC_POSITIONS card"""
        option, body = cls.split_card_string(s)
        symbols = [line[0] for line in body]
        positions = [line[1:] for line in body]
        return cls(option, symbols, positions)


class CellParametersCard(InputCard):
    """CELL_PARAMETERS card"""

    class CellParametersOptions(CardOptions):
        alat = "alat"
        bohr = "bohr"
        angstrom = "angstrom"

    name = "cell_parameters"
    required = False
    opts = CellParametersOptions
    default_option = opts.alat
    default_deprecated = True

    def __init__(self, option, a1, a2, a3):
        self.option = option
        self.a1, self.a2, self.a3 = a1, a2, a3

    @classmethod
    def from_string(cls, s: str):
        """Parse a string containing an ATOMIC_SPECIES card"""
        option, body = cls.split_card_string(s)
        a1, a2, a3 = body[0], body[1], body[2]
        return cls(option, a1, a2, a3)
