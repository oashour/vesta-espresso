"""
Classes for reading/manipulating PWscf xml files.
"""


from copy import copy
import itertools
import pathlib

import f90nml

from vesta_espresso.utils import ibrav_to_lattice
from vesta_espresso.cards import (
    AtomicPositionsCard,
    CellParametersCard,
)

bohr_to_ang = 0.529177210903


class PWinStructure:
    """
    Class for PWscf input files
    """

    # Names of all possible cards
    card_names = [
        "atomic_species",
        "atomic_positions",
        "k_points",
        "additional_k_points",
        "cell_parameters",
        "constraints",
        "occupations",
        "atomic_velocities",
        "atomic_forces",
        "solvents",
        "hubbard",
    ]

    def __init__(self, system, atomic_positions, cell_parameters):
        # self.system = system
        # self.cell_parameters = cell_parameters
        # self.atomic_positions = atomic_positions

        ibrav, alat, celldm = self.get_lattice_info(system, atomic_positions, cell_parameters)
        self.lattice_matrix = self.get_lattice_matrix(cell_parameters, ibrav, alat, celldm)
        self.species, self.coords, self.coords_are_cartesian = self.get_species_coords(
            atomic_positions, alat
        )

    @classmethod
    def from_file(cls, filename):
        """
        Reads an inputfile from file

        Args:
            filename: path to file

        Returns:
            PWin object
        """
        parser = f90nml.Parser()
        parser.comment_tokens += "#"

        pwi_str = pathlib.Path(filename).read_text()
        system = parser.reads(pwi_str).get("system", None)
        if system is None:
            raise PWinParserError("No &SYSTEM namelist found.")

        atomic_positions, cell_parameters = cls._parse_cards(pwi_str)
        if atomic_positions is None:
            raise PWinParserError("No ATOMIC_POSITIONS card found.")

        return cls(system, atomic_positions, cell_parameters)

    @classmethod
    def _parse_cards(cls, pwi_str):
        card_strings = pwi_str.rsplit("/", 1)[1].split("\n")
        card_strings = [c for c in card_strings if c]
        card_idx = [
            i
            for i, string in enumerate(card_strings)
            if string.split()[0].lower() in cls.card_names
        ]

        atomic_positions, cell_parameters = None, None
        for i, j in zip(card_idx, card_idx[1:] + [None]):  # type: ignore
            card_name = card_strings[i].split()[0].lower()
            if card_name == "atomic_positions":
                card_string = "\n".join(card_strings[i:j])
                atomic_positions = AtomicPositionsCard.from_string(card_string)
            elif card_name == "cell_parameters":
                card_string = "\n".join(card_strings[i:j])
                cell_parameters = CellParametersCard.from_string(card_string)

        return atomic_positions, cell_parameters

    @staticmethod
    def get_lattice_info(system, atomic_positions, cell_parameters):
        """
        Returns:
            ibrav (int): Bravais lattice type
            alat (float): lattice constant
            celldm (list): list of lattice parameters
        """

        def get_celldm_from_ABC():
            """Converts ABC to celldm"""
            # A is already in angstrom
            B = system.get("B", 0)
            C = system.get("C", 0)
            cosAB = system.get("cosAB", 0)
            cosAC = system.get("cosAC", 0)
            cosBC = system.get("cosBC", 0)
            return [A, B / A, C / A, cosBC, cosAC, cosAB]

        try:
            ibrav = system["ibrav"]
        except KeyError as e:
            raise PWinParserError("ibrav must be set in system namelist") from e
        if ibrav == 0 and cell_parameters is None:
            raise PWinParserError("cell_parameters must be set if ibrav=0")

        alat_required = (
            cell_parameters is not None and cell_parameters.option == CellParametersCard.opts.alat
        ) or (atomic_positions.option == AtomicPositionsCard.opts.alat)

        celldm = copy(system.get("celldm", None))
        A = system.get("A", None)
        if celldm is None and A is None and alat_required:
            raise PWinParserError(
                "either celldm(1) or A must be set if any cards options are alat."
            )
        if celldm is not None and A is not None:
            raise PWinParserError("celldm(1) and A cannot both be set.")
        alat = celldm[0] * bohr_to_ang if celldm is not None else A

        if celldm is not None:
            celldm[0] *= bohr_to_ang  # celldm(1) is originally in bohr
            # Get it to the right length since not all 6 are required in input
            celldm = celldm + [0] * (6 - len(celldm))
        elif A is not None:
            celldm = get_celldm_from_ABC()

        return ibrav, alat, celldm

    @staticmethod
    def get_species_coords(atomic_positions, alat):
        """
        Returns:
            species (list): list of species symbols
            coords (list): list of atomic coordinates
            coords_are_cartesian (bool): True if coordinates are cartesian (in Angstrom)
        """
        option = atomic_positions.option
        species = atomic_positions.symbols
        coords = atomic_positions.positions
        if option == AtomicPositionsCard.opts.alat:
            coords = [[xi * alat for xi in c] for c in coords]
            coords_are_cartesian = True
        elif option == AtomicPositionsCard.opts.bohr:
            coords = [[xi * bohr_to_ang for xi in c] for c in coords]
            coords_are_cartesian = True
        elif option == AtomicPositionsCard.opts.angstrom:
            coords_are_cartesian = True
        elif option == AtomicPositionsCard.opts.crystal:
            coords_are_cartesian = False
        elif option == AtomicPositionsCard.opts.crystal_sg:
            raise PWinParserError("Atomic positions with crystal_sg option are not supported.")

        return species, coords, coords_are_cartesian

    @staticmethod
    def get_lattice_matrix(cell_parameters, ibrav, alat, celldm):
        """
        Returns:
            Lattice object (in ANGSTROM no matter what's in the input file)
        """
        if ibrav != 0:
            return ibrav_to_lattice(ibrav, celldm)

        lattice_matrix = [
            cell_parameters.a1,
            cell_parameters.a2,
            cell_parameters.a3,
        ]

        option = cell_parameters.option
        if option == CellParametersCard.opts.alat:
            lattice_matrix = [[x * alat for x in ai] for ai in lattice_matrix]
        elif option == CellParametersCard.opts.bohr:
            lattice_matrix = [[x * bohr_to_ang for x in ai] for ai in lattice_matrix]
        elif option != CellParametersCard.opts.angstrom:
            raise PWinParserError(
                "cell_parameters option must be one of 'alat', 'bohr', or 'angstrom'. "
                f"{option} is not supported."
            )
        return lattice_matrix

    def __str__(self):
        """
        Returns a POSCAR string from the QE input file.

        Returns:
            String representation of POSCAR.
        """
        site_symbols = list(set(self.species))
        natoms = [self.species.count(i) for i in site_symbols]
        species, coords = zip(
            *sorted(zip(self.species, self.coords), key=lambda x: site_symbols.index(x[0]))
        )
        # natoms = [len(tuple(s[1])) for s in itertools.groupby(self.species)]
        # site_symbols = [s[0] for s in itertools.groupby(self.species)]

        format_str = "{:21.16f}"
        lines = ["Temporary POSCAR generated by vesta-espresso", "1.0"]
        lines.extend(" ".join(format_str.format(c) for c in v) for v in self.lattice_matrix)
        lines.extend(
            (
                " ".join(site_symbols),
                " ".join(map(str, natoms)),
                "cartesian" if self.coords_are_cartesian else "direct",
            )
        )
        for s, c in zip(species, coords):
            line = " ".join(format_str.format(xi) for xi in c)
            line += f" {s}"
            lines.append(line)

        return "\n".join(lines) + "\n"


class PWinParserError(Exception):
    """
    Exception class for PWin parsing.
    """
