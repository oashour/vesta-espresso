import argparse
import logging
import os
import shlex
import shutil
import subprocess
import sys
from argparse import RawTextHelpFormatter
from importlib.metadata import metadata
from inspect import cleandoc
from logging import FileHandler, StreamHandler

from vesta_espresso import PWinStructure

# Metadata
metadata = metadata("vesta-espresso")
__author__ = metadata["Author"]
__version__ = metadata["Version"]
__email__ = metadata["Author-email"]
__maintainer__ = metadata["Maintainer"]
__maintainer_email__ = metadata["Maintainer-email"]
__summary__ = metadata["Summary"]
__URL__ = metadata["Home-page"]


def _setup_tmpdir():
    """
    Set up a temporary directory $HOME/.vesta-espresso
    to store temporary files
    """
    tmpdir = os.path.join(os.getenv("HOME"), ".vesta-espresso")
    if not os.path.exists(tmpdir):
        logging.debug(f"Creating temporary directory {tmpdir}")
        os.makedirs(tmpdir)
    scratch = os.path.join(tmpdir, "scratch")
    if not os.path.exists(scratch):
        logging.debug(f"Creating scratch directory {scratch}")
        os.makedirs(scratch)
    return tmpdir


def _get_parser():
    """Sets up the command line parser"""
    parser = argparse.ArgumentParser(
        description=f"""{__summary__}""",
        epilog=cleandoc(
            f"""
    Author: {__author__} ({__email__}), 
    Version: {__version__}, 
    Maintainer: {__maintainer__} ({__email__}), 
    URL: {__URL__}"""
        ),
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument("filenames", nargs="*", help="Files to open in VESTA")

    parser.add_argument(
        "--version",
        action="store_true",
        help="Return the version of vesta-espresso",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output",
    )

    return parser


def _get_vesta():
    """Finds the VESTA executable"""
    if (vesta := shutil.which("VESTA")) is not None:
        logging.debug(f"VESTA found at {vesta}")
    else:
        # TODO: make it work for other OSes?
        vesta = "/Applications/VESTA/VESTA.app/Contents/MacOS/VESTA"
        logging.debug(f"VESTA not found in path. Using {vesta}")
    return vesta


def _setup_logger(verbose):
    """Sets up file and stream logger"""
    logger = logging.getLogger()
    formatter = logging.Formatter("%(asctime)s |  %(levelname)s: %(message)s")
    logger.setLevel(logging.DEBUG)

    stream_handler = StreamHandler()
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)

    logFilePath = ".vesta-espresso.log"
    file_handler = FileHandler(filename=logFilePath, mode="w")
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


def create_poscars(filenames, tmpdir):
    """
    Creates temporary POSCARs from .in or .pwi files
    Doesn't touch non-espresso files
    """
    temp_poscars = []
    for i, f in enumerate(filenames):
        if os.path.splitext(f)[1] in [".in", ".pwi"]:
            poscar = os.path.join(tmpdir, "scratch", f".{os.path.basename(f)}.vasp")
            temp_poscars.append(poscar)
            logging.debug("Creating %s from %s", poscar, f)
            pwin = PWinStructure.from_file(f)
            with open(poscar, "w") as f:
                f.write(str(pwin))
            filenames[i] = poscar

    # Quote them to handle spaces in filenames
    filenames = [f'"{f}"' for f in filenames]
    return filenames, temp_poscars


def main():
    args = _get_parser().parse_args()

    _setup_logger(args.verbose)

    if args.version:
        sys.exit(f"vesta-espresso version {__version__}")

    vesta = _get_vesta()
    tmpdir = _setup_tmpdir()

    if filenames := args.filenames:
        # Get extensions of all filenames
        extensions = [os.path.splitext(filename)[1] for filename in filenames]
        # Check if .in or .pwi files are present
        if ".in" in extensions or ".pwi" in extensions:
            logging.debug("Detected some QE files. Creating temporary POSCARs.")
            filenames, temp_poscars = create_poscars(filenames, tmpdir)
        else:
            logging.debug('Running "vesta" with provided files as is.')
        command = f"{vesta} {' '.join(filenames)}"
    else:
        command = vesta

    logging.debug(f"Command: {command}")
    command = shlex.split(command)

    subprocess.Popen(command)


if __name__ == "__main__":
    main()
