# vesta-espresso

[VESTA](https://jp-minerals.org/vesta/en/) is an extremely powerful and useful tool heavily used in the condensed matter physics, quantum chemistry and materials science communities. However, it does not natively support viewing Quantum ESPRESSO (QE) input files (which contain the structural information). Unfortunately, since VESTA is closed-source, third-party developers cannot add this functionality to VESTA.

This repository contains a simple and ultra-lightweight Python package that allows you to open Quantum ESPRESSO input files in VESTA. (We consider files with `.pwi` or `.in` extensions to be QE files). If you're on a Mac, you can simply download the app from the [releases](https://github.com/oashour/vesta-espresso/releases) page and use VESTA to open QE files from Finder like you would with a CIF or VASP POSCAR. If not, you can install the command line tool with `pip` and open the files in VESTA by simply typing `vesta-espresso <input file>` in your terminal. See below for installation instructions.

The script and app can accept both Quantum ESPRESSO input files and other VESTA-supported formats (which will be passed to VESTA as is without any processing), or any combination thereof. QE files are converted to a hidden temporary VASP POSCAR file, which itself is opened in VESTA. We do not perform any symmetry analysis or reduction, simply format and unit conversion.

Currently, `vesta-espresso` supports all types and units of structure specification (`alat`, `crystal`, `angstrom`, etc.), including all `ibrav != 0`, except for `crystal_sg` (i.e., specifying the space group and inequivalent atoms).

Our QE input to POSCAR conversion is based on a stripped-down and extremely minimalist version of [pymatgen-io-espresso](https://github.com/oashour/pymatgen-io-espresso), with [f90nml](https://github.com/marshallward/f90nml), an ultralight Fortran namelist parser, as the sole dependency. Using `pymatgen-io-espresso` with this project (as in v0.0.1) lead to a huge package/app due to the Pymatgen dependency (700+ MB). On the other hand, installing `vesta-espresso` with `pip` takes up about 300 kB, and the entire Mac app is 2.5 MB and has no external dependencies, just drag and drop.

## Requirements

1. [VESTA](https://jp-minerals.org/vesta/en/). I tested with version 3.5.6, but any decently modern version should work.
2. If using the command-line package only, you need Python 3.8 or later. If using the Mac app, you don't need python.

## Installation Instructions

# Mac

These instructions assume you have already downloaded [VESTA](https://jp-minerals.org/vesta/en/) and installed it in the default location (`/Applications/VESTA.app`).

Download the `.dmg` file from the [releases](https://github.com/oashour/vesta-espresso/releases) page and install the app as usual. If you add this to your Dock, you'll be able to drag Quantum ESPRESSO input files with the extensions `.in` or `.pwi` onto the icon to open them in VESTA, or right-click them in Finder and select 'Open With -> Vesta (ESPRESSO).app'. You can also change the default app to VESTA (ESPRESSO).app as well.

Note that you can't use the `open` dialogue within VESTA to open Quantum ESPRESSO files. The app only acts as a wrapper around VESTA and doesn't add any functionality to the program itself.

If you'd like to use `vesta-espresso` from the command line, simply add the following to your `~/.bashrc` or equivalent:

```bash
alias vesta="/Applications/VESTA/VESTA\ \(ESPRESSO\).app/Contents/Resources/vesta-espresso.dist/vesta-espresso"
```

Now you can call VESTA from the command line with both Quantum ESPRESSO files and any other natively supported files:

```bash
vesta file1.in file2.pwi file3.cif file4.vasp
```

## Command Line Installation

If you're not using Linux or you just want to use the script without the Mac App (i.e., no Finder integration), you can install it with `pip`. I haven't tested this on Windows.

```bash
# Activate venv, or user pip install --user, or whatever you prefer
pip install git+https://github.com/oashour/vesta-espresso.git
```

For Linux, you need to add the VESTA binary to your path (it should have the default name `VESTA`)

```bash
export PATH=/path/to/vesta/bin:$PATH
```

On Mac, you don't need to do anything if VESTA is installed in the default location (`/Applications/VESTA/VESTA.app`). If it's not, adding it to the path won't work (see this [bug](https://groups.google.com/g/vesta-discuss/c/Cq1_QJwrvhU/m/bU_GYBemBgAJ)), so you'll need to create an Alias

```bash
alias VESTA="/path/to/VESTA.app/Contents/MacOS/VESTA"
```
