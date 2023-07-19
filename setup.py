# Copyright (c) Omar A. Ashour
# Distributed under the terms of the MIT License.

from setuptools import setup, find_packages

import os

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(SETUP_PTH, "README.md"), "r") as f:
    desc = f.read()


setup(
    name="vesta-espresso",
    packages=find_packages(),
    version="0.0.1",
    install_requires=[
        "f90nml",
    ],
    extras_require={},
    package_data={},
    author="Omar A. Ashour",
    author_email="ashour@berkeley.edu",
    maintainer="Omar A. Ashour",
    url="https://github.com/oashour/vesta-espresso.git",
    license="MIT",
    description="An extremely tiny tool for using VESTA with Quantum Espresso",
    long_description=desc,
    keywords=["quantum-espresso", "vasp", "DFT", "quantum-chemistry", "materials-science"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License" "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    entry_points={
        "console_scripts": [
            "vesta-espresso = vesta_espresso.cli:main",
        ],
    },
)
