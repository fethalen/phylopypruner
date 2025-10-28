"Build script for setuptools."

from __future__ import absolute_import

import os
from distutils.util import convert_path

from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

here = os.path.abspath(os.path.dirname(__file__))
version_path = os.path.join(here, "phylopypruner", "VERSION")
with open(version_path) as version_file:
    version_no = version_file.read().strip()

setup(
    name="phylopypruner",
    version=version_no,
    author="Felicia Sandberg",
    author_email="felicia.sandberg@hest.ethz.ch",
    license="GPL 3",
    description="Tree-based orthology inference",
    # scripts=['scripts/orthofinder2phylopypruner'],
    packages=find_packages(),
    entry_points={"console_scripts": ["phylopypruner = phylopypruner.__main__:entry"]},
    install_requires=["setuptools>=30.3.0", "wheel", "setuptools_scm"],
    python_requires=">=3.6",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/fethalen/phylopypruner",
    keywords=[
        "orthology inference",
        "orthologs",
        "tree-based",
        "phylogenetics",
        "phylogenomics",
        "orthology",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
