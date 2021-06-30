"Build script for setuptools."

from __future__ import absolute_import
from setuptools import setup, find_packages
from distutils.util import convert_path

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

with open(convert_path("phylopypruner/VERSION")) as version_file:
    version_no = version_file.read().strip()

setup(
    name="phylopypruner",
    version=version_no,
    author="Felix Thalen",
    author_email="felix.thalen@uni-goettingen.de",
    license="GPL 3",
    description="tree-based orthology inference",
    # scripts=['scripts/orthofinder2phylopypruner'],
    packages=find_packages(),
    entry_points={"console_scripts": ["phylopypruner = phylopypruner.__main__:entry"]},
    install_requires=[
        "setuptools>=30.3.0",
        "wheel",
        "setuptools_scm"
    ],
    python_requires='>=3.6',
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/fethalen/phylopypruner",
    keywords=["orthology inference", "orthologs", "tree-based",
              "phylogenetics", "phylogenomics", "orthology"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
