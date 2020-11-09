"Build script for setuptools."

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
    # scripts=['scripts/phylopypruner'],
    entry_points={
        "console_scripts": [
            "phylopypruner=phylopypruner.__main__:phylopypruner",
        ],
    },
    install_requires=[
        "setuptools>=30.3.0",
        "wheel",
        "setuptools_scm"
    ],
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/fethalen/phylopypruner",
    keywords=["orthology inference", "orthologs", "tree-based",
              "phylogenetics", "phylogenomics", "orthology"],
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
