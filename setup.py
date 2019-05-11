"Build script for setuptools."

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name="phylopypruner",
    version="0.6.5",
    author="Felix Thalen",
    author_email="fe1430th-s@student.lu.se",
    license="MIT",
    description="tree-based orthology inference",
    scripts=['scripts/phylopypruner'],
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/fethalen/phylopypruner",
    keywords=["orthology inference", "orthologs", "tree-based",
              "phylogenetics", "phylogenomics", "orthology"],
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
