#!/usr/bin/env bash
#
# Script for running and uploading new versions of PhyloPyPruner.

rm -rf build dist phylopypruner.*
python3 setup.py sdist bdist_wheel
# clear cache
curl -X PURGE https://pypi.python.org/pypi/canonicalwebteam-yaml-redirects/
python3 -m twine upload dist/*
