#!/usr/bin/env bash
#
# Script for running and uploading new versions of PhyloPyPruner.

rm -rf build dist venv ../phylopypruner.*
python3 setup.py sdist bdist_wheel
curl -X PURGE https://pypi.python.org/pypi/canonicalwebteam-yaml-redirects/
curl -X PURGE https://pypi.python.org/simple/phylopypruner
python3 -m twine upload dist/*
