#!/usr/bin/env python

"Module for running PhyloPyPruner by just invoking this script."

import os
import sys

os.chdir("..")
os.system("python -m phylopypruner {}".format(" ".join(sys.argv[1:])))
