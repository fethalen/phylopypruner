#!/usr/bin/env python

"Module for running PhyloPyPruner by just invoking this script."

import os
import sys

ARGUMENTS = []
for argument in sys.argv[1:]:
    if os.path.exists(argument):
        argument = os.path.abspath(argument)
    ARGUMENTS.append(argument)

os.chdir("..")
os.system("{} -m phylopypruner {}".format(sys.executable, " ".join(ARGUMENTS)))
