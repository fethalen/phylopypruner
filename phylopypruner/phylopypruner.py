#!/usr/bin/env python3

"Module for running PhyloPyPruner by just invoking this script."

import os
import sys

arguments = []
for argument in sys.argv[1:]:
    if os.path.exists(argument):
        argument = os.path.abspath(argument)
    arguments.append(argument)

os.chdir("..")
os.system("{} -m phylopypruner {}".format(sys.executable, " ".join(arguments)))
