# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
#
# """
# modeling_morphogenesis.py
# A short description of the project.
#
# Handles the primary functions
# """
#
# import sys
# import argparse
#
#
# def warning(*objs):
#     """Writes a message to stderr."""
#     print("WARNING: ", *objs, file=sys.stderr)
#
#
# def canvas(with_attribution=True):
#     """
#     Placeholder function to show example docstring (NumPy format)
#
#     Replace this function and doc string for your own project
#
#     Parameters
#     ----------
#     with_attribution : bool, Optional, default: True
#         Set whether or not to display who the quote is from
#
#     Returns
#     -------
#     quote : str
#         Compiled string including quote and optional attribution
#     """
#
#     quote = "The code is but a canvas to our imagination."
#     if with_attribution:
#         quote += "\n\t- Adapted from Henry David Thoreau"
#     return quote
#
#
# def parse_cmdline(argv):
#     """
#     Returns the parsed argument list and return code.
#     `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
#     """
#     if argv is None:
#         argv = sys.argv[1:]
#
#     # initialize the parser object:
#     parser = argparse.ArgumentParser()
#     # parser.add_argument("-i", "--input_rates", help="The location of the input rates file",
#     #                     default=DEF_IRATE_FILE, type=read_input_rates)
#     parser.add_argument("-n", "--no_attribution", help="Whether to include attribution",
#                         action='store_false')
#     args = None
#     try:
#         args = parser.parse_args(argv)
#     except IOError as e:
#         warning("Problems reading file:", e)
#         parser.print_help()
#         return args, 2
#
#     return args, 0
#
#
# def main(argv=None):
#     args, ret = parse_cmdline(argv)
#     if ret != 0:
#         return ret
#     print(canvas(args.no_attribution))
#     return 0  # success
#
#
# if __name__ == "__main__":
#     status = main()
#     sys.exit(status)

#Ensuring compliance of code with both python2 and python3

from __future__ import division, print_function
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_context('talk')

import pyNetLogo

#Import the sampling and analysis modules for a Sobol variance-based
#sensitivity analysis
from SALib.sample import saltelli
from SALib.analyze import sobol

problem = {
    'num_vars': 6,
    'names': ['random-seed',
              'grass-regrowth-time',
              'sheep-gain-from-food',
              'wolf-gain-from-food',
              'sheep-reproduce',
              'wolf-reproduce'],
    'bounds': [[1, 100000],
               [20., 40.],
               [2., 8.],
               [16., 32.],
               [2., 8.],
               [2., 8.]]
}

n = 1000
param_values = saltelli.sample(problem, n, calc_second_order=True)

param_values.shape

