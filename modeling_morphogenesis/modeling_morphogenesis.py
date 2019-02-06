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
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_context('talk')
import jpype

import pyNetLogo
import os
# os.chdir('/Users/agnesresto/Documents/NetLogo 6.0.4')

# print(os.environ)
# os.environ["JAVA\_HOME"] = "//Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java python setup.py install"

# print(os.environ["JAVA\_HOME"]) #makes sure that java home system variable was set correctly

#Import the sampling and analysis modules for a Sobol variance-based
#sensitivity analysis
from SALib.sample import saltelli
from SALib.analyze import sobol

import pyNetLogo
netlogo = pyNetLogo.NetLogoLink(gui=False,netlogo_home = '/Users/agnesresto/Documents/NetLogo 6.0.4')
netlogo.load_model('./models/Wolf Sheep Predation_v6.nlogo')

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

results = pd.DataFrame(columns=['Avg. sheep', 'Avg. wolves'])

# import ipyparallel
#
# client = ipyparallel.Client()
# client.ids
#
# direct_view = client[:]

# import os
#
# #Push the current working directory of the notebook to a "cwd" variable on the engines that can be accessed later
# direct_view.push(dict(cwd=os.getcwd()))
#
# #Push the "problem" variable from the notebook to a corresponding variable on the engines
# direct_view.push(dict(problem=problem))
#
# %%px
#
# import pyNetLogo
# netlogo = pyNetLogo.NetLogoLink(gui=False,netlogo_home = '/Users/agnesresto/Documents/NetLogo 6.0.4') #jvm_home = '/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home/
#netlogo = pyNetLogo.NetLogoLink(gui=False, netlogo_home = '/Users/agnesresto/Documents/NetLogo 6.0.4', netlogo_version = '6') #, jvm_home = '/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home'#


# def simulation(experiment):
#
#     #Set the input parameters
#     for i, name in enumerate(problem['names']):
#         if name == 'random-seed':
#             #The NetLogo random seed requires a different syntax
#             netlogo.command('random-seed {}'.format(experiment[i]))
#         else:
#             #Otherwise, assume the input parameters are global variables
#             netlogo.command('set {0} {1}'.format(name, experiment[i]))
#
#     netlogo.command('setup')
#     #Run for 100 ticks and return the number of sheep and wolf agents at each time step
#     counts = netlogo.repeat_report(['count sheep','count wolves'], 100)
#
#     results = pd.Series([counts['count sheep'].values.mean(),
#                          counts['count wolves'].values.mean()],
#                         index=['Avg. sheep', 'Avg. wolves'])
#
#     return results

# lview = client.load_balanced_view()

# results = pd.DataFrame(lview.map_sync(simulation, param_values))

import time

t0=time.time()

for run in range(param_values.shape[0]):

    #Set the input parameters
    for i, name in enumerate(problem['names']):
        if name == 'random-seed':
            #The NetLogo random seed requires a different syntax
            netlogo.command('random-seed {}'.format(param_values[run,i]))
        else:
            #Otherwise, assume the input parameters are global variables
            netlogo.command('set {0} {1}'.format(name, param_values[run,i]))

    netlogo.command('setup')
    #Run for 100 ticks and return the number of sheep and wolf agents at each time step
    counts = netlogo.repeat_report(['count sheep','count wolves'], 100)

    #For each run, save the mean value of the agent counts over time
    results.loc[run, 'Avg. sheep'] = counts['count sheep'].values.mean()
    results.loc[run, 'Avg. wolves'] = counts['count wolves'].values.mean()

elapsed=time.time()-t0 #Elapsed runtime in seconds

elapsed


# results = pd.DataFrame(simulation, param_values)

# results.to_csv('./data/Sobol_parallel.csv')

results.head(5)