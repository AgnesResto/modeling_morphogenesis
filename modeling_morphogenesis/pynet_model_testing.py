from __future__ import division, print_function
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import time
from shapely import geometry
from shapely.geometry import Polygon
import seaborn as sns
sns.set_style('white')
sns.set_context('talk')


# Import the sampling and analysis modules for a Sobol variance-based
# sensitivity analysis
from SALib.sample import saltelli
from SALib.analyze import sobol

import pyNetLogo
netlogo = pyNetLogo.NetLogoLink(gui=False,netlogo_home = '/Users/agnesresto/Documents/NetLogo 6.0.4')
netlogo.load_model('/Users/agnesresto/modeling_morphogenesis/modeling_morphogenesis/modeling_morphogenesis/Morphogenesis_3Denv.nlogo')

problem = {
    'num_vars': 2,
    'names': ['random-seed',
              'num-cells'],
    'bounds': [[1, 100000],
               [5, 15]]
}

n = 1

param_values_noround = saltelli.sample(problem, n, calc_second_order=True)
param_values = np.around(param_values_noround, decimals=0)
param_values.shape

print(param_values)

results = pd.DataFrame(columns=['Num Cysts', 'Num Differentiated', 'Num Pluripotent'])


t0 = time.time()

for run in range(param_values.shape[0]):

    # Set the input parameters
    for i, name in enumerate(problem['names']):
        if name == 'random-seed':
            # The NetLogo random seed requires a different syntax
            netlogo.command('random-seed {}'.format(param_values[run, i]))
        else:
            # Otherwise, assume the input parameters are global variables
            netlogo.command('set {0} {1}'.format(name, param_values[run, i]))

    netlogo.command('setup')
    netlogo.repeat_command('go', 300)
    # Run for 300 ticks to ensure model has finished running.

    cyst_num = netlogo.report('num-cysts')
    cyst_num_int = np.array(cyst_num, dtype=int, ndmin=2)

    differentiated = netlogo.report('count cells with [color = green]')
    pluripotent = netlogo.report('count cells with [color = red]')
    total_cells = netlogo.report('count cells')
    total_cells_int = np.array(total_cells, dtype=int, ndmin=1)

    print(total_cells_int)
    xcor = netlogo.report('map [s -> [xcor] of s] sort cells')
    ycor = netlogo.report('map [s -> [ycor] of s] sort cells')
    group_id = netlogo.report('map [s -> [group-id] of s] sort cells')

    cir_x = np.zeros((len(xcor), int(cyst_num_int[0])-1))
    cir_x.fill(np.nan)
    # -1 because the last value of cyst_num in netlogo does not count
    cir_y = np.zeros((len(xcor), int(cyst_num_int[0])-1))
    cir_y.fill(np.nan)

    for j in range(int(cyst_num_int[0])):
        n = 0
        pointList = []
        for k, num in enumerate(group_id):
            if j == num:
                cir_x[n, j-1] = xcor[k]
                cir_y[n, j-1] = ycor[k]
                pnt = geometry.Point(xcor[k], ycor[k])
                pointList.append(pnt)
                n = n + 1
        # cyst_x = cir_x
        # cyst_y
        circle_x = cir_x[:, j-1]
        circle_y = cir_y[:, j-1]
        tup_corx = circle_x[~np.isnan(circle_x)]
        tup_cory = circle_y[~np.isnan(circle_y)]
        tup_cor = tuple(map(tuple, (tup_corx, tup_cory)))
            #cyst_area = Polygon([(tup_corx[:], tup_cory[:])]).area


    poly = geometry.Polygon([[p.x, p.y] for p in pointList])
    print(poly.wkt)
    print(poly.area)

        #tup_cor = tuple(map(tuple, (cir_x[:,j-1], cir_y[:,j-1])))
        #poly = geometry.Polygon([[p.x, p.y] for p in pointList])
        #print(poly.wkt)
        #cyst_area = Polygon([(cir_x[:,:], cir_y[:,:])]).area
    #(Polygon([(0,0), (4,0), (2,4)]).area)
    #print(cir_x)
    #tuple(map(tuple, cir_x))

    #area=feature.geometry().area()

    # For each model run, save the number of cysts and number of
    # differentiated and undifferentiated
    results.loc[run, 'Num Cysts'] = netlogo.report('num-cysts')
    results.loc[run, 'Num Differentiated'] = netlogo.report('count cells with [color = green]')
    results.loc[run, 'Num Pluripotent'] = netlogo.report('count cells with [color = red]')

elapsed=time.time()-t0 # Elapsed runtime in seconds

elapsed

# results = pd.DataFrame(simulation, param_values)

# results.to_csv('./data/Sobol_parallel.csv')

print(results.head(5))

# sns.set_style('white')
# sns.set_context('talk')
# fig, ax = plt.subplots(1,len(results.columns), sharey=True)
#
# for i, n in enumerate(results.columns):
#     ax[i].hist(results[n], 20)
#     ax[i].set_xlabel(n)
#
# ax[0].set_ylabel('Counts')
#
# fig.set_size_inches(10,4)
# fig.subplots_adjust(wspace=0.1)
# #plt.savefig('JASSS figures/SA - Output distribution.pdf', bbox_inches='tight')
# #plt.savefig('JASSS figures/SA - Output distribution.png', dpi=300, bbox_inches='tight')
# plt.show()
#
# # %matplotlib
# import scipy
#
# nrow=2
# ncol=3
# fig, ax = plt.subplots(nrow, ncol, sharey=True)
# sns.set_context('talk')
# y = results['Num Cysts']
#
# for i, a in enumerate(ax.flatten()):
#     x = param_values[:, i]
#     sns.regplot(x, y, ax=a, ci=None, color='k',scatter_kws={'alpha': 0.2, 's': 4, 'color': 'gray'})
#     pearson = scipy.stats.pearsonr(x, y)
#     a.annotate("r: {:6.3f}".format(pearson[0]), xy=(0.15, 0.85), xycoords='axes fraction',fontsize=13)
#     if divmod(i,ncol)[1]>0:
#         a.get_yaxis().set_visible(False)
#     a.set_xlabel(problem['names'][i])
#     a.set_ylim([0,1.1*np.max(y)])
#
# fig.set_size_inches(9,9,forward=True)
# fig.subplots_adjust(wspace=0.2, hspace=0.3)
# #plt.savefig('JASSS figures/SA - Scatter.pdf', bbox_inches='tight')
# #plt.savefig('JASSS figures/SA - Scatter.png', dpi=300, bbox_inches='tight')
# plt.show()
#
# Si = sobol.analyze(problem, results['Avg. sheep'].values, calc_second_order=True, print_to_console=False)
#
# Si_filter = {k:Si[k] for k in ['ST','ST_conf','S1','S1_conf']}
# Si_df = pd.DataFrame(Si_filter, index=problem['names'])
#
# Si_df
#
# sns.set_style('white')
# fig, ax = plt.subplots(1)
#
# indices = Si_df[['S1','ST']]
# err = Si_df[['S1_conf','ST_conf']]
#
# indices.plot.bar(yerr=err.values.T,ax=ax)
# fig.set_size_inches(8,4)
#
# #plt.savefig('JASSS figures/SA - Indices.pdf', bbox_inches='tight')
# #plt.savefig('JASSS figures/SA - Indices.png', dpi=300, bbox_inches='tight')
#
# plt.show()
#
# import itertools
# from math import pi
#
#
# def normalize(x, xmin, xmax):
#     return (x-xmin)/(xmax-xmin)
#
#
# def plot_circles(ax, locs, names, max_s, stats, smax, smin, fc, ec, lw,
#                  zorder):
#     s = np.asarray([stats[name] for name in names])
#     s = 0.01 + max_s * np.sqrt(normalize(s, smin, smax))
#
#     fill = True
#     for loc, name, si in zip(locs, names, s):
#         if fc=='w':
#             fill=False
#         else:
#             ec='none'
#
#         x = np.cos(loc)
#         y = np.sin(loc)
#
#         circle = plt.Circle((x,y), radius=si, ec=ec, fc=fc, transform=ax.transData._b,
#                             zorder=zorder, lw=lw, fill=True)
#         ax.add_artist(circle)
#
#
# def filter(sobol_indices, names, locs, criterion, threshold):
#     if criterion in ['ST', 'S1', 'S2']:
#         data = sobol_indices[criterion]
#         data = np.abs(data)
#         data = data.flatten() # flatten in case of S2
#         # TODO:: remove nans
#
#         filtered = ([(name, locs[i]) for i, name in enumerate(names) if
#                      data[i]>threshold])
#         filtered_names, filtered_locs = zip(*filtered)
#     elif criterion in ['ST_conf', 'S1_conf', 'S2_conf']:
#         raise NotImplementedError
#     else:
#         raise ValueError('unknown value for criterion')
#
#     return filtered_names, filtered_locs
#
#
# def plot_sobol_indices(sobol_indices, criterion='ST', threshold=0.01):
#     '''plot sobol indices on a radial plot
#
#     Parameters
#     ----------
#     sobol_indices : dict
#                     the return from SAlib
#     criterion : {'ST', 'S1', 'S2', 'ST_conf', 'S1_conf', 'S2_conf'}, optional
#     threshold : float
#                 only visualize variables with criterion larger than cutoff
#
#     '''
#     max_linewidth_s2 = 15#25*1.8
#     max_s_radius = 0.3
#
#     # prepare data
#     # use the absolute values of all the indices
#     #sobol_indices = {key:np.abs(stats) for key, stats in sobol_indices.items()}
#
#     # dataframe with ST and S1
#     sobol_stats = {key:sobol_indices[key] for key in ['ST', 'S1']}
#     sobol_stats = pd.DataFrame(sobol_stats, index=problem['names'])
#
#     smax = sobol_stats.max().max()
#     smin = sobol_stats.min().min()
#
#     # dataframe with s2
#     s2 = pd.DataFrame(sobol_indices['S2'], index=problem['names'],
#                       columns=problem['names'])
#     s2[s2<0.0]=0. #Set negative values to 0 (artifact from small sample sizes)
#     s2max = s2.max().max()
#     s2min = s2.min().min()
#
#     names = problem['names']
#     n = len(names)
#     ticklocs = np.linspace(0, 2*pi, n+1)
#     locs = ticklocs[0:-1]
#
#     filtered_names, filtered_locs = filter(sobol_indices, names, locs,
#                                            criterion, threshold)
#
#     # setup figure
#     fig = plt.figure()
#     ax = fig.add_subplot(111, polar=True)
#     ax.grid(False)
#     ax.spines['polar'].set_visible(False)
#     ax.set_xticks(ticklocs)
#
#     ax.set_xticklabels(names)
#     ax.set_yticklabels([])
#     ax.set_ylim(ymax=1.4)
#     legend(ax)
#
#     # plot ST
#     plot_circles(ax, filtered_locs, filtered_names, max_s_radius,
#                  sobol_stats['ST'], smax, smin, 'w', 'k', 1, 9)
#
#     # plot S1
#     plot_circles(ax, filtered_locs, filtered_names, max_s_radius,
#                  sobol_stats['S1'], smax, smin, 'k', 'k', 1, 10)
#
#     # plot S2
#     for name1, name2 in itertools.combinations(zip(filtered_names, filtered_locs), 2):
#         name1, loc1 = name1
#         name2, loc2 = name2
#
#         weight = s2.ix[name1, name2]
#         lw = 0.5+max_linewidth_s2*normalize(weight, s2min, s2max)
#         ax.plot([loc1, loc2], [1,1], c='darkgray', lw=lw, zorder=1)
#
#     return fig
#
#
# from matplotlib.legend_handler import HandlerPatch
# class HandlerCircle(HandlerPatch):
#     def create_artists(self, legend, orig_handle,
#                        xdescent, ydescent, width, height, fontsize, trans):
#         center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
#         p = plt.Circle(xy=center, radius=orig_handle.radius)
#         self.update_prop(p, orig_handle, legend)
#         p.set_transform(trans)
#         return [p]
#
# def legend(ax):
#     some_identifiers = [plt.Circle((0,0), radius=5, color='k', fill=False, lw=1),
#                         plt.Circle((0,0), radius=5, color='k', fill=True),
#                         plt.Line2D([0,0.5], [0,0.5], lw=8, color='darkgray')]
#     ax.legend(some_identifiers, ['ST', 'S1', 'S2'],
#               loc=(1,0.75), borderaxespad=0.1, mode='expand',
#               handler_map={plt.Circle: HandlerCircle()})
#
#
# sns.set_style('whitegrid')
# fig = plot_sobol_indices(Si, criterion='ST', threshold=0.005)
# fig.set_size_inches(7,7)
# #plt.savefig('JASSS figures/Figure 8 - Interactions.pdf', bbox_inches='tight')
# #plt.savefig('JASSS figures/Figure 8 - Interactions.png', dpi=300, bbox_inches='tight')
# plt.show()
