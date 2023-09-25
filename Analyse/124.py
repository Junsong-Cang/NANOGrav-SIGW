Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/2_GW_Neff/2_GW_Neff_'
Root_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'

ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/param_posteriors.pdf'
ResultFile='/Users/cangtao/Desktop/param_posteriors_tmp.pdf'

# ---- Initialise ----
import getdist, os
from getdist import plots
import matplotlib.pyplot as plt

# ---- Getdist Plot ----
plt.rcParams.update({'font.family':'Times'})

samples_1 = getdist.mcsamples.loadMCSamples(Root_1)
samples_2 = getdist.mcsamples.loadMCSamples(Root_2)
samples_3 = getdist.mcsamples.loadMCSamples(Root_3)

p = samples_1.getParams()
g = plots.getSubplotPlotter(subplot_size = 3)
g.settings.axes_fontsize=15
g.settings.title_limit_fontsize = 15
g.settings.lab_fontsize =15

g.triangle_plot(
    [samples_3, samples_2, samples_1],
    # ['LgK', 'LgA', 'S', 'LgMc', 'LgFbh', 'Sbh', 'dNeff'],
    ['LgK', 'LgA', 'S', 'dNeff'],
    width_inch=12,
    contour_colors=['g', 'r', 'b'],
    # legend_labels=['GW', 'GW + $\Delta N_{\mathrm{eff}}$', 'GW + $\Delta N_{\mathrm{eff}}$ + PBH'],
    filled = True,
    line_args=[
        {'lw':1.5,'ls':'-', 'color':'g'},
        {'lw':1.5,'ls':'-', 'color':'r'},
        {'lw':1.5,'ls':'-', 'color':'b'},
        ],
    title_limit = 2,
    markers={'dNeff': 0.175},
    param_limits = {'dNeff': [0, 0.5]}
    )

g.export(ResultFile)
print(ResultFile)

from PyLab import *

Root = Root_3

LgA = Getdist_Marg_Stat(Root, 'LgA')
LgK = Getdist_Marg_Stat(Root, 'LgK')
S = Getdist_Marg_Stat(Root, 'S')
Fbh = Getdist_Marg_Stat(Root, 'LgFbh')
Mc = Getdist_Marg_Stat(Root, 'LgMc')
Sbh = Getdist_Marg_Stat(Root, 'Sbh')
dNeff = Getdist_Marg_Stat(Root, 'dNeff')

print('LgA -- ', LgA['low_95'], LgA['upper_95'])
print('LgK -- ', LgK['low_95'], LgK['upper_95'])
print('S -- ', S['low_95'], S['upper_95'])
print('LgFbh -- ', Fbh['low_95'], Fbh['upper_95'])
print('LgMc -- ', Mc['low_95'], Mc['upper_95'])
print('Sbh -- ', Sbh['low_95'], Sbh['upper_95'])
print('dNeff -- ', dNeff['low_95'], dNeff['upper_95'])

samples = samples_3
best_fit = samples.getLikeStats()
print(best_fit)
