Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/2_GW_Neff/2_GW_Neff_'
Root_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'

ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/PBH_posteriors_Full.pdf'
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
    [samples_3, samples_3, samples_3],
    ['LgMc', 'LgFbh', 'Sbh'],
    width_inch=12,
    contour_colors=['g', 'r', 'b'],
    legend_labels=['GW', 'GW + $\Delta N_{\mathrm{eff}}$', 'GW + $\Delta N_{\mathrm{eff}}$ + PBH'],
    filled = True,
    line_args=[
        {'lw':1.5,'ls':'-', 'color':'g'},
        {'lw':1.5,'ls':'-', 'color':'r'},
        {'lw':1.5,'ls':'-', 'color':'b'},

        ],
    title_limit=2,
    param_limits = {
        'LgFbh': [-20, 0]
        })
g.export(ResultFile)
print(ResultFile)

from PyLab import *
r = Getdist_Marg_Stat(Root_2, 'LgFbh')
print(r)
