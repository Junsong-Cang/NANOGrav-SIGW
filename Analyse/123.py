Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/2_GW_Neff/2_GW_Neff_'
Root_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'

ResultFile='/Users/cangtao/Desktop/param_posteriors.pdf'

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
    [samples_1, samples_2, samples_3],
    # ['LgK', 'LgA', 'S', 'LgMc', 'LgFbh', 'Sbh', 'dNeff'],
    ['LgK', 'LgA', 'S', 'dNeff'],
    width_inch=12,
    contour_colors=['g', 'r', 'b'],
    legend_labels=['GW', 'GW + $\Delta N_{\mathrm{eff}}$', 'GW + $\Delta N_{\mathrm{eff}}$ + PBH'],
    filled = True,
    line_args=[
        {'lw':1.5,'ls':'-', 'color':'g'},
        {'lw':1.5,'ls':'-', 'color':'r'},
        {'lw':1.5,'ls':'-', 'color':'b'},
        ],
    title_limit=1,
    markers={'dNeff': 0.175},
    param_limits = {'dNeff': [0, 0.5]}
    )

g.export(ResultFile)
print(ResultFile)
