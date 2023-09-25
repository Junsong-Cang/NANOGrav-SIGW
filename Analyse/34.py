Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_Merger_Rate_1D'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'

ResultFile='/Users/cangtao/Desktop/tmp.pdf'

# ---- Initialise ----
import getdist, os
from getdist import plots
import matplotlib.pyplot as plt

# ---- Getdist Plot ----
plt.rcParams.update({'font.family':'Times'})

samples_1 = getdist.mcsamples.loadMCSamples(Root_1)
samples_2 = getdist.mcsamples.loadMCSamples(Root_2)

p = samples_1.getParams()
g = plots.getSubplotPlotter(subplot_size = 3)
g.settings.axes_fontsize=14
g.settings.title_limit_fontsize = 14

g.triangle_plot(
    [samples_1, samples_2],
    ['LgMc', 'LgFbh', 'Sbh'],
    width_inch=12,
    contour_colors=['g', 'b'],
    legend_labels=['$f_{\mathrm{bh}} \in [10^{-20}, 1]$', '$f_{\mathrm{bh}} \in [10^{-20}, 0.1]$'],
    filled = True,
    line_args=[
        {'lw':1.5,'ls':'-', 'color':'g'},
        {'lw':1.5,'ls':'-', 'color':'b'},
        ],
    title_limit=1,
    markers={'dNeff': 0.175},
    param_limits = {'dNeff': [0, 0.5]}
    )

g.export(ResultFile)
print(ResultFile)
