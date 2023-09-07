Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/0_NG15/0_NG15_'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'

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
    # ['LgK', 'LgA', 'S', 'dNeff'],
    ['LgK', 'LgA', 'S'],
    width_inch=12,
    contour_colors=['g', 'r', 'b'],
    legend_labels=['NG15', 'NG15 + other'],
    filled = True,
    line_args=[
        {'lw':1.5,'ls':'-', 'color':'g'},
        {'lw':1.5,'ls':'-', 'color':'r'},
        ],
    title_limit=1,
    markers={'dNeff': 0.175},
    param_limits = {'dNeff': [0, 0.5]}
    )

g.export(ResultFile)
print(ResultFile)
