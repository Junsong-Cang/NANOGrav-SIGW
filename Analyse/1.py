FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'
ResultFile='/Users/cangtao/Desktop/tmp.png'

# ---- Initialise ----
import getdist, os
from getdist import plots
import matplotlib.pyplot as plt

# ---- Getdist Plot ----
plt.rcParams.update({'font.family':'Times'})
samples = getdist.mcsamples.loadMCSamples(FileRoot)
p = samples.getParams()
g = plots.getSubplotPlotter(subplot_size = 3)
g.settings.axes_fontsize=14
g.settings.title_limit_fontsize = 14

g.triangle_plot(
    samples,
    # ['LgF', 'LgA', 'S'],
    ['LgK', 'LgA', 'S'],
    width_inch=12,
    contour_colors=['blue'],
    filled = True,
    line_args=[{'lw':1.5,'ls':'-', 'color':'k'}],
    title_limit=1,
    #markers={'LgA' : 6.77, 'LgK':6.63, 'S':-0.04}
    markers={'dNeff': 0.175}
    )

g.export(ResultFile)
print(ResultFile)
