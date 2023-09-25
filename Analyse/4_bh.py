FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/PBH_posteriors.pdf'

# ---- Initialise ----
import getdist, os
from getdist import plots
import matplotlib.pyplot as plt
import numpy as np

# ---- Getdist Plot ----
plt.rcParams.update({'font.family':'Times New Roman'})
samples = getdist.mcsamples.loadMCSamples(FileRoot)
p = samples.getParams()
g = plots.getSubplotPlotter(subplot_size = 3)
g.settings.axes_fontsize=15
g.settings.title_limit_fontsize = 15
g.settings.lab_fontsize =15

g.triangle_plot(
    samples,
    ['LgMc', 'LgFbh', 'Sbh'],
    width_inch=12,
    contour_colors=['blue'],
    filled = True,
    line_args=[{'lw':1.5,'ls':'-', 'color':'k'}],
    title_limit=1,
    param_limits = {
        'LgFbh': [-20, -1]
        })

g.export(ResultFile)
print(ResultFile)

# Now get some GW plots

best_fit = samples.getLikeStats()
print(best_fit)

fbh = np.array([1e-2, 1e-3, 1e-1])
mc = np.array([1e-3, 1e-4, 1e-2])
sbh = np.array([0.3, 0.1, 1])

print('    fbh         mc       sigma')
print("{0:.4E}".format(fbh[0]), "{0:.4E}".format(mc[0]), "  {0:.4f}".format(sbh[0]))

for fid in np.arange(1,3):
    for mid in np.arange(1,3):
        for sid in np.arange(1,3):
            f_ = fbh[fid]
            m = mc[mid]
            s = sbh[sid]
            print("{0:.4E}".format(f_), "{0:.4E}".format(m), "  {0:.4f}".format(s))

print('Plot saved to :')
print(ResultFile)
