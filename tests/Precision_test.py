from src.merger_module import *
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

reload = 0

fbh = 1e-3
mc = 30
sbh = 1

v = np.logspace(-6, 4, 40)

ncpu = 12
LineWidth = 2
FontSize = 18

if reload:
    t1 = PyLab.TimeNow()
    g = Get_dOmGW_dlnv(v = v, show_status = 1, ncpu = ncpu, sbh = sbh, fbh = fbh, mc = mc, nm = 300, nz = 300, mf_model = 0)
    PyLab.Timer(t1)
    np.savez('tmp.npz', g = g)

g = np.load('tmp.npz')['g']

t1 = PyLab.TimeNow()
v2 = np.logspace(-6, 4, 15)
g2 = Get_dOmGW_dlnv(
    v = v2, 
    mf_model = 0,
    sbh = sbh, 
    fbh = fbh, 
    mc = mc,
    show_status = 1,
    ncpu = 1,
    nm = 30,
    nz = 30,
    sbh_width = 5,
    zmax = 25)
PyLab.Timer(t1)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

# fig.set_size_inches(8, 8)
plt.loglog(v, g, 'k', linewidth=LineWidth, label = 'HiRes')
plt.loglog(v2, g2, '--r', linewidth=LineWidth, label = 'LowRes')

plt.xlabel('$v$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('GW',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower left')

#plt.xlim([1e-6, 1e4])
#plt.ylim([1e-14,1e-5])
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
plt.show()
