# compare my merger rate with 2012.02786
from src.merger_module_2 import *

reload = 1
f = np.logspace(-4, 0, 50)

LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
from joblib import Parallel, delayed

f0, r0 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2012.02786.fig1.black.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

def model(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.6, z = 0, sbh_width = 10, nm = 200, show_status = 0)
    PyLab.SaySomething()
    return r

if reload:
    t1 = PyLab.TimeNow()
    r = Parallel(n_jobs = 12)(delayed(model)(x) for x in f)
    PyLab.Timer(t1)
    np.savez('data/merger_rate.npz', r = r)

r = np.load('data/merger_rate.npz')['r']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(f0, r0, '-k', linewidth=LineWidth, label = '2012.02786')
plt.loglog(f, r, '-r', linewidth=LineWidth, label = 'my codes')

plt.xlabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[\mathrm{Gpc^{-3} yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
