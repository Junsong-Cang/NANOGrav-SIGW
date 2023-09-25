# compare my merger rate with 1812 paper
from src.merger import *

reload = 1
f = np.logspace(-4, 0, 50)

LineWidth = 2
FontSize = 18
nm = 100
sbh_width = 5

import matplotlib.pyplot as plt
from joblib import Parallel, delayed

f1, R1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_solid.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f2, R2 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1812.01930.fig3.black_dashed.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

def model_1(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 0.6, z = 0, sbh_width = sbh_width, nm = nm, show_status = 0, mf_model = 0, S1_method = 0)
    PyLab.SaySomething()
    return r

def model_2(fbh):
    r = Merger_Rate(fbh = fbh, mc = 20, sbh = 2, z = 0, sbh_width = sbh_width, nm = nm, show_status = 0, mf_model = 0, S1_method = 0)
    PyLab.SaySomething()
    return r

if reload:
    t1 = PyLab.TimeNow()
    r1 = Parallel(n_jobs = 12)(delayed(model_1)(x) for x in f)
    r2 = Parallel(n_jobs = 12)(delayed(model_2)(x) for x in f)
    PyLab.Timer(t1)
    np.savez('tmp.npz', r1 = r1, r2 = r2)

r = np.load('tmp.npz')
r1 = r['r1']
r2 = r['r2']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(f1, R1, '-k', linewidth=LineWidth, label = '$\sigma = 0.6$')
plt.loglog(f2, R2, '-b', linewidth=LineWidth, label = '$\sigma = 2$')
plt.loglog(f, r1, '--k', linewidth=LineWidth, label = '$\sigma = 0.6$, my codes')
plt.loglog(f, r2, '--b', linewidth=LineWidth, label = '$\sigma = 2$, my codes')

plt.xlabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[\mathrm{Gpc^{-3} yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower right')
#plt.xlim([-1,7])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
