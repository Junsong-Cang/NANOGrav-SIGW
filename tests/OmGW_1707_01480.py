from src.merger_module import *

reload = 1
v = np.logspace(-6, 4, 100)
LineWidth = 2
FontSize = 15
nm = 50
nz = 50
model = 0

import matplotlib.pyplot as plt
from joblib import Parallel, delayed

v1, g1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1707.01480.fig2.red_solid_top.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

v2, g2 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1707.01480.fig2.blue_solid_top.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

v3, g3 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1707.01480.fig2.red_solid_lower.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

v4, g4 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1707.01480.fig2.blue_solid_lower.txt',
    nx = 100,
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

if reload:
    t1 = PyLab.TimeNow()
    g1_ = Get_dOmGW_dlnv(mc = 30, v = v, mf_model = model, ncpu = 12, nm = nm, nz = nz, fbh = 1e-2, sbh = 1)
    print('g1_ done')
    g2_ = Get_dOmGW_dlnv(mc = 30, v = v, mf_model = model, ncpu = 12, nm = nm, nz = nz, fbh = 1e-2, sbh = 0.05)
    print('g2_ done')
    g3_ = Get_dOmGW_dlnv(mc = 30, v = v, mf_model = model, ncpu = 12, nm = nm, nz = nz, fbh = 1e-3, sbh = 1)
    print('g3_ done')
    g4_ = Get_dOmGW_dlnv(mc = 30, v = v, mf_model = model, ncpu = 12,nm = nm, nz = nz,  fbh = 1e-3, sbh = 0.05)
    print('g4_ done')
    np.savez('data/OmGW_1707_01480.npz', g1_ = g1_, g2_ = g2_, g3_ = g3_, g4_ = g4_)
    PyLab.Timer(t1)

r = np.load('data/OmGW_1707_01480.npz')
g1_ = r['g1_']
g2_ = r['g2_']
g3_ = r['g3_']
g4_ = r['g4_']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()
plt.loglog(v1, g1, 'k', linewidth=LineWidth, label = '1707.01480')
plt.loglog(v, g1_, '--k', linewidth=LineWidth, label = 'my codes')
plt.loglog(v2, g2, 'r', linewidth=LineWidth)
plt.loglog(v, g2_, '--r', linewidth=LineWidth)
plt.loglog(v3, g3, 'b', linewidth=LineWidth)
plt.loglog(v, g3_, '--b', linewidth=LineWidth)
plt.loglog(v4, g4, 'g', linewidth=LineWidth)
plt.loglog(v, g4_, '--g', linewidth=LineWidth)

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper left')
plt.xlim([1e-6, 1e4])
plt.ylim([1e-14, 1e-7])

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
