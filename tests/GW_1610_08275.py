from src.merger import *

reload = 1
v = np.logspace(0, 7, 100)

LineWidth = 2
FontSize = 15
nm = 50
nz = 100
zmax = 10000
Use_S2 = 0
S1_method = 0

import matplotlib.pyplot as plt

v1, g1 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1610.08275.fig2.red.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

v2, g2 = PyLab.Read_Curve(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/1610.08275.fig2.blue.txt',
    nx = 100,
    model = 2,
    Convert_x = 1,
    Convert_y = 1)

if reload:
    t1 = PyLab.TimeNow()
    g1_ = Get_dOmGW_dlnv(
        fbh = 0.05,
        mc = 1,
        sbh = 0.02, 
        v = v,
        mf_model = 2,
        sbh_width = 10, 
        nm = nm,
        nz = nz,
        zmax = zmax,
        show_status = 1,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        ncpu = 12)
    g2_ = Get_dOmGW_dlnv(
        fbh = 0.06,
        mc = 0.1,
        sbh = 0.02, 
        v = v,
        mf_model = 2,
        sbh_width = 10, 
        nm = nm,
        nz = nz,
        zmax = zmax,
        show_status = 1,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        ncpu = 12)
    np.savez('tmp.npz', g1_ = g1_, g2_ = g2_)
    PyLab.Timer(t1)

r = np.load('tmp.npz')
g1_ = r['g1_']
g2_ = r['g2_']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()
plt.loglog(v1, g1, 'k', linewidth=LineWidth, label = '1610.08275')
plt.loglog(v, g1_, '--k', linewidth=LineWidth, label = 'my codes')
plt.loglog(v2, g2, 'r', linewidth=LineWidth)
plt.loglog(v, g2_, '--r', linewidth=LineWidth)

plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower left')
plt.xlim([1, 1e6])
plt.ylim([1e-12, 1e-6])

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
