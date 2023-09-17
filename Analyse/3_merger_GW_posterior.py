from src.merger import *

reload = 0
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_GW_posteriors.pdf'
v = np.logspace(-3, 10, 80)
LineWidth = 2
FontSize = 18
ncpu = 12

import matplotlib.pyplot as plt

def model_1(theta):
    Use_S2 = 1
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    g = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        Fast = 0,
        mf_model = 0, 
        sbh_width = 6, 
        nm = 50,
        nz = 50,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething()
    return g

def model_2(theta):
    Use_S2 = 0
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    g = Get_dOmGW_dlnv(
        fbh = fbh,
        mc = mc,
        sbh = sbh, 
        v = v,
        Fast = 0,
        mf_model = 0, 
        sbh_width = 6, 
        nm = 50,
        nz = 50,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething()
    return g

if reload:
    t1 = PyLab.TimeNow()
    r1 = PyLab.mcmc_derived_stat(model_function = model_1, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    r2 = PyLab.mcmc_derived_stat(model_function = model_2, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    np.savez('data/3_merger_GW_posterior.npz', r1 = r1, r2 = r2)
    PyLab.Timer(t1)

r = np.load('data/3_merger_GW_posterior.npz')
r1 = r['r1']
r2 = r['r2']
print(np.shape(r1))
Curve_Path = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'
f1, h1 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.levitated_sensors.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f2, h2 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.BAW.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f3, h3 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.holometer.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f4, h4 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.EDGES.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f5, h5 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.ADMX.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f6, h6 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.SQMS.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f7, h7 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.ARCADE.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f8, h8 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.UHF_GW_Landscape.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f9, h9 = PyLab.Read_Curve(
    File = Curve_Path + '2202.00695.fig1.DMRadio8.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

g1 = h2_OmGW(h = h1, f = f1)
g2 = h2_OmGW(h = h2, f = f2)
g3 = h2_OmGW(h = h3, f = f3)
g4 = h2_OmGW(h = h4, f = f4)
g5 = h2_OmGW(h = h5, f = f5)
g6 = h2_OmGW(h = h6, f = f6)
g7 = h2_OmGW(h = h7, f = f7)
g8 = h2_OmGW(h = h8, f = f8)
g9 = h2_OmGW(h = h9, f = f9)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.ylim([1e-13, 1e1])

Top = 100.0
'''
plt.fill_between(f1, g1, Top,color = 'b', alpha=0.3, label='Levitated Sensors')
plt.fill_between(f2, g2, Top,color = 'r', alpha=0.3, label='BAW')
plt.fill_between(f3, g3, Top,color = 'g', alpha=0.3, label='holometer')
plt.fill_between(f4, g4, Top,color = 'y', alpha=0.3, label='EDGES')
plt.fill_between(f5, g5, Top,color = 'm', alpha=0.3, label='ADMX')
plt.fill_between(f6, g6, Top,color = 'c', alpha=0.3, label='SQMS')
plt.fill_between(f7, g7, Top,color = 'grey', alpha=0.3, label='ARCADE')
plt.fill_between(f8, g8, Top,color = 'brown', alpha=0.3, label='UHF GW Landscape')
plt.fill_between(f9, g9, Top,color = 'purple', alpha=0.3, label='DMRadio8')
'''

plt.fill_between(f1, g1, Top, color = 'grey', alpha=0.3, label='Experimental Reach')
plt.fill_between(f2, g2, Top, color = 'grey', alpha=0.3)
plt.fill_between(f3, g3, Top, color = 'grey', alpha=0.3)
plt.fill_between(f4, g4, Top, color = 'grey', alpha=0.3)
plt.fill_between(f5, g5, Top, color = 'grey', alpha=0.3)
plt.fill_between(f6, g6, Top, color = 'grey', alpha=0.3)
plt.fill_between(f7, g7, Top, color = 'grey', alpha=0.3)
plt.fill_between(f8, g8, Top, color = 'grey', alpha=0.3)
plt.fill_between(f9, g9, Top, color = 'grey', alpha=0.3)

plt.fill_between(v, r2[3,:], r2[4,:],color = 'g',alpha=0.5,label = '$2 \sigma$, no S2')
# plt.plot(v, r2[0,:], 'b', linewidth = LineWidth)

plt.fill_between(v, r1[3,:], r1[4,:],color = 'b',alpha=0.6,label = '$2 \sigma$, with S2')

# Also plot best-fit
v2 = np.logspace(0, 10, 100)
r2 = np.load('data/3_merger.npz')['r']
plt.plot(v2, r2[0,:], 'k', linewidth = LineWidth, label = 'best-fit')

plt.xscale('log')
plt.yscale('log')

# plt.title('Posterior for merger GW',fontsize=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('${\mathrm{d}}\Omega_{\mathrm{GW}}/{\mathrm{d}}\ln f$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1e-3, 1e10])
plt.ylim([1e-13, 1e1])


plt.tight_layout()
plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)
