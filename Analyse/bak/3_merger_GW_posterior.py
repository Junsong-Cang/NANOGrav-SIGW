from src.merger import *

reload = 0
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_GW_posteriors.pdf'
v = np.logspace(-3, 10, 80)
LineWidth = 2
FontSize = 16
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

r = np.load('data/3_merger_GW_posterior_np_zp.npz')
r4 = np.load('data/4_merger_GW_posterior_no_zp.npz')
r1 = r['r1']
r2 = r['r2']
r41 = r4['r1']
r42 = r4['r2']

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

f1, f2, f3, f4, f5, f6, f7, f8, f9 = 1e6*f1, 1e6*f2, 1e6*f3, 1e6*f4, 1e6*f5, 1e6*f6, 1e6*f7, 1e6*f8, 1e6*f9

g1 = h2_OmGW(h = h1, f = f1)
g2 = h2_OmGW(h = h2, f = f2)
g3 = h2_OmGW(h = h3, f = f3)
g4 = h2_OmGW(h = h4, f = f4)
g5 = h2_OmGW(h = h5, f = f5)
g6 = h2_OmGW(h = h6, f = f6)
g7 = h2_OmGW(h = h7, f = f7)
g8 = h2_OmGW(h = h8, f = f8)
g9 = h2_OmGW(h = h9, f = f9)
g10 = h2_OmGW(h = h10, f = f10)
g11 = h2_OmGW(h = h11, f = f11)
g12 = h2_OmGW(h = h12, f = f12)
g13 = h2_OmGW(h = h13, f = f13)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

Top = 1e20
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

#plt.fill_between(f1, g1, Top, color = 'r', alpha=0.3, label='Experimental Reach')
plt.fill_between(f1, g1, Top, color = 'r', alpha=0.3)
plt.fill_between(f2, g2, Top, color = 'r', alpha=0.3)
plt.fill_between(f3, g3, Top, color = 'r', alpha=0.3)
plt.fill_between(f4, g4, Top, color = 'r', alpha=0.3)
plt.fill_between(f5, g5, Top, color = 'r', alpha=0.3)
plt.fill_between(f6, g6, Top, color = 'r', alpha=0.3)
plt.fill_between(f7, g7, Top, color = 'r', alpha=0.3)
#plt.fill_between(f8, g8, Top, color = 'r', alpha=0.3)
#plt.fill_between(f9, g9, Top, color = 'r', alpha=0.3)
plt.fill_between(f10, g10, Top, color = 'b', alpha=0.3, linestyle = 'dashed')
plt.fill_between(f13, g13, Top, color = 'r', alpha=0.3, linestyle = 'dashed')

plt.fill_between(f11, g11, Top, color = 'g', alpha=0.55, linestyle = 'dashed')
plt.fill_between(f12, g12, Top, color = 'm', alpha=0.5, linestyle = 'dashed')

plt.fill_between(v, r42[3,:], r42[4,:], color = 'g', alpha=0.5, label = 'No S2')
plt.fill_between(v, r41[3,:], r41[4,:],color = 'b',alpha=0.6,label = 'With S2')
plt.fill_between(v, r42[4,:], r2[4,:], color = 'grey',alpha=0.5, linestyle = 'dashed')

# Also plot best-fit
GW_1d = np.load('data/4_merger_GW.npz')
v4 = GW_1d['v']
r4 = GW_1d['r']

plt.plot(v4, r4[0,:], 'k', linewidth = LineWidth)
plt.plot(v4, r4[1,:], '--k', linewidth = LineWidth)
plt.plot(v4, r4[2,:], color = 'k', linestyle = '-.', linewidth = LineWidth)

'''
plt.plot(v2, r2[0,:], 'k', linewidth = LineWidth, label = 'Best-fit')
plt.plot(v2, r2[1,:], '--k', linewidth = LineWidth, label = '$\sigma = 0$')
plt.plot(v2, r2[2,:], 'r', linewidth = LineWidth, label = '$\sigma = 1$')
'''

# Plot Neff
Neff_GW = 2.1e-6*np.ones(len(v))
plt.plot(v, Neff_GW, color = 'grey', linestyle = ':', linewidth = LineWidth)

# Now add texts
plt.text(3e-3, 2e-12, "$f_{\mathrm{bh}} \in [10^{-20}, 1]$", size=FontSize, rotation = 25,color='k')
plt.text(1e6, 3e-6, "$\Delta N_{\mathrm{eff}} < 0.175$", size=FontSize, rotation = 0, color='k')
#plt.text(6e-2, 5e-6, "Levitated Sensors", size=FontSize, rotation = 90, color='grey')
#plt.text(1.5e0, 4e-3, "BAW", size=FontSize, rotation = 90, color='grey')
#plt.text(2.5e1, 2.5e-2, "EDGES", size=FontSize, rotation = 90, color='grey')
#plt.text(2.3e2, 4e-2, "ADMX", size=FontSize, rotation = 90, color='grey')
#plt.text(8e2, 1e-3, "SQMS", size=FontSize, rotation = 90, color='grey')
#plt.text(3.5e3, 1e-3, "ARCADE", size=FontSize, rotation = 90, color='grey')
plt.text(2e0, 1e-4, "ET", size=FontSize/1.3, rotation = 90, color='grey')
plt.text(1.6e1, 3e-6, "A+", size=FontSize/1.3, rotation = 0, color='k')
plt.text(2e1, 1e-2, "aLIGO", size = FontSize/1.3, rotation = 0, color='grey')
plt.text(1.5e2, 1e-7, "CE", size=FontSize/1.3, rotation = 60, color='grey')

# plt.legend(fontsize = FontSize/1.5, loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(fontsize=FontSize,loc = 'upper right')

plt.xscale('log')
plt.yscale('log')

# plt.title('Posterior for merger GW',fontsize=FontSize)
plt.xlabel('$f$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('${\mathrm{d}}\Omega_{\mathrm{GW}}/{\mathrm{d}}\ln f$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1e-3, 1e10])
plt.ylim([1e-13, 1e2])

plt.tight_layout()
plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)
