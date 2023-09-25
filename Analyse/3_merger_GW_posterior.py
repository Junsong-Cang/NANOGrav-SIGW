from src.merger import *

reload = 0
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_GW_posteriors.pdf'
v = np.logspace(-3, 10, 80)
LineWidth = 2
FontSize = 16
ncpu = 12
nm = 50
nz = 50

import matplotlib.pyplot as plt
import os

LogFile_0 = FileRoot + 'GW_0.txt'
LogFile_1 = FileRoot + 'GW_1.txt'

def model_0(theta):
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
        nm = nm,
        nz = nz,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething('tmp.txt')

    # Print results to a file
    
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    nv = len(v)
    for idx in np.arange(0, nv):
        Str = Str + '  {0:.5E}'.format(g[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    return g

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
        nm = nm,
        nz = nz,
        show_status = 0,
        Use_interp = 1,
        S1_method = S1_method,
        S_Tab_Len = 200,
        Use_S2 = Use_S2,
        Precision = 1e-2,
        ncpu = 1)
    PyLab.SaySomething('tmp.txt')
    
    # Print results to a file
    
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    nv = len(v)
    for idx in np.arange(0, nv):
        Str = Str + '  {0:.5E}'.format(g[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()
    
    return g

if reload:
    try:
        os.remove(LogFile_0)
    except:
        pass
    try:
        os.remove(LogFile_1)
    except:
        pass
    t1 = PyLab.TimeNow()
    r0 = PyLab.mcmc_derived_stat(model_function = model_0, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    r1 = PyLab.mcmc_derived_stat(model_function = model_1, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    # np.savez('data/3_merger_GW_posterior.npz', r0 = r0, r1 = r1)
    PyLab.Timer(t1)

r = np.load('data/3_merger_GW_posterior.npz')
r4 = np.load('data/4_merger_GW_posterior_no_zp.npz')
r1 = r['r1']
r2 = r['r2']
r41 = r4['r1']
r42 = r4['r2']

Curve_Path = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

f1, g1 = PyLab.Read_Curve(
    File = Curve_Path + '2109_11376_fig1_LISA.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f2, g2 = PyLab.Read_Curve(
    File = Curve_Path + '2109_11376_fig1_ET.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f3, g3 = PyLab.Read_Curve(
    File = Curve_Path + '2109_11376_fig1_CE.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f4, g4 = PyLab.Read_Curve(
    File = Curve_Path + '2109_11376_fig1_aLIGO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f5, g5 = PyLab.Read_Curve(
    File = Curve_Path + '2109_11376_fig1_A+.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f6, g6 = PyLab.Read_Curve(
    File = Curve_Path + '1310_5300_fig11_BBO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f7, g7 = PyLab.Read_Curve(
    File = Curve_Path + '2109_01398_fig1_DECIGO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)
g7 = g7/0.6766**2


# Get plot now
Top = 1e20
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

# Confidence regions
'''
plt.fill_between(v, r42[3,:], r42[4,:], color = 'g', alpha=0.5, label = 'No S2')
plt.fill_between(v, r41[3,:], r41[4,:],color = 'b',alpha = 0.6,label = 'With S2')
plt.fill_between(v, r42[4,:], r2[4,:], color = 'grey',alpha=0.5, linestyle = 'dashed')
'''

idx = 6
v_ = v[::idx]

plt.fill_between(v_, r42[3,::idx], r42[4,::idx], color = 'g', alpha=0.5, label = 'No S2')
plt.fill_between(v_, r41[3,::idx], r41[4,::idx],color = 'b',alpha = 0.6,label = 'With S2')
plt.fill_between(v_, r42[4,::idx], r2[4,::idx], color = 'grey',alpha=0.5, linestyle = 'dashed')

# Best-fit
GW_1d = np.load('data/4_merger_GW.npz')
v4 = GW_1d['v']
r4 = GW_1d['r']
plt.plot(v4, r4[0,:], 'k', linewidth = LineWidth)
plt.plot(v4, r4[1,:], '--k', linewidth = LineWidth)
plt.plot(v4, r4[2,:], color = 'k', linestyle = '-.', linewidth = LineWidth)

# Experimental reach
plt.fill_between(f7, g7, Top, color = 'r', alpha = 0.4, linestyle = 'dashed')
plt.fill_between(f1, g1, Top, color = 'b', alpha = 0.3, linestyle = 'dashed')
plt.fill_between(f2, g2, Top, color = 'r', alpha = 0.4, linestyle = 'dashed')
plt.fill_between(f3, g3, Top, color = 'y', alpha = 0.4, linestyle = 'dashed')
plt.fill_between(f4, g4, Top, color = 'g', alpha = 0.4, linestyle = 'dashed')
#plt.fill_between(f5, g5, Top, color = 'm', alpha = 0.4, linestyle = 'dashed')
plt.fill_between(f6, g6, Top, color = 'r', alpha = 0.4, linestyle = 'dashed')

# Plot Neff
Neff_GW = 2.1e-6*np.ones(len(v))
plt.plot(v, Neff_GW, color = 'grey', linestyle = ':', linewidth = LineWidth)

# Now add texts
#plt.text(4e3, 1e-8, "$f_{\mathrm{bh}} \in [10^{-20}, 1]$", size=FontSize/1.3, rotation = 25,color='k')
plt.text(4e3, 1e-8, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = 25,color='k')

plt.text(2e3, 3e-6, "$\Delta N_{\mathrm{eff}} < 0.175$", size=FontSize/1.3, rotation = 0, color='k')
plt.text(2e-2, 1e-14, "DECIGO", size=FontSize/1.3, rotation = 0, color='k')
plt.text(1.2e-3, 1e-11, "LISA", size=FontSize/1.3, rotation = 0, color='k')
plt.text(1.2e1, 3e-12, "ET", size=FontSize/1.3, rotation = 0, color='k')
plt.text(6e1, 5e-13, "CE", size=FontSize/1.3, rotation = 0, color='k')
plt.text(1.5e1, 2e-9, "aLIGO", size=FontSize/1.3, rotation = 0, color='k')
plt.text(4e-2, 1e-16, "BBO", size=FontSize/1.3, rotation = 0, color='k')

plt.legend(fontsize=FontSize,loc = 'upper right')

plt.xscale('log')
plt.yscale('log')

# plt.title('Posterior for merger GW',fontsize=FontSize)
plt.xlabel('$f$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1e-3, 1e10])
plt.ylim([1e-17, 1e-1])

plt.tight_layout()
plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)

print(len(v_))
