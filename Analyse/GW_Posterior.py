from PyLab import *

Use_Log = 0
reload = 0
FileRoot_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
FileRoot_4 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_GW_posteriors.pdf'
v = np.logspace(-3, 10, 80)

LineWidth = 2
FontSize = 16

Tab_30 = np.loadtxt(FileRoot_3 + 'GW_0.txt')
Tab_31 = np.loadtxt(FileRoot_3 + 'GW_1.txt')
Tab_40 = np.loadtxt(FileRoot_4 + 'GW_0.txt')
Tab_41 = np.loadtxt(FileRoot_4 + 'GW_1.txt')
nv = len(v)

def kernel(theta, type):
    '''
    Read GW from computed tables
    '''

    lf, lm, sbh = theta[5], theta[6], theta[7]

    if type == 30:
        Tab = Tab_30
    elif type == 31:
        Tab = Tab_31
    elif type == 40:
        Tab = Tab_40
    else:
        Tab = Tab_41
    
    lf_vec = Tab[:,0]
    lm_vec = Tab[:,1]
    sbh_vec = Tab[:,2]
    dist = (lf - lf_vec)**2 + (lm - lm_vec)**2 + (sbh - sbh_vec)**2
    idx = np.argmin(dist)
    r = Tab[idx, 3 : 3+nv]
    if Use_Log:
        r = np.log10(r)
    return r

def model_30(theta):
    return(kernel(theta, 30))

def model_31(theta):
    return(kernel(theta, 31))

def model_40(theta):
    return(kernel(theta, 40))

def model_41(theta):
    return(kernel(theta, 41))

if reload:
    # Get posterior
    P30 = mcmc_derived_stat(model_function = model_30, FileRoot = FileRoot_3, ncpu = 10, print_status=1)
    P31 = mcmc_derived_stat(model_function = model_31, FileRoot = FileRoot_3, ncpu = 10, print_status=1)
    P40 = mcmc_derived_stat(model_function = model_40, FileRoot = FileRoot_4, ncpu = 10, print_status=1)
    P41 = mcmc_derived_stat(model_function = model_41, FileRoot = FileRoot_4, ncpu = 10, print_status=1)
    np.savez('data/GW_Posterior.npz', P30 = P30, P31 = P31, P40 = P40, P41 = P41)

R = np.load('data/GW_Posterior.npz')
if Use_Log:
    P30 = 10**R['P30']
    P31 = 10**R['P31']
    P40 = 10**R['P40']
    P41 = 10**R['P41']
else:
    P30 = R['P30']
    P31 = R['P31']
    P40 = R['P40']
    P41 = R['P41']

Curve_Path = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

f1, g1 = Read_Curve(
    File = Curve_Path + '2109_11376_fig1_LISA.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f2, g2 = Read_Curve(
    File = Curve_Path + '2109_11376_fig1_ET.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f3, g3 = Read_Curve(
    File = Curve_Path + '2109_11376_fig1_CE.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f4, g4 = Read_Curve(
    File = Curve_Path + '2109_11376_fig1_aLIGO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f5, g5 = Read_Curve(
    File = Curve_Path + '2109_11376_fig1_A+.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f6, g6 = Read_Curve(
    File = Curve_Path + '1310_5300_fig11_BBO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)

f7, g7 = Read_Curve(
    File = Curve_Path + '2109_01398_fig1_DECIGO.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)
g7 = g7/0.6766**2

# Now the strain

def h2_OmGW(h, f):
    H0 = 2.192695336552484e-18
    C = 4 * np.pi**2 /(3 * H0**2)
    OmGW = C * f**2 * h**2
    return OmGW

f8, h8 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.LSD.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f9, h9 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.DMR.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f10, h10 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.BAW.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f11, h11 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.EDGES.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f12, h12 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.ADMX.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f13, h13 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.SQMS.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f14, h14 = Read_Curve(
    File = Curve_Path + '2205.02153.fig5.ARCADE.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)

f17, g17 = Read_Curve(
    File = Curve_Path + '2306.16219.fig4.HLV.txt',
    model = 1,
    Convert_x = 1,
    Convert_y = 1)
g17 = g17/0.6766**2

f18, g18 = Read_Curve(
    File = Curve_Path + '2306.16219.fig4.HLVK.txt',
    model = 3,
    Convert_x = 1,
    Convert_y = 1)
g18 = g18/0.6766**2

g8 = h2_OmGW(h = h8, f = f8)
g9 = h2_OmGW(h = h9, f = f9)
g10 = h2_OmGW(h = h10, f = f10)
g11 = h2_OmGW(h = h11, f = f11)
g12 = h2_OmGW(h = h12, f = f12)
g13 = h2_OmGW(h = h13, f = f13)
g14 = h2_OmGW(h = h14, f = f14)

# Get plot now
Top = 1e20
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

'''
plt.loglog(f8, h8, 'r')
plt.loglog(f9, h9, 'r')
plt.loglog(f10, h10, 'r')
plt.loglog(f11, h11, 'r')
plt.loglog(f12, h12, 'b', label = 'ADMX')
plt.loglog(f13, h13, 'k', label = 'SQMS')
plt.loglog(f14, h14, 'r')
plt.legend(fontsize=FontSize,loc = 'upper left')

plt.show()
'''

# Confidence regions

idx = 6
v_ = v[::idx]

plt.fill_between(v_, P40[3,::idx], P40[4,::idx], color = 'g', alpha=0.5, label = 'No $S_2$')
plt.fill_between(v_, P41[3,::idx], P41[4,::idx],color = 'b',alpha = 0.6,label = 'With $S_2$')
plt.fill_between(v_, P40[4,::idx], P30[4,::idx], color = 'grey',alpha=0.5, linestyle = 'dashed')

# Best-fit
GW_1d = np.load('data/4_merger_GW.npz')
v4 = GW_1d['v']
r4 = GW_1d['r']
plt.plot(v4, r4[0,:], 'k', linewidth = LineWidth)
plt.plot(v4, r4[1,:], '--k', linewidth = LineWidth)
plt.plot(v4, r4[2,:], color = 'k', linestyle = '-.', linewidth = LineWidth)

# Experimental reach
plt.fill_between(f7, g7, Top, color = 'r', alpha = 0.4, linestyle = 'dashed')
#plt.plot(f7, g7, '--r', linewidth = LineWidth)

plt.fill_between(f1, g1, Top, color = 'b', alpha = 0.3, linestyle = 'dashed')
plt.fill_between(f2, g2, Top, color = 'm', alpha = 0.3, linestyle = 'dashed')
plt.fill_between(f3, g3, Top, color = 'y', alpha = 0.3, linestyle = 'dashed')
plt.fill_between(f4, g4, Top, color = 'g', alpha = 0.45, linestyle = 'solid')
plt.fill_between(f6, g6, Top, color = 'r', alpha = 0.3, linestyle = 'dashed')
plt.fill_between(f18, g18, Top, color = 'c', alpha = 0.4, linestyle = 'dashed')
plt.fill_between(f17, g17, Top, color = 'b', alpha = 0.3, linestyle = 'solid')


# Now add texts
# plt.text(4e8, 1e-11, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = -60,color='k')
plt.text(6e8, 3e-12, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = -60,color='k')

#plt.text(2e3, 4e-6, "$\Delta N_{\mathrm{eff}} < 0.175$", size=FontSize/1.3, rotation = 0, color='k')
plt.text(2e-2, 1e-14, "DECIGO", size=FontSize/1.3, rotation = 0, color='k')
plt.text(1.2e-3, 1e-11, "LISA", size=FontSize/1.3, rotation = 0, color='k')
#plt.text(1.2e1, 3e-12, "ET", size=FontSize/1.3, rotation = 0, color='k')
plt.text(6e1, 5e-13, "CE", size=FontSize/1.3, rotation = 0, color='k')
plt.text(4e1, 5e-10, "aLIGO", size=FontSize/1.3, rotation = 50, color='k')

plt.text(2.3e2, 1e-9, "ET", size=FontSize/1.5, rotation = 50, color='k')

plt.text(4e-2, 1e-16, "BBO", size=FontSize/1.3, rotation = 0, color='k')
plt.text(2e1, 1e-7, "HLV", size=FontSize/1.3, rotation = 0, color='k')

plt.text(8e0, 2e-9, "HLVK", size=FontSize/1.5, rotation = -80, color='k')

plt.xscale('log')
plt.yscale('log')

# plt.legend(fontsize=FontSize,loc = 'upper right')
plt.xlabel('$f$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1e-3, 1e10])
plt.ylim([1e-17, 1e0])

# inset for MHz experiments

left, bottom, width, height = [0.773, 0.765, 0.2, 0.2]
# left, bottom, width, height = [0.76, 0.765, 0.2, 0.2]
left, bottom, width, height = [0.724, 0.715, 0.25, 0.25]

# Strain converted reach
ax_MHz = fig.add_axes([left, bottom, width, height])
ax_MHz.fill_between(f8, g8, Top, color = 'r', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f9, g9, Top, color = 'g', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f10, g10, Top, color = 'm', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f11, g11, Top, color = 'm', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f12, g12, Top, color = 'b', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f13, g13, Top, color = 'g', alpha = 0.4, linestyle = 'solid')
ax_MHz.fill_between(f14, g14, Top, color = 'y', alpha = 0.4, linestyle = 'solid')
ax_MHz.set_xscale('log')
ax_MHz.set_yscale('log')
ax_MHz.set_xlim(1e4, 8e9)
ax_MHz.set_ylim(1e2, 4e12)

plt.text(4e4, 7e3, "Levitated Sensors", size=FontSize/2, rotation = 90, color='k')
plt.text(1e6, 2e8, "BAW", size=FontSize/2, rotation = 40, color='k')
plt.text(3.5e7, 3e9, "EDGES", size=FontSize/2, rotation = 90, color='k')
plt.text(2.5e8, 7e9, "ADMX", size=FontSize/2, rotation = 90, color='k')
plt.text(1e9, 2e9, "SQMS", size=FontSize/2.3, rotation = 90, color='k')
plt.text(3e9, 0.6e8, "ARCADE", size=FontSize/2.0, rotation = 90, color='k')
plt.text(4e5, 1.6e6, "DMRadio", size=FontSize/2.0, rotation = 0, color='k')

plt.tight_layout()
plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)

print(len(v_))

