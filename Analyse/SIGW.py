from src.main import *

reload = 1
nv = 100
v = np.logspace(-9, -7.5, nv)
v = np.logspace(-9, -3, nv)

LineWidth = 2
FontSize = 15

Root_1 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/1_GW/1_GW_'
Root_2 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/2_GW_Neff/2_GW_Neff_'
Root_4 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/SIGW_posteriors.pdf'
h2 = 0.6766**2

def model(theta):
    LgFc, LgA, S = theta[0], theta[1], theta[2]
    r = dOmGWh2_dlnk(
        f = v,
        fc = 10**LgFc,
        k = v,
        kc = 1.0,
        A = 10**LgA,
        Sigma = S,
        Use_Freq_Domain = True,
        Use_Today_Value = True,
        Use_Delta_Approx = False,
        Use_Fast_Mode = True
    )
    return np.log10(r/h2)

def Get_Errors(GW_Data):
    f = GW_Data['freq']
    middle = GW_Data['LgGW_middle']
    top = GW_Data['LgGW_top']
    low = GW_Data['LgGW_lower']

    mean = 10**middle/h2
    sl = mean - 10**low/h2
    st = 10**top/h2 - mean
    
    r = {'f' : f, 'mean' : mean, 'sl' : sl, 'st' : st}
    return r

if reload:
    t1 = TimeNow()
    P1 = mcmc_derived_stat(model_function = model, FileRoot = Root_1, ncpu = 11, print_status=1)
    P2 = mcmc_derived_stat(model_function = model, FileRoot = Root_2, ncpu = 11, print_status=1)
    P4 = mcmc_derived_stat(model_function = model, FileRoot = Root_4, ncpu = 11, print_status=1)    
    Timer(t1)
    np.savez('data/SIGW_Posterior.npz', P1 = P1, P2 = P2, P4 = P4)
    
r = np.load('data/SIGW_Posterior.npz')
P1 = r['P1']
P2 = r['P2']
P4 = r['P4']

# Get plot now
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

# Confidence regions
plt.fill_between(v, 10**P1[1,:], 10**P1[2,:], color = 'g', alpha=0.5, label = 'GW')
plt.fill_between(v, 10**P2[1,:], 10**P2[2,:], color = 'r', alpha=0.6, label = 'GW + $\Delta N_{\mathrm{eff}}$')
plt.fill_between(v, 10**P4[1,:], 10**P4[2,:], color = 'b', alpha=0.5, label = 'GW + $\Delta N_{\mathrm{eff}}$ + PBH')

'''
NG = Get_Errors(NG15_conservative)
plt.errorbar(NG['f'], NG['mean'], [NG['sl'], NG['st']], color = 'k', linewidth=LineWidth, label = 'NG15', fmt='+')

IPTA = Get_Errors(IPTA_conservative)
plt.errorbar(IPTA['f'], IPTA['mean'], [IPTA['sl'], IPTA['st']], color = 'r', linewidth=LineWidth, label = 'IPTA', fmt='+')

PPTA = Get_Errors(PPTA)
plt.errorbar(PPTA['f'], PPTA['mean'], [PPTA['sl'], PPTA['st']], color = 'b', linewidth=LineWidth, label = 'PPTA', fmt='+')
'''

# now the nasty part


# NG15
NG_violin = Read_Curve_Pro(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2306_17136_fig1a_NG15.txt')
Violin = NG_violin
n = len(Violin)/2

for idx in np.arange(0, n):
    x1 = Violin[int(idx*2)]['x']
    y = Violin[int(idx*2)]['y']
    x2 = 2*x1[0] - x1
    
    y = 10**y/h2
    x1 = 10**x1
    x2 = 10**x2
    if idx == 0:
        plt.fill_betweenx(y, x1, x2, color = 'm', alpha=0.5, label = 'NG15')
    else:
        plt.fill_betweenx(y, x1, x2, color = 'm', alpha=0.5)

# PPTA
Violin = Read_Curve_Pro()
n = len(Violin)/2

for idx in np.arange(0, n):
    x1 = Violin[int(idx*2)]['x']
    y = Violin[int(idx*2)]['y']
    x2 = 2*x1[0] - x1
    y = 10**y/h2
    x1 = 10**x1
    x2 = 10**x2
    if idx == 0:
        plt.fill_betweenx(y, x1, x2, color = 'y', alpha=0.5, label = 'IPTA')
    else:
        plt.fill_betweenx(y, x1, x2, color = 'y', alpha=0.5)

# PPTA
Violin = Read_Curve_Pro(
    File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2306_17124_fig1a_PPTA.txt')
n = len(Violin)

for idx in np.arange(0, n):
    
    x1 = Violin[idx]['x']
    y = Violin[idx]['y']
    nx = len(x1)
    x2 = 2*x1[0] - x1

    y = 10**y/h2
    x1 = 10**x1
    x2 = 10**x2
    
    if idx == 0:
        plt.fill_betweenx(y, x1, x2, color = 'c', alpha=0.5, label = 'PPTA')
    else:
        plt.fill_betweenx(y, x1, x2, color = 'c', alpha=0.5)
    #plt
    #plt.plot(x2, y, '+k')



plt.legend(fontsize=FontSize,loc = 'upper left')
#plt.legend(fontsize=FontSize,loc = 'lower right')

plt.xscale('log')
plt.yscale('log')

# plt.title('Posterior for merger GW',fontsize=FontSize)
plt.xlabel('$f$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
#plt.xlim([1e-3, 1e10])
plt.ylim([1e-13, 1e-3])

plt.tight_layout()
plt.savefig(ResultFile, dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')
print('Plot saved to :')
print(ResultFile)

# plt.show()
