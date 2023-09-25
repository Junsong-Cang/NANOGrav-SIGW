from src.merger import *

reload = 1

FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_rate_posteriors.pdf'

LineWidth = 2
FontSize = 18
nm = 100
S_Tab_Len = 200
ncpu = 12
z = np.logspace(0, 2.2, 40) - 1
LogFile_0 = FileRoot + 'Merger_Rate_0.txt'
LogFile_1 = FileRoot + 'Merger_Rate_1.txt'

# params for best-fit
best_params = [
    [0.1, 10**-3.589, 0.54],
    [0.1, 10**-3.589, 0.01],
    [0.1, 10**-3.589, 1],
    [0.1, 10**-3.589, 2]]

import matplotlib.pyplot as plt
import os

nz = len(z)
N = len(best_params)

def model_0(theta):
    Use_S2 = 0
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    r = Merger_Rate(
        fbh = fbh, 
        mc = mc, 
        sbh = sbh, 
        z = z,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm, 
        Use_S2 = Use_S2, 
        S1_method = S1_method,
        Use_interp = 1,
        S_Tab_Len = S_Tab_Len,
        show_status = 0)
    
    # Print results to a file
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    
    for idx in np.arange(0, nz):
        Str = Str + '  {0:.5E}'.format(r[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    PyLab.SaySomething()
    return r

def model_1(theta):
    Use_S2 = 1
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    r = Merger_Rate(
        fbh = fbh, 
        mc = mc, 
        sbh = sbh, 
        z = z,
        mf_model = 0, 
        sbh_width = 6, 
        nm = nm, 
        Use_S2 = Use_S2, 
        S1_method = S1_method,
        Use_interp = 1,
        S_Tab_Len = S_Tab_Len,
        show_status = 0)
    
    # Print results to a file
    if Use_S2:
        File = LogFile_1
    else:
        File = LogFile_0
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    
    for idx in np.arange(0, nz):
        Str = Str + '  {0:.5E}'.format(r[idx])
    F = open(File, 'a')
    print(Str, file = F)
    F.close()

    PyLab.SaySomething()
    return r

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
    
    # Get posterior
    # Get posterior
    P0 = PyLab.mcmc_derived_stat(model_function = model_0, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    P1 = PyLab.mcmc_derived_stat(model_function = model_1, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    
    # Get data for the best-fit
    best_fit = np.zeros((N, nz))
    for idx in np.arange(0, N):
        params = [0, 1, 2, 3, 4] + best_params[idx]
        r = model_0(params)
        best_fit[idx,:] = r[:]

    np.savez('data/3_merger_rate_posterior.npz', P0 = P0, P1 = P1, best_fit = best_fit)
    PyLab.Timer(t1)

R3 = np.load('data/3_merger_rate_posterior.npz')
p3 = R3['Posterior']
best_fit = R3['best_fit']

R4 = np.load('data/4_merger_rate_posterior.npz')
p4 = R4['Posterior']
p4b = R4['Posterior_2']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

Smooth = 4

plt.fill_between(1+z[::Smooth], p4[3,::Smooth], p4[4,::Smooth],color = 'g',alpha=0.4,label = 'No $S_2$')
plt.fill_between(1+z[::Smooth], p4b[3,::Smooth], p4b[4,::Smooth],color = 'b',alpha=0.6,label = 'With $S_2$')
plt.fill_between(1+z[::Smooth], p3[4,::Smooth], p4[4,::Smooth], color = 'grey',alpha=0.4, linestyle = 'dashed')

# plt.loglog(1+z[::Smooth], best_fit[0,::Smooth], 'k', linewidth = LineWidth)
plt.loglog(1+z, best_fit[0,:], 'k', linewidth = LineWidth)

plt.text(1e1, 3e9, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = 25,color='k')

plt.xscale('log')
plt.yscale('log')

#plt.title('Posterior for merger rate',fontsize=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[{\mathrm{Gpc^{-3}}yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$R \ [{\mathrm{Gpc^{-3}}yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1, 100])
plt.ylim([1e6, 1e13])

plt.tight_layout()
plt.savefig(ResultFile)
print('Plot saved to :')
print(ResultFile)
