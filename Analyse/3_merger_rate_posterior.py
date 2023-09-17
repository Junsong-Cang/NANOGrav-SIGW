from src.merger import *

reload = 0

FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_rate_posteriors.pdf'

LineWidth = 2
FontSize = 18
nm = 100
S_Tab_Len = 200
ncpu = 12
z = np.logspace(0, 2.2, 40) - 1
LogFile = FileRoot + 'Merger_Rate.txt'

import matplotlib.pyplot as plt
import os

Tab = np.loadtxt(LogFile)
def model(theta):
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
    Str = '{0:.5E}'.format(lf) + '  {0:.5E}'.format(lm) + '  {0:.5E}'.format(sbh)
    nz = len(z)
    
    for idx in np.arange(0, nz):
        Str = Str + '  {0:.5E}'.format(r[idx])
    F = open(LogFile, 'a')
    print(Str, file = F)
    F.close()

    PyLab.SaySomething()
    return r

def model_2(theta):
    # By the time this is called model_1 should be complete so we can just find results in previous table
    Use_S2 = 1
    S1_method = 0
    lf, lm, sbh = theta[5], theta[6], theta[7]
    # lf = min(lf,-1)
    if lf > -1:
        lf = -1
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
        PyLab.SaySomething()
    else:
        nz = len(z)
        lf_axis = Tab[:,0]
        lm_axis = Tab[:,1]
        sbh_axis = Tab[:,2]
        dist = (lf - lf_axis)**2 + (lm - lm_axis)**2 + (sbh - sbh_axis)**2
        idx = np.argmin(dist)
        r = Tab[idx, 3:nz+3]
    return r

if reload:
    # os.remove(LogFile)
    t1 = PyLab.TimeNow()
    r = PyLab.mcmc_derived_stat(model_function = model, FileRoot = FileRoot, ncpu = ncpu, print_status=1)
    r2 = PyLab.mcmc_derived_stat(model_function = model_2, FileRoot = FileRoot, ncpu = ncpu, print_status = 1)
    # Get a plot for the best-fit
    lf = -1
    lm = -3
    sbh = 0.3
    params = [0,1,2,3,4,lf, lm, sbh]
    r0 = model(params)
    np.savez('data/3_merger_rate_posterior.npz', r = r, r0 = r0, r2 = r2)
    PyLab.Timer(t1)

R = np.load('data/3_merger_rate_posterior.npz')
r = R['r']
r0 = R['r0']
r2 = R['r2']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

# plt.fill_between(1+z, r2[3,:], r2[4,:],color = 'g',alpha=0.9,label = 'Reliable')
plt.fill_between(1+z, r[3,:], r[4,:],color = 'b',alpha=0.5,label = '$2 \sigma$')
plt.loglog(1+z, r0, 'k', linewidth = LineWidth, label = 'best-fit')

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
plt.ylim([1e5, 1e12])

plt.tight_layout()
plt.savefig(ResultFile)
print('Plot saved to :')
print(ResultFile)
