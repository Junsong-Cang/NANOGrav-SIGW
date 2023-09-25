from PyLab import *

reload = 0
Use_Log = 0
FileRoot_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
FileRoot_4 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_rate_posteriors.pdf'

LineWidth = 2
FontSize = 18
ncpu = 12
z = np.logspace(0, 2.2, 40) - 1
# params for best-fit
best_params = [
    [0.1, 10**-3.589, 0.54],
    [0.1, 10**-3.589, 0.01],
    [0.1, 10**-3.589, 1],
    [0.1, 10**-3.589, 2]]

# Initialising
Log30 = np.loadtxt(FileRoot_3 + 'Merger_Rate_0.txt')
Log31 = np.loadtxt(FileRoot_3 + 'Merger_Rate_1.txt')
Log40 = np.loadtxt(FileRoot_4 + 'Merger_Rate_0.txt')
Log41 = np.loadtxt(FileRoot_4 + 'Merger_Rate_1.txt')
nz = len(z)
N = len(best_params)

def kernel(theta, type):
    lf, lm, sbh = theta[5], theta[6], theta[7]
    # fbh = 10**lf
    # mc = 10**lm
    if type == 30:
        Tab = Log30
    elif type == 31:
        Tab = Log31
    elif type == 40:
        Tab = Log40
    else:
        Tab = Log41
    lf_vec = Tab[:,0]
    lm_vec = Tab[:,1]
    sbh_vec = Tab[:,2]
    dist = (lf - lf_vec)**2 + (lm - lm_vec)**2 + (sbh - sbh_vec)**2
    idx = np.argmin(dist)
    r = Tab[idx, 3 : 3+nz]
    if Use_Log:
        r = np.log10(r)
    return r

def model_30(theta):
    return kernel(theta, 30)

def model_31(theta):
    return kernel(theta, 31)

def model_40(theta):
    return kernel(theta, 40)

def model_41(theta):
    return kernel(theta, 41)

if reload:

    t1 = TimeNow()
    
    # Get posterior
    P30 = mcmc_derived_stat(model_function = model_30, FileRoot = FileRoot_3, ncpu = ncpu, print_status=1)
    P31 = mcmc_derived_stat(model_function = model_31, FileRoot = FileRoot_3, ncpu = ncpu, print_status=1)
    P40 = mcmc_derived_stat(model_function = model_40, FileRoot = FileRoot_4, ncpu = ncpu, print_status=1)
    P41 = mcmc_derived_stat(model_function = model_41, FileRoot = FileRoot_4, ncpu = ncpu, print_status=1)
    
    # Get data for the best-fit
    best_fit = np.zeros((N, nz))
    for idx in np.arange(0, N):
        params = [0, 1, 2, 3, 4] + best_params[idx]
        r = model_30(params)
        if Use_Log:
            best_fit[idx,:] = 10**r[:]
        else:
            best_fit[idx,:] = r[:]

    Timer(t1)
    np.savez('data/Merger_Rate_Posterior.npz', P30 = P30, P31 = P31, P40 = P40, P41 = P41, best_fit = best_fit)

R = np.load('data/Merger_Rate_Posterior.npz')
if Use_Log:
    R30 = 10**R['P30']
    R31 = 10**R['P31']
    R40 = 10**R['P40']
    R41 = 10**R['P41']
else:
    R30 = R['P30']
    R31 = R['P31']
    R40 = R['P40']
    R41 = R['P41']

best_fit = R['best_fit']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

Smooth = 4

plt.fill_between(1+z[::Smooth], R40[3,::Smooth], R40[4,::Smooth],color = 'g',alpha=0.4,label = 'No $S_2$')
plt.fill_between(1+z[::Smooth], R41[3,::Smooth], R41[4,::Smooth],color = 'b',alpha=0.6,label = 'With $S_2$')
plt.fill_between(1+z[::Smooth], R40[4,::Smooth], R30[4,::Smooth], color = 'grey',alpha=0.4, linestyle = 'dashed')

# plt.loglog(1+z[::Smooth], best_fit[0,::Smooth], 'k', linewidth = LineWidth)
plt.loglog(1+z, best_fit[0,:], 'k', linewidth = LineWidth)

plt.text(1e1, 3e9, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = 25,color='k')

plt.xscale('log')
plt.yscale('log')

#plt.title('Posterior for merger rate',fontsize=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Merger Rate $[{\mathrm{Gpc^{-3}}yr^{-1}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$R \ [{\mathrm{Gpc^{-3}yr^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([1, 100])
plt.ylim([1e6, 1e13])

plt.tight_layout()
plt.savefig(ResultFile)
print('Plot saved to :')
print(ResultFile)

x = np.log(1+z)
y = np.log(best_fit[0,:])
P = np.polyfit(x = x, y = y, deg = 1)
print(P)

plt.show()
