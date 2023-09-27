from PyLab import *

reload_posterior = 0
reload_samples = 0
reload_best_fit = 0

Use_Log = 0

FileRoot_3 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
FileRoot_4 = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/4_GW_Neff_fbh_01/4_GW_Neff_fbh_01_'
ResultFile='/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/LaTex/figs/Merger_rate_posteriors.pdf'

LineWidth = 2
FontSize = 16
ncpu = 12
z = np.logspace(0, 2.2, 40) - 1

best_params = [
    [0.1, 10**-3.589, 0.54],
    [0.1, 10**-3.589, 0.01],
    [0.1, 10**-3.589, 1],
    [0.1, 10**-3.589, 2]]

from src.merger import *
nz = len(z)

# -------- samples --------
# samples at z=0
if reload_samples:

    def Get_params(idx, sample):
        if sample == 3:
            File = FileRoot_3 + '.txt'
        elif sample == 4:
            File = FileRoot_4 + '.txt'
        Tab = np.loadtxt(File)
        chain = Tab[idx,:]
        fbh = 10**chain[7]
        mc = 10**chain[8]
        sbh = chain[9]
        r = fbh, mc, sbh
        return r

    def sample_kernel(idx, sample):
        fbh, mc, sbh = Get_params(idx = idx, sample = sample)
        R = Merger_Rate(
            fbh = fbh, 
            mc = mc, 
            sbh = sbh, 
            z = 0, 
            mf_model = 0, 
            sbh_width = 7, 
            nm = 400,
            Use_S2 = 0, 
            S1_method = 0,
            Use_interp = 1,
            S_Tab_Len = 500,
            show_status = 0)
        PyLab.SaySomething('tmp.txt')
        lf, lm = np.log10(fbh), np.log10(mc)
        r = np.array([lf, lm, sbh, R])
        return r
    
    n3 = len(np.loadtxt(FileRoot_3 + '.txt'))
    n4 = len(np.loadtxt(FileRoot_4 + '.txt'))
    t1 = PyLab.TimeNow()
    r3 = Parallel(n_jobs = ncpu)(delayed(sample_kernel)(idx = idx, sample = 3) for idx in np.arange(0, n3))
    r4 = Parallel(n_jobs = ncpu)(delayed(sample_kernel)(idx = idx, sample = 4) for idx in np.arange(0, n4))
    np.savez('data/Merger_Rate_Samples.npz', r3 = r3, r4 = r4)
    PyLab.Timer(t1)

# -------- posterior --------

if reload_posterior:

    Template = np.load('data/Merger_Rate_Samples.npz')
    Tab_3 = Template['r3']
    Tab_4 = Template['r4']

    def kernel(theta, sample, Use_S2):
        lf, lm, sbh = theta[5], theta[6], theta[7]
        if sample == 3:
            Tab = Tab_3
        else:
            Tab = Tab_4
        lf_vec = Tab[:,0]
        lm_vec = Tab[:,1]
        sbh_vec = Tab[:,2]
        R0_vec = Tab[:,3]
        
        dist = (lf - lf_vec)**2 + (lm - lm_vec)**2 + (sbh - sbh_vec)**2
        idx = np.argmin(dist)
        R0 = R0_vec[idx]
        if Use_S2:
            # Samples has no S2, add it now
            t0 = 4.3523939036782445e+17
            S2 = S_Factor(fbh = 10**lf, sbh=0,mc=0,M=0,m_ave=0,m2_ave=0,mf_model=0,t = t0, Use_S2=1,S1_method=2)
            R0 = R0*S2
        
        # now do scaling
        r = np.zeros(nz)
        for idx in np.arange(0, nz):
            z_ = z[idx]
            S = Merger_Rate_Scale_Factor(z = z_, fbh = 10**lf, Use_S2 = Use_S2)
            r[idx] = R0*S
        if Use_Log:
            r = np.log10(r)
        return r
    
    def model_30(theta):
        return kernel(theta, 3, 0)
    def model_31(theta):
        return kernel(theta, 3, 1)
    def model_40(theta):
        return kernel(theta, 4, 0)
    def model_41(theta):
        return kernel(theta, 4, 1)
    
    # Get posterior
    t1 = TimeNow()
    P30 = mcmc_derived_stat(model_function = model_30, FileRoot = FileRoot_3, ncpu = ncpu, print_status=1)
    P31 = mcmc_derived_stat(model_function = model_31, FileRoot = FileRoot_3, ncpu = ncpu, print_status=1)
    P40 = mcmc_derived_stat(model_function = model_40, FileRoot = FileRoot_4, ncpu = ncpu, print_status=1)
    P41 = mcmc_derived_stat(model_function = model_41, FileRoot = FileRoot_4, ncpu = ncpu, print_status=1)
    Timer(t1)
    np.savez('data/Merger_Rate_Posterior.npz', P30 = P30, P31 = P31, P40 = P40, P41 = P41)

# -------- best-fit --------

if reload_best_fit:
    t1 = TimeNow()
    N = len(best_params)
    best_fit = np.zeros((N, nz))
    for idx in np.arange(0, N):
        print('status : ', idx/N)
        fbh, mc, sbh = best_params[idx]
        R = Merger_Rate(
            fbh = fbh,
            mc = mc, 
            sbh = sbh, 
            z = z, 
            mf_model = 0, 
            sbh_width = 8, 
            nm = 500,
            Use_S2 = 0,
            S1_method = 0,
            Use_interp = 1,
            S_Tab_Len = 1000,
            show_status = 0)
        best_fit[idx,:] = R[:]
    Timer(t1)
    np.savez('data/Merger_Rate_Best_Fit.npz', best_fit = best_fit)

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

best_fit = np.load('data/Merger_Rate_Best_Fit.npz')['best_fit']

# -------- plot --------
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

Smooth = 4

plt.fill_between(1+z[::Smooth], R40[3,::Smooth], R40[4,::Smooth],color = 'g',alpha=0.4,label = 'No $S_2$')
plt.fill_between(1+z[::Smooth], R41[3,::Smooth], R41[4,::Smooth],color = 'b',alpha=0.6,label = 'With $S_2$')
plt.fill_between(1+z[::Smooth], R40[4,::Smooth], R30[4,::Smooth], color = 'grey',alpha=0.4, linestyle = 'dashed')

plt.loglog(1+z, best_fit[0,:], 'k', linewidth = LineWidth)

plt.text(1e1, 3e9, "$f_{\mathrm{bh}} < 1$", size=FontSize/1.3, rotation = 25,color='k')

plt.xscale('log')
plt.yscale('log')

#plt.title('Posterior for merger rate',fontsize=FontSize)
plt.legend(fontsize=FontSize/1.2,loc = 'upper left')
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
