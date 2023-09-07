from src.merger_module import *

reload = 0
FileRoot = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/3_GW_Neff_fbh/3_GW_Neff_fbh_'
v = np.logspace(0, 8, 20)
LineWidth = 2
FontSize = 20

import matplotlib.pyplot as plt

def gw_posterior_model(theta):
    lf, lm, sbh = theta[5], theta[6], theta[7]
    fbh = 10**lf
    mc = 10**lm
    g = Get_dOmGW_dlnv(
        v = v, 
        mf_model = 0,
        sbh = sbh, 
        fbh = fbh, 
        mc = mc,
        show_status = 0,
        ncpu = 1,
        nm = 30,
        nz = 30,
        sbh_width = 5,
        zmax = 25)
    PyLab.SaySomething()
    return g

if reload:
    t1 = PyLab.TimeNow()
    r = PyLab.mcmc_derived_stat(model_function = gw_posterior_model, FileRoot = FileRoot, ncpu = 12, print_status=1)
    np.savez('3_merger_GW_posterior.npz', r = r)
    PyLab.Timer(t1)

r = np.load('3_merger_GW_posterior.npz')['r']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.fill_between(v, r[3,:], r[4,:],color = 'b',alpha=0.5,label='95 \%C.L.')
plt.plot(v, r[0,:], 'k', linewidth = LineWidth, label='Mean')
plt.xscale('log')
plt.yscale('log')

plt.legend(fontsize=FontSize,loc = 'lower right')
plt.xlabel('$\\nu$ [Hz]',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\Omega_{\mathrm{GW}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
# plt.ylim([1e0, 1e3])

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight', dpi = 1000)
# plt.savefig('/Users/cangtao/Desktop/tmp.eps',bbox_inches='tight')