
reload = 0
d1 = 0.25
d2 = 1/3
d3 = 0.45
ncpu = 12
map_nx = 500

LineWidth = 2
FontSize = 18

from src.src_1_Lite import *

kc = np.logspace(1, 20, 30)
f1 = lambda x: Find_Pmax(Sigma = 1, kc = x, fbh_max = 1, map_nx = map_nx, show_status = True, DeltaC = d1)
f2 = lambda x: Find_Pmax(Sigma = 1, kc = x, fbh_max = 1, map_nx = map_nx, show_status = True, DeltaC = d2)
f3 = lambda x: Find_Pmax(Sigma = 1, kc = x, fbh_max = 1, map_nx = map_nx, show_status = True, DeltaC = d3)

if reload:

    p1 = Parallel(n_jobs=ncpu)(delayed(f1)(x) for x in kc)
    p2 = Parallel(n_jobs=ncpu)(delayed(f2)(x) for x in kc)
    p3 = Parallel(n_jobs=ncpu)(delayed(f3)(x) for x in kc)
    np.savez('tmp.npz', p1 = p1, p2 = p2, p3 = p3)

r = np.load('tmp.npz')
p1 = r['p1']
p2 = r['p2']
p3 = r['p3']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.plot(kc, p3, '-k',linewidth=LineWidth, label = '$\delta_{\mathrm{c}} = 0.45$')
plt.plot(kc, p2, '-r',linewidth=LineWidth, label = '$\delta_{\mathrm{c}} = 1/3$')
plt.plot(kc, p1, '-b',linewidth=LineWidth, label = '$\delta_{\mathrm{c}} = 0.25$')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('$k\ [{\mathrm{Mpc^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$P_{\zeta}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.title('$P_{\zeta}$ amplitude for $f_{\mathrm{bh}} = 1$, $\Delta = 1$',fontsize=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper right')
# plt.ylim([-1.2,1.2])
plt.tight_layout()

# Save plot to a file
plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=500)
# plt.show()
