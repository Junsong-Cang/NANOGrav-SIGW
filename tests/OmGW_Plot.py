from src.merger_module import *
import matplotlib.pyplot as plt

reload = 1

fbh = 2E-2
mc = 1e9
sbh = 0.05

v = np.logspace(-11, 0, 40)
ncpu = 12
LineWidth = 2
FontSize = 18

if reload:
    t1 = PyLab.TimeNow()
    g = Get_dOmGW_dlnv(v = v, show_status = 1, ncpu = ncpu, sbh = sbh, fbh = fbh, mc = mc, nm = 200, nz = 200, mf_model = 0)
    PyLab.Timer(t1)
    np.savez('tmp.npz', g = g)

g = np.load('tmp.npz')['g']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

# fig.set_size_inches(8, 8)

plt.loglog(v, g, 'k', linewidth=LineWidth, label = 'Our codes')
plt.xlabel('$v$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('GW',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower left')

#plt.xlim([1e-6, 1e4])
#plt.ylim([1e-14,1e-5])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
plt.show()
