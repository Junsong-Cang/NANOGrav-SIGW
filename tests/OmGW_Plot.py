from src.merger import *
import matplotlib.pyplot as plt

reload = 0

fbh = 1
mc = 1e5
sbh = 0.05

v = np.logspace(-11, 0, 400)
ncpu = 12
LineWidth = 2
FontSize = 18

if reload:
    t1 = PyLab.TimeNow()
    g = Get_dOmGW_dlnv(v = v, show_status = 0, ncpu = ncpu, sbh = sbh, fbh = fbh, mc = mc, nm = 100, nz = 70, mf_model = 2)
    PyLab.Timer(t1)
    np.savez('tmp.npz', g = g)

g = np.load('tmp.npz')['g']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

# fig.set_size_inches(8, 8)

plt.loglog(v, g*0.6766**2, 'k', linewidth=LineWidth, label = 'Our codes')
plt.xlabel('$v$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('GW',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower left')

plt.xlim([1e-11, 1e0])
plt.ylim([1e-11,1e-5])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
plt.show()
