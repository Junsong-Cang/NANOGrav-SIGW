from src.merger import *
import matplotlib.pyplot as plt

reload = 0

fbh = 1e-1
mc = 2.6e-4
sbh = 0.54

v = np.logspace(-4, 10, 50)
ncpu = 12
LineWidth = 2
FontSize = 18

if reload:
    t1 = PyLab.TimeNow()
    g = Get_dOmGW_dlnv(fbh = fbh, mc = mc, sbh = sbh, v = v, mf_model = 0, nm = 50, nz = 1000, zmax = 1000, Fast = 1, Use_S2 = 0)
    PyLab.Timer(t1)
    np.savez('tmp.npz', g = g)

g = np.load('tmp.npz')['g']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

# fig.set_size_inches(8, 8)

GW_1d = np.load('/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/Analyse/data/4_merger_GW.npz')
v4 = GW_1d['v']
r4 = GW_1d['r'][0,:]

plt.loglog(v, g, 'k', linewidth=LineWidth, label = '$z_{max} = 1000$')
plt.loglog(v4, r4, 'r', linewidth=LineWidth, label = '$z_{max} = \infty$')

plt.xlabel('$v$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('GW',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'lower left')

plt.xlim([1e-1, 1e9])
plt.ylim([1e-17,1e-6])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
tmp = np.interp(x = -1, xp = np.log(v), fp = g)
print(tmp)
if tmp<4e-16:
    print('Null ')
# plt.show()
