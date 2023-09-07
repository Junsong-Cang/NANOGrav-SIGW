from src.merger_module import *

reload = 0
v = np.logspace(0, 8, 40)

LineWidth = 2
FontSize = 18

fbh = np.array([1e-2, 1e-3, 1e-1])
mc = np.array([1e-3, 1e-4, 1e-2])
sbh = np.array([0.3, 0.1, 1])

import matplotlib.pyplot as plt
def model(f, m, s):
    r = Get_dOmGW_dlnv(fbh = f, mc = m, sbh = s, v = v, mf_model = 0, sbh_width = 6, nm = 200, nz = 200, ncpu = 12)
    return r

if reload:
    
    r = np.zeros((9, len(v)))

    idx = 0
    t1 = PyLab.TimeNow()
    r[0, :] = model(f = fbh[0], m = mc[0], s = sbh[0])

    for fid in np.arange(1,3):
        for mid in np.arange(1,3):
            for sid in np.arange(1,3):
                print(idx)
                idx = idx + 1
                f_ = fbh[fid]
                m_ = mc[mid]
                s_ = sbh[sid]
                r[idx, :] = model(f = f_, m = m_, s = s_)
    np.savez('tmp.npz', r = r)
    PyLab.Timer(t1)

r = np.load('tmp.npz')['r']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(v, r[0,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[1,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[2,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[3,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[4,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[5,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[6,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[7,:], 'k', linewidth=LineWidth)
plt.loglog(v, r[8,:], 'k', linewidth=LineWidth)
plt.xlabel('$\\nu$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('GW',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

#plt.xlim([-1,7])
#plt.ylim([-1.2,1.2])
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
